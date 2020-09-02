#include "Simulator.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "spdlog/spdlog.h"
#include "easy_profiler_stub.h"

#include "file_format/yarnRepr.h"
#include "./SimulatorParams.h"
#include "./Helper.h"
#include "./threading/threading.h"

namespace simulator {

//////////////////////////////////////////////
//
// Mass Matrix
//

// Mass contribution matrix of Catmull-Rom curve with tightness 0.5.
// b0(s) = (1/2)(-s + 2s^2 - s^3)
// b1(s) = (1/2)(2 - 5s^2 + 3s^3)
// b2(s) = (1/2)(s + 4s^2 - 3s^3)
// b3(s) = (1/2)(-s^2 + s^3)
//
// MASS_CONTRIBUTION[i][j] = integrate(bi(s)*bj(s), (s, 0, 1))
// TODO: check correctness
static const double MASS_CONTRIBUTION[4][4] =
{
  {1. / 420, -47. / 1680, -1. / 56, 1. / 560},
  {-47. / 1680, 17. / 42, 307. / 1680, -1. / 56},
  {-1. / 56, 307. / 1680, 17. / 42, -47. / 1680},
  {1. / 560, -1. / 56, -47. / 1680, 1. / 420}
};

void Simulator::constructMassMatrix() {
  M = Eigen::SparseMatrix<double>(3ll * m, 3ll * m);

  // Iterate though each segment
  int N = m - 3;

  double mUnit = params.m;
  for (int i = 0; i < N; i++) {
    // Iterate through combination of control points to find mass contribution.
    for (int j = 0; j < 4; j++) {
      for (int k = j; k < 4; k++) {
        // M[(i+j)*3][(i+k)*3] = mUnit * l[i] * integrate(bj(s), bk(s), (s, 0, 1))
        double contribution = mUnit * catmullRomLength[i] * MASS_CONTRIBUTION[j][k];
        M.coeffRef((i + j) * 3, (i + k) * 3) += contribution;
        M.coeffRef((i + j) * 3 + 1, (i + k) * 3 + 1) += contribution;
        M.coeffRef((i + j) * 3 + 2, (i + k) * 3 + 2) += contribution;

        // transpose
        if (j != k) {
          M.coeffRef((i + k) * 3, (i + j) * 3) += contribution;
          M.coeffRef((i + k) * 3 + 1, (i + j) * 3 + 1) += contribution;
          M.coeffRef((i + k) * 3 + 2, (i + j) * 3 + 2) += contribution;
        }
      }
    }
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  if (solver.info() != Eigen::Success) {
    SPDLOG_ERROR("Decomposition Failed");
    return;
  }

  Eigen::SparseMatrix<double> I(M.rows(), M.cols());
  I.setIdentity();

  invM = solver.solve(I);
  if (solver.info() != Eigen::Success) {
    SPDLOG_ERROR("Solve Failed");
    return;
  }
}

//////////////////////////////////////////////
//
// Energy Gradient
//

void Simulator::calculate(void (Simulator::* func)(int), int start, int end) {
  auto calculateTask = [this, &func](int thread_id, int start_index, int end_index) {
    for (int i = start_index; i < end_index; i++) {
      (this->*func)(i);
    }
  };

  threading::runSequentialJob(thread_pool, calculateTask, start, end);
}

void Simulator::calculateGravity() {
	F += Eigen::Vector3d(0, -params.gravity, 0).replicate(m, 1);
}

//////////////////////////////////////////////
//
// Simulator implementation
//

Simulator::Simulator(file_format::YarnRepr yarns, SimulatorParams params_) :
  BaseSimulator(yarns, params_)
{
  initialize();
}

void Simulator::stepImpl(const StateGetter& cancelled) {
  F.setZero();

  // update gradE, gradD, f
  if (params.debug)
    SPDLOG_INFO("Calculating Gradient");
  int N = m - 3;

  // Energy
  if (params.kLen != 0) {
    if (cancelled()) return;
    if (params.debug)
      SPDLOG_INFO("- Length Energy");
    EASY_BLOCK("Length Energy");
    calculate(&Simulator::calculateLengthEnergy, 0, N);
    EASY_END_BLOCK;
  }

  if (params.kBend != 0) {
    if (cancelled()) return;
    if (params.debug)
      SPDLOG_INFO("- Bending Energy");
    EASY_BLOCK("Bending Energy");
    calculate(&Simulator::calculateBendingEnergy, 0, N);
    EASY_END_BLOCK;
  }

  if (cancelled()) return;
  if (params.debug)
    SPDLOG_INFO("- Collision Energy");
  applyContactForce(cancelled);

  // Damping
  if (params.kGlobal != 0) {
    if (cancelled()) return;
    if (params.debug)
      SPDLOG_INFO("- Global Damping");
    EASY_BLOCK("Global Damping");
    calculate(&Simulator::calculateGlobalDamping, 0, N);
    EASY_END_BLOCK;
  }

  if (cancelled()) return;
  if (params.debug)
    SPDLOG_INFO("- Gravity");
  EASY_BLOCK("Gravity");
  calculateGravity();
  EASY_END_BLOCK;

  // unconstrained step
  const double& h = params.h;
  Q += h * dQ + (h * h) * (invM * F);
}

void Simulator::setUpConstraints() {
  int N = m - 3;

  for (int i = 0; i < N; i++) {
    addCatmullRomLengthConstraint(i);
  }

  addSegmentLengthConstraint(0);
  addSegmentLengthConstraint(m - 2);

  std::vector<int> pin = {
    159, 171, 183, 195
  };

  for (int p : pin) {
    addPinConstraint(p, pointAt(Q, p));
  }

  //addPinConstraint(0, pointAt(Q, 0));
  //addPinConstraint(m - 1, pointAt(Q, m - 1));
}

};  // namespace simulator
