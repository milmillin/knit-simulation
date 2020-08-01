#include "Simulator.h"

#include "./SimulatorParams.h"
#include "./Helper.h"
#include "../file_format/yarnRepr.h"
#include "../easy_profiler_stub.h"
#include "./threading/threading.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <algorithm>

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
static const float MASS_CONTRIBUTION[4][4] =
{
  {1.f / 420, -47.f / 1680, -1.f / 56, 1.f / 560},
  {-47.f / 1680, 17.f / 42, 307.f / 1680, -1.f / 56},
  {-1.f / 56, 307.f / 1680, 17.f / 42, -47.f / 1680},
  {1.f / 560, -1.f / 56, -47.f / 1680, 1.f / 420}
};

void Simulator::constructMassMatrix() {
  M = Eigen::SparseMatrix<float>(3ll * m, 3ll * m);

  // Iterate though each segment
  int N = m - 3;

  float mUnit = params.m;
  for (int i = 0; i < N; i++) {
    // Iterate through combination of control points to find mass contribution.
    for (int j = 0; j < 4; j++) {
      for (int k = j; k < 4; k++) {
        // M[(i+j)*3][(i+k)*3] = mUnit * l[i] * integrate(bj(s), bk(s), (s, 0, 1))
        float contribution = mUnit * catmullRomLength[i] * MASS_CONTRIBUTION[j][k];
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

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
  solver.compute(M);
  if (solver.info() != Eigen::Success) {
    log() << "Decomposition Failed" << std::endl;
    return;
  }

  Eigen::SparseMatrix<float> I(M.rows(), M.cols());
  I.setIdentity();

  invM = solver.solve(I);
  if (solver.info() != Eigen::Success) {
    log() << "Solve Failed" << std::endl;
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
	F += Eigen::Vector3f(0, -params.gravity, 0).replicate(m, 1);
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
  IFDEBUG log() << "Calculating Gradient" << std::endl;
  int N = m - 3;

  /*
  // Energy
  if (cancelled()) return;
  IFDEBUG log() << "- Length Energy" << std::endl;
  EASY_BLOCK("Length Energy");
  calculate(&Simulator::calculateLengthEnergy, 0, N);
  EASY_END_BLOCK;
  */


  if (cancelled()) return;
  IFDEBUG log() << "- Bending Energy" << std::endl;
  EASY_BLOCK("Bending Energy");
  calculate(&Simulator::calculateBendingEnergy, 0, N);
  EASY_END_BLOCK;

  /*
  if (cancelled()) return;
  IFDEBUG log() << "- Collision Energy" << std::endl;
  applyContactForce(cancelled);

  // Damping
  if (cancelled()) return;
  IFDEBUG log() << "- Global Damping" << std::endl;
  EASY_BLOCK("Global Damping");
  calculate(&Simulator::calculateGlobalDamping, 0, N);
  EASY_END_BLOCK;
  */

  /*
  if (cancelled()) return;
  IFDEBUG log() << "- Gravity" << std::endl;
  EASY_BLOCK("Gravity");
  calculateGravity();
  EASY_END_BLOCK;
  */

  // unconstrained step
  const float& h = params.h;
  Q += h * dQ + (h * h) * (invM * F);
}

void Simulator::setUpConstraints() {
  int N = m - 3;

  for (int i = 0; i < N; i++) {
    addCatmullRomLengthConstraint(i);
  }

  addSegmentLengthConstraint(0);
  addSegmentLengthConstraint(m - 2);

  addPinConstraint(0, pointAt(Q, 0));
  addPinConstraint(m - 1, pointAt(Q, m - 1));
}

};  // namespace simulator
