#include "DiscreteSimulator.h"

#include <iostream>
#include <functional>

#include <Eigen/SparseLU>
#include "spdlog/spdlog.h"
#include "easy_profiler_stub.h"

#include "./threading/threading.h"
#include "./macros.h"
#include "./Helper.h"

namespace simulator {

static constexpr double pi = 3.14159265358979323846;

DiscreteSimulator::DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params) :
  BaseSimulator(yarns, params)
{
  initialize();

  // Initialize bending force calculation
  initBendingForceMetadata();
}

void DiscreteSimulator::stepImpl(const StateGetter& cancelled) {
  EASY_FUNCTION();

  // Calculate acceleration
  F.setZero();

  if (params.gravity >= 1e-6) applyGravity();
  checkNaN(F);
  if (params.kContact >= 1e-6) applyContactForce(cancelled);
  checkNaN(F);
  if (params.kBend >= 1e-6) applyBendingForce();
  checkNaN(F);
  if (params.kLen >= 1e-6) applyLengthSpringForce();
  checkNaN(F);

  Eigen::MatrixXd ddQ = F; // invM * F

  // Calculate velocity
  dQ += ddQ * params.h;
  if (params.enableGround) applyGroundVelocityFilter();
  checkNaN(dQ);
  if (params.kGlobalDamping >= 1e-6) applyGlobalDamping();
  checkNaN(dQ);

  // Calculate position
  Q += dQ * params.h;
  checkNaN(Q);
}

void DiscreteSimulator::postStep(const StateGetter& cancelled) {
  if (params.kBend >= 1e-6) updateBendingForceMetadata();
}

void DiscreteSimulator::applyGravity() {
  EASY_FUNCTION();

  F += Eigen::Vector3d(0, -params.gravity, 0).replicate(nControlPoints, 1);
}

void DiscreteSimulator::applyGroundVelocityFilter() {
  EASY_FUNCTION();

  for (const auto& yarn : yarns.yarns) {
    threading::runParallelFor(thread_pool, yarn.begin, yarn.end, [this](int thread_id, size_t i) {
      if (coordAt(Q, i, 1) < params.groundHeight && coordAt(dQ, i, 1) < 0) {
        coordAt(dQ, i, 0) *= params.groundFriction;
        coordAt(dQ, i, 1) = 0;
        coordAt(dQ, i, 2) *= params.groundFriction;
        coordAt(Q, i, 1) = params.groundHeight;
      }
      });
  }
}

// CHECK
void DiscreteSimulator::curvatureBinormalTask(int thread_id, size_t i) {
  Eigen::Vector3d a = 2 * e.row(i - 1).cross(e.row(i));
  double b = segmentLength[i - 1] * segmentLength[i];
  double c = e.row(i - 1).dot(e.row(i));
  curvatureBinormal.row(i) = a.transpose() / (b + c);
}

// CHECK
Eigen::Vector2d DiscreteSimulator::omega(size_t i, size_t j) {
  Eigen::Vector3d kb = curvatureBinormal.row(i);
  Eigen::Vector2d result;
  result(0) = kb.dot(m2.row(j));
  result(1) = -kb.dot(m1.row(j));
  return result;
}

void DiscreteSimulator::initBendingForceMetadata() {
  EASY_FUNCTION();

  // No bending energy
  if (nControlPoints <= 2) {
    return;
  }

  // Allocate memory
  e.resize(nControlPoints - 1, Eigen::NoChange);
  m1.resize(nControlPoints - 1, Eigen::NoChange);
  m2.resize(nControlPoints - 1, Eigen::NoChange);
  restOmega.resize(nControlPoints - 1, Eigen::NoChange);
  restOmega_1.resize(nControlPoints - 1, Eigen::NoChange);
  curvatureBinormal.resize(nControlPoints - 1, Eigen::NoChange);
  gradCurvatureBinormal.resize(nControlPoints - 1,
    std::vector<Eigen::Matrix3d>(3, Eigen::Matrix3d::Zero()));

  e.setZero();
  m1.setZero();
  m2.setZero();
  restOmega.setZero();
  restOmega_1.setZero();
  curvatureBinormal.setZero();

  // For each yarn
  for (const auto& yarn : yarns.yarns) {
    // Initialize tangent
    for (size_t i = yarn.begin; i < yarn.end - 1; i++) {
      e.row(i) = pointAt(Q, i + 1) - pointAt(Q, i);
    }

    // Initialize direction 1 with arbitrary vector that's normal to u.row(0)
    auto firstPoint = pointAt(Q, yarn.begin);
    if (std::abs(firstPoint(0) + firstPoint(1)) < std::abs(firstPoint(2))) {
      m1.row(yarn.begin) = Eigen::Vector3d(1, 0, 0).cross(e.row(yarn.begin));
    }
    else {
      m1.row(yarn.begin) = Eigen::Vector3d(0, 0, 1).cross(e.row(yarn.begin));
    }
    m1.row(yarn.begin).normalize();

    // Initialize direction 2
    m2.row(yarn.begin) = e.row(yarn.begin).cross(m1.row(yarn.begin)).normalized();

    // Fill in all frames
    for (size_t i = yarn.begin + 1; i < yarn.end - 1; i++) {
      m1.row(i) = parallelTransport(m1.row(i - 1), e.row(i - 1), e.row(i)).normalized();
      m2.row(i) = e.row(i).cross(m1.row(i)).normalized();
    }

    // FIXME: Need to initialize curvatureBinormal before calling omega?

    // Calculate rest omega
    for (size_t i = yarn.begin + 1; i < yarn.end - 1; i++) {
      restOmega.row(i) = omega(i, i).transpose();
      restOmega_1.row(i) = omega(i, i - 1).transpose();
    }
  }

  // Initialize theta
  u = m1;
  v = m2;
  theta = std::vector<double>(nControlPoints - 1, 0.0);
  thetaHat = std::vector<double>(nControlPoints - 1, 0.0);
  thetaHatOffset = std::vector<int>(nControlPoints - 1, 0);
}

static inline Eigen::Matrix3d crossMatrix(Eigen::Vector3d e) {
  Eigen::Matrix3d result;
  result <<
    0, -e(2), e(1),
    e(2), 0, -e(0),
    -e(1), e(0), 0;
  return result;
}

void DiscreteSimulator::gradCurvatureBinormalTask(int thread_id, size_t i) {
  Eigen::Matrix3d a = 2 * crossMatrix(e.row(i)) + 2 * crossMatrix(e.row(i - 1ull));
  Eigen::RowVector3d b = (e.row(i) - e.row(i - 1ull));
  double c = segmentLength[i - 1ull] * segmentLength[i] + e.row(i - 1ull).dot(e.row(i));
  // FIXME: assert(i - 1 >= 0 && i + 1 < nControlPoints - 1)
  for (size_t k = i - 1; k <= i + 1; k++) {
    gradCurvatureBinormal[i][k - (i - 1ull)] =
      -(a + curvatureBinormal.row(k) * b) / c;
  }
}

// CHECK
static double newThetaHat(RowMatrixX3d& e, RowMatrixX3d& u, size_t i) {
  Eigen::Vector3d newU = parallelTransport(u.row(i - 1ull), e.row(i - 1ull), e.row(i));
  return std::atan2(newU.cross(u.row(i)).norm(), newU.dot(u.row(i)));
}

void DiscreteSimulator::bendingForceTask(int thread_id, size_t i) {
  Eigen::Vector3d force = Eigen::Vector3d::Zero();

  // FIXME: assert(i - 1 >= 0 && i + 1 < nControlPoints - 1)
  for (size_t k = i - 1; k <= i + 1; k++) {
    double l = segmentLength[k - 1ull] + segmentLength[k];
    Eigen::Vector3d kb = curvatureBinormal.row(k);
    for (int j = k - 1; j <= k; j++) {
      Eigen::MatrixXd coeff(2, 3);
      coeff.row(0) = m2.row(j);
      coeff.row(1) = -m1.row(j);
      auto omegaBar = (j == k) ?
        restOmega.row(k).transpose() :
        restOmega_1.row(k).transpose();
      force -= (1.0 / l)
        * (coeff * gradCurvatureBinormal[i][k - (i - 1ull)]).transpose()
        * (omega(k, j) - omegaBar);
    }
  }
  pointAt(F, i) += params.kBend * force;
}

// CHECK
void DiscreteSimulator::twistingForceTask(int thread_id, size_t i) {
  Eigen::Vector3d dTheta_dQi_1 = curvatureBinormal.row(i) / 2 / segmentLength[i - 1ull];
  Eigen::Vector3d dTheta_dQi = curvatureBinormal.row(i) / 2 / segmentLength[i];
  double coeff = 2 * (theta[i] - theta[i + 1ull] - thetaHat[i]) / (segmentLength[i - 1ull] + segmentLength[i]);
  pointAt(F, i) += params.kTwist * coeff * (dTheta_dQi_1 - dTheta_dQi);
}

void DiscreteSimulator::applyBendingForce() {
  EASY_FUNCTION();
  using namespace std::placeholders;

  {
    EASY_BLOCK("curvatureBinormalTask");
    for (const auto& yarn : yarns.yarns) {
      auto task = std::bind(&DiscreteSimulator::curvatureBinormalTask, this, _1, _2);
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, task);
    }
  }

  {
    EASY_BLOCK("gradCurvatureBinormalTask");
    for (const auto& yarn : yarns.yarns) {
      auto task = std::bind(&DiscreteSimulator::gradCurvatureBinormalTask, this, _1, _2);
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, task);
    }
  }

  {
    EASY_BLOCK("bendingForceTask");
    for (const auto& yarn : yarns.yarns) {
      auto task = std::bind(&DiscreteSimulator::bendingForceTask, this, _1, _2);
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, task);
    }
  }

  {
    EASY_BLOCK("twistingForceTask");
    for (const auto& yarn : yarns.yarns) {
      auto task = std::bind(&DiscreteSimulator::twistingForceTask, this, _1, _2);
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, task);
    }
  }
}

void DiscreteSimulator::updateBendingForceMetadata() {
  EASY_FUNCTION();
  if (nControlPoints <= 2) return;

  {
    EASY_BLOCK("Update frame");
    for (const auto& yarn : yarns.yarns) {
      // FIXME: Start from 0?
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, [this](int thread_id, size_t i) {
          Eigen::Vector3d newE = pointAt(Q, i + 1) - pointAt(Q, i);
          u.row(i) = parallelTransport(u.row(i), e.row(i), newE).normalized();
          v.row(i) = newE.cross(u.row(i)).normalized();
          e.row(i) = newE.transpose();
        });
    }
  }

  {
    EASY_BLOCK("Solve theta");
    bool shouldContinue = true;
    std::vector<double> thetaUpdate(nControlPoints - 1);
    for (int iter = 0; shouldContinue && iter < params.materialFrameMaxUpdate; iter++) {
      shouldContinue = false;
      std::fill(thetaUpdate.begin(), thetaUpdate.end(), 0);

      for (const auto& yarn : yarns.yarns) {
        // NOTE: We are using the stress-free boundary condition, so theta[0] and theta[end] are always 0.
        threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 2, [this, &thetaUpdate](int thread_id, size_t i) {
          double li = segmentLength[i] + segmentLength[i - 1];
          thetaUpdate[i] = params.kTwist * 2 * (theta[i] - theta[i - 1] - thetaHat[i]) / li;
          double li_1 = segmentLength[i] + segmentLength[i + 1];
          thetaUpdate[i] -= params.kTwist * 2 * (theta[i + 1] - theta[i] - thetaHat[i + 1]) / li_1;
          });
      }

      double maxUpdate = 0;
      for (const auto& yarn : yarns.yarns) {
        for (size_t i = yarn.begin; i < yarn.end - 1; i++) {
          theta[i] -= thetaUpdate[i];
          maxUpdate = std::max(maxUpdate, thetaUpdate[i]);
        }
      }

      if (maxUpdate > params.materialFrameTolerance) shouldContinue = true;

      if (params.debug)
        SPDLOG_INFO("Solve for material frame iteration {}, maximum update is {}", iter, maxUpdate);

      statistics.materialFrameUpdateCount++;
    }
  }

  {
    EASY_BLOCK("Update material frame");
    for (const auto& yarn : yarns.yarns) {
      // FIXME: Should start from 0?
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, [this](int thread_id, size_t i) {
        m1.row(i) = std::cos(theta[i]) * u.row(i) + std::sin(theta[i]) * v.row(i);
        m2.row(i) = -std::sin(theta[i]) * u.row(i) + std::cos(theta[i]) * v.row(i);
        });
    }
  }

  {
    EASY_BLOCK("updateThetaHat");
    for (const auto& yarn : yarns.yarns) {
      threading::runParallelFor(thread_pool, yarn.begin + 1, yarn.end - 1, [this](int thread_id, size_t i) {
        double newValue = newThetaHat(e, u, i);
        if (newValue + pi * thetaHatOffset[i] - thetaHat[i] < -pi / 2) {
          thetaHatOffset[i] += 1;
        }
        else if (newValue + pi * thetaHatOffset[i] - thetaHat[i] > pi / 2) {
          thetaHatOffset[i] -= 1;
        }
        assert(std::abs(newValue + pi * thetaHatOffset[i] - thetaHat[i]) < pi / 2);
        thetaHat[i] = newValue + pi * thetaHatOffset[i];
        // FIXME: Should adjust on 2pi ?
        });
    }
  }
}

void DiscreteSimulator::applyGlobalDamping() {
  EASY_FUNCTION();

  for (const auto& yarn : yarns.yarns) {
    threading::runParallelFor(thread_pool, yarn.begin, yarn.end, [this](int thread_id, size_t i) {
      double v = pointAt(dQ, i).norm();
      pointAt(dQ, i) *= std::max(0.0, v - params.kGlobalDamping) / v;
      });
  }
}

void DiscreteSimulator::setUpConstraints() {
  // Length Constraints
  if (params.enableLengthConstrain) {
    for (const auto& yarn : yarns.yarns) {
      for (size_t i = yarn.begin; i < yarn.end - 1; i++) {
        addSegmentLengthConstraint(i);
      }
    }
  }

  // TODO: remove hard-coded pin
  // std::vector<int> pins = { 1203,
    // 1195, 1179, 1183, 1169, 1165, 1151, 1155, 1137, 1141, 1123, 1127, 1109, 1113, 1095, 1099,
    // 1089, 1024, 798, 578, 382, 210, 62,
    // 909, 687, 479, 295, 135, 0, 1, 10,
    // 12, 24, 36, 48 };

  // for (auto i : pins) {
  //   addPinConstraint(i, pointAt(Q, i));
  // }
}

}  // namespace Simulator
