#include "DiscreteSimulator.h"

#include "macros.h"
#include "Helper.h"

#include <iostream>

#include <Eigen/SparseLU>
#include "../easy_profiler_stub.h"

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
  applyGravity();
  applyContactForce(cancelled);
  applyBendingForce();

  Eigen::MatrixXf ddQ = F; // invM * F

  // Calculate velocity
  dQ += ddQ * params.h;
  applyGroundVelocityFilter();
  applyGlobalDamping();

  // Calculate position
  Q += dQ * params.h;
}

void DiscreteSimulator::postStep(const StateGetter& cancelled) {
  updateBendingForceMetadata();
}

void DiscreteSimulator::applyGravity() {
  EASY_FUNCTION();

	F += Eigen::Vector3f(0, -params.gravity, 0).replicate(m, 1);
}

void DiscreteSimulator::applyGroundVelocityFilter() {
  EASY_FUNCTION();

  for (int i = 0; i < m; i++) {
    if (coordAt(Q, i, 1) < params.groundHeight && coordAt(dQ, i, 1) < 0) {
      coordAt(dQ, i, 0) *= params.groundFriction;
      coordAt(dQ, i, 1) = 0;
      coordAt(dQ, i, 2) *= params.groundFriction;
      coordAt(Q, i, 1) = params.groundHeight;
    }
  }
}

inline static Eigen::Vector3f parallelTransport(const Eigen::Vector3f& u, const Eigen::Vector3f& e1, const Eigen::Vector3f& e2) {
    Eigen::Vector3f t1 = e1 / e1.norm();
    Eigen::Vector3f t2 = e2 / e2.norm();
    Eigen::Vector3f n = t1.cross(t2);
    if (n.norm() < 1e-10)
        return u;
    n /= n.norm();
    Eigen::Vector3f p1 = n.cross(t1);
    Eigen::Vector3f p2 = n.cross(t2);
    return u.dot(n)*n + u.dot(t1)*t2 + u.dot(p1)*p2;
}

static inline Eigen::Vector3f vec(Eigen::MatrixXf &v, int index) {
  return Eigen::Vector3f(v.row(index).transpose());
}

static inline Eigen::Vector3f curvatureBinormal(Eigen::MatrixXf &e, int i) {
  Eigen::Vector3f a = 2 * vec(e, i - 1).cross(vec(e, i));
  float b = e.row(i-1).norm() * e.row(i).norm();
  float c = e.row(i-1).dot(e.row(i));
  return a / (b + c);
  // return 2 * vec(u, i - 1).cross(vec(u, i))
  //   / (u.row(i-1).norm() * u.row(i).norm() + u.row(i-1).dot(u.row(i)));
}

static inline Eigen::Vector2f omega(Eigen::MatrixXf &e, Eigen::MatrixXf &m1, Eigen::MatrixXf &m2,
  int i, int j) {
  Eigen::Vector3f kb = curvatureBinormal(e, i);
  Eigen::Vector2f result;
  result(0) = kb.dot(m2.row(j));
  result(1) = -kb.dot(m1.row(j));
  return result;
}

void DiscreteSimulator::initBendingForceMetadata() {
  EASY_FUNCTION();

  // No bending energy
  if (m <= 2) {
    return;
  }

  // Allocate memory
  e.resize(m - 1, 3);
  m1.resize(m - 1, 3);
  m2.resize(m - 1, 3);
  restOmega.resize(m, 2);
  restOmega_1.resize(m, 2);

  // Initialize tangent
  for (int i = 0; i < m - 1; i++) {
    e.row(i) = pointAt(Q, i + 1) - pointAt(Q, i);
  }

  // Initialize direction 1 with arbitrary vector that's normal to u.row(0)
  if (std::abs(Q(0, 0) + Q(1, 0)) < std::abs(Q(2, 0))) {
    m1.row(0) = Eigen::Vector3f(1, 0, 0).cross(vec(e, 0));
  } else {
    m1.row(0) = Eigen::Vector3f(0, 0, 1).cross(vec(e, 0));
  }
  m1.row(0).normalize();

  // Initialize direction 2
  m2.row(0) = vec(e, 0).cross(vec(m1, 0)).normalized();

  // Fill in all frames
  for (int i = 1; i < m - 1; i++) {
    m1.row(i) = parallelTransport(m1.row(i-1), e.row(i-1), e.row(i)).normalized();
    m2.row(i) = vec(e, i).cross(vec(m1, i)).normalized();
  }


  // Calculate rest omega
  for (int i = 1; i < m - 1; i++) {
    restOmega.row(i) = omega(e, m1, m2, i, i).transpose();
    restOmega_1.row(i) = omega(e, m1, m2, i, i - 1).transpose();
  }

  // Initialize theta
  u = m1;
  v = m2;
  theta = std::vector<float>(m - 1, 0.0f);
  thetaHat = std::vector<float>(m - 1, 0.0f);
  thetaHatOffset = std::vector<int>(m - 1, 0);
}

static inline Eigen::Matrix3f crossMatrix(Eigen::Vector3f e) {
  Eigen::Matrix3f result;
  result <<
        0, -e(2),  e(1),
     e(2),     0, -e(0),
    -e(1),  e(0),     0;
  return result;
}

static inline Eigen::Matrix3f gradCurvatureBinormal(Eigen::MatrixXf &e, Eigen::Vector3f kb, int i) {
  return - (2 * crossMatrix(e.row(i)) + 2 * crossMatrix(e.row(i-1))
      + kb * (e.row(i) - e.row(i-1)))
    / (e.row(i-1).norm() * e.row(i).norm() + e.row(i-1).dot(e.row(i)));
}

static float newThetaHat(Eigen::MatrixXf &e, Eigen::MatrixXf &u, int i) {
  Eigen::Vector3f newU = parallelTransport(u.row(i-1), e.row(i - 1), e.row(i));
  return std::atan2(newU.cross(vec(u, i)).norm(), newU.dot(u.row(i)));
}

void DiscreteSimulator::applyBendingForce() {
  EASY_FUNCTION();

  for (int i = 1; i < m - 1; i++) {
    Eigen::Vector3f force;
    force.setZero();
    for (int k = std::max(1, i - 1); k <= std::min(m - 2, i + 1); k++) {
      float l = e.row(k - 1).norm() + e.row(k).norm();
      Eigen::Vector3f kb = curvatureBinormal(e, k);
      for (int j = k - 1; j <= k; j++) {
        Eigen::MatrixXf coeff(2, 3);
        coeff.row(0) = m2.row(j);
        coeff.row(1) = -m1.row(j);
        auto omegaBar = (j == k) ?
          restOmega.row(k).transpose() :
          restOmega_1.row(k).transpose();
        force -= (1.0f / l)
          * (coeff * gradCurvatureBinormal(e, kb, i)).transpose()
          * (omega(e, m1, m2, k, j) - omegaBar);
      }
    }
    
    pointAt(F, i) += params.kBend * force;
    Eigen::Vector3f dTheta_dQi_1 = curvatureBinormal(e, i) / 2 / e.row(i - 1).norm();
    Eigen::Vector3f dTheta_dQi = curvatureBinormal(e, i) / 2 / e.row(i).norm();
    float coeff = 2 * (theta[i] - theta[i + 1] - thetaHat[i]) / (e.row(i - 1).norm() + e.row(i).norm());
    pointAt(F, i) += params.kTwist * coeff * (dTheta_dQi_1 - dTheta_dQi);
  }
}

void DiscreteSimulator::updateBendingForceMetadata() {
  EASY_FUNCTION();

  if (m <= 2) {
    return;
  }

  // Update frames
  for (int i = 0; i < m - 1; i++) {
    Eigen::Vector3f newE = pointAt(Q, i + 1) - pointAt(Q, i);
    u.row(i) = parallelTransport(u.row(i), e.row(i), newE).normalized();
    v.row(i) = newE.cross(vec(u, i)).normalized();
    e.row(i) = newE.transpose();
  }

  bool shouldContinue = true;
  for (int iter = 0; shouldContinue && iter < 100; iter++) {
    shouldContinue = false;

    std::vector<float> thetaUpdate(m - 1, 0.0f);
    for (int i = 1; i < m - 1; i++) {
      float li = e.row(i).norm() + e.row(i-1).norm();
      thetaUpdate[i] = params.kTwist * 2 * (theta[i] - theta[i - 1] - thetaHat[i]) / li;
      if (i != m - 2) {
        float li_1 = e.row(i).norm() + e.row(i+1).norm();
        thetaUpdate[i] -= params.kTwist * 2 * (theta[i+1] - theta[i] - thetaHat[i]) / li_1;
      }
    }


    float maxUpdate = 0;
    for (int i = 1; i < m - 1; i++) {
      theta[i] -= thetaUpdate[i];
      maxUpdate = std::max(maxUpdate, thetaUpdate[i]);
    }

    if (maxUpdate > 1e-4) {
      shouldContinue = true;
    }

    if (params.debug) {
      std::cout << "Solve for material frame iteration " << iter << " " << maxUpdate << std::endl;
    }
  }

  for (int i = 0; i < m - 1; i++) {
    m1.row(i) = std::cos(theta[i]) * u.row(i) + std::sin(theta[i]) * v.row(i);
    m2.row(i) = std::sin(theta[i]) * u.row(i) + std::cos(theta[i]) * v.row(i);
  }

  for (int i = 1; i < m - 1; i++) {
    float newValue = newThetaHat(e, u, i);
    if (newValue + pi * thetaHatOffset[i] - thetaHat[i] < - pi / 2) {
      thetaHatOffset[i] += 1;
    } else if (newValue + pi * thetaHatOffset[i] - thetaHat[i] > pi / 2) {
      thetaHatOffset[i] -= 1;
    }
    assert(std::abs(newValue + pi * thetaHatOffset[i] - thetaHat[i]) < pi / 2);
    thetaHat[i] = newValue + pi * thetaHatOffset[i];
  }
}

void DiscreteSimulator::applyGlobalDamping() {
  EASY_FUNCTION();

  dQ *= exp(-params.kGlobal * params.h);
}

void DiscreteSimulator::setUpConstraints() {
  // Length Constraints
  for (int i = 0; i < m - 1; i++) {
    addSegmentLengthConstraint(i);
  }

  // Add pin constraints. TODO: remove hard-coded pin
  addPinConstraint(0, pointAt(Q, 0));
  addPinConstraint(81-59, pointAt(Q, 81-59));
  addPinConstraint(88-59, pointAt(Q, 88-59));
  addPinConstraint(74-59, pointAt(Q, 74-59));
}

}  // namespace Simulator
