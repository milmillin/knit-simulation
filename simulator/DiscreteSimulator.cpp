#include "DiscreteSimulator.h"

#include <iostream>
#include <functional>

#include <Eigen/SparseLU>
#include "../easy_profiler_stub.h"

#include "macros.h"
#include "Helper.h"
#include "./threading/threading.h"

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

static inline Eigen::Vector3f vec(RowMatrixX3f &v, int index) {
  return Eigen::Vector3f(v.row(index).transpose());
}

void DiscreteSimulator::curvatureBinormalTask
    (int thread_id, int start_index, int end_index) {
  for (int i = start_index; i < end_index; i++) {
    Eigen::Vector3f a = 2 * vec(e, i - 1).cross(vec(e, i));
    float b = segmentLength[i-1] * segmentLength[i];
    float c = e.row(i-1).dot(e.row(i));
    curvatureBinormal.row(i) = a.transpose() / (b + c);
  }
}

Eigen::Vector2f DiscreteSimulator::omega(int i, int j) {
  Eigen::Vector3f kb = vec(curvatureBinormal, i);
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
  e.resize(m - 1, Eigen::NoChange);
  m1.resize(m - 1, Eigen::NoChange);
  m2.resize(m - 1, Eigen::NoChange);
  restOmega.resize(m - 1, Eigen::NoChange);
  restOmega_1.resize(m - 1, Eigen::NoChange);
  curvatureBinormal.resize(m - 1, Eigen::NoChange);

  e.setZero();
  m1.setZero();
  m2.setZero();
  restOmega.setZero();
  restOmega_1.setZero();
  curvatureBinormal.setZero();


  for (int i = 0; i < m - 1; i++) {
    Eigen::Matrix3f empty;
    empty.setZero();
    gradCurvatureBinormal.push_back(std::vector<Eigen::Matrix3f>(3, empty));
  }

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
    restOmega.row(i) = omega(i, i).transpose();
    restOmega_1.row(i) = omega(i, i - 1).transpose();
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

void DiscreteSimulator::gradCurvatureBinormalTask
    (int thread_id, int start_index, int end_index) {
  Eigen::Matrix3f a;
  Eigen::RowVector3f b;
  float c;
  for (int i = start_index; i < end_index; i++) {
    a = 2 * crossMatrix(e.row(i)) + 2 * crossMatrix(e.row(i-1));
    b = (e.row(i) - e.row(i-1));
    c = segmentLength[i-1] * segmentLength[i] + e.row(i-1).dot(e.row(i));
    for (int k = std::max(1, i - 1); k <= std::min(m - 2, i + 1); k++) {
      gradCurvatureBinormal[i][k - (i-1)] =
        - (a + vec(curvatureBinormal, k) * b) / c;
    }
  }
}

static float newThetaHat(RowMatrixX3f &e, RowMatrixX3f &u, int i) {
  Eigen::Vector3f newU = parallelTransport(u.row(i-1), e.row(i - 1), e.row(i));
  return std::atan2(newU.cross(vec(u, i)).norm(), newU.dot(u.row(i)));
}

void DiscreteSimulator::bendingForceTask
    (int thread_id, int start_index, int end_index) {
  for (int i = start_index; i < end_index; i++) {
    Eigen::Vector3f force;
    force.setZero();

    for (int k = std::max(1, i - 1); k <= std::min(m - 2, i + 1); k++) {
      float l = segmentLength[k-1] + segmentLength[k];
      Eigen::Vector3f kb = vec(curvatureBinormal, k);
      for (int j = k - 1; j <= k; j++) {
        Eigen::MatrixXf coeff(2, 3);
        coeff.row(0) = m2.row(j);
        coeff.row(1) = -m1.row(j);
        auto omegaBar = (j == k) ?
          restOmega.row(k).transpose() :
          restOmega_1.row(k).transpose();
        force -= (1.0f / l)
          * (coeff * gradCurvatureBinormal[i][k - (i - 1)]).transpose()
          * (omega(k, j) - omegaBar);
      }
    }
    pointAt(F, i) += params.kBend * force;
  }
}

void DiscreteSimulator::twistingForceTask
    (int thread_id, int start_index, int end_index) {
  for (int i = start_index; i < end_index; i++) {
    Eigen::Vector3f dTheta_dQi_1 = vec(curvatureBinormal, i) / 2 / e.row(i - 1).norm();
    Eigen::Vector3f dTheta_dQi = vec(curvatureBinormal, i) / 2 / e.row(i).norm();
    float coeff = 2 * (theta[i] - theta[i + 1] - thetaHat[i]) / (e.row(i - 1).norm() + e.row(i).norm());
    pointAt(F, i) += params.kTwist * coeff * (dTheta_dQi_1 - dTheta_dQi);
  }
}

void DiscreteSimulator::applyBendingForce() {
  EASY_FUNCTION();
  using namespace std::placeholders;

  int step =
    (m < 200) ?
    (1 + (m / thread_pool.size()))
    : 200;

  {
    EASY_BLOCK("curvatureBinormalTask");
    auto task = std::bind(&DiscreteSimulator::curvatureBinormalTask,
                          this, _1, _2, _3);
    threading::runSequentialJob(thread_pool, task, 1, m-1, step);
  }

  {
    EASY_BLOCK("gradCurvatureBinormalTask");
    auto task = std::bind(&DiscreteSimulator::gradCurvatureBinormalTask,
                          this, _1, _2, _3);
    threading::runSequentialJob(thread_pool, task, 1, m-1, step);
  }

  {
    EASY_BLOCK("bendingForceTask");
    auto task = std::bind(&DiscreteSimulator::bendingForceTask,
                          this, _1, _2, _3);
    threading::runSequentialJob(thread_pool, task, 1, m-1, step);
  }

  {
    EASY_BLOCK("twistingForceTask");
    auto task = std::bind(&DiscreteSimulator::twistingForceTask,
                          this, _1, _2, _3);
    threading::runSequentialJob(thread_pool, task, 1, m-2, step);
  }
}

void DiscreteSimulator::updateBendingForceMetadata() {
  EASY_FUNCTION();

  if (m <= 2) {
    return;
  }

  int step = m < 500 ?
    1 + (m / thread_pool.size())
    : 500;

  // Update frames
  EASY_BLOCK("Update frame")
  threading::runSequentialJob(thread_pool,
    [this](int thread_id, int start, int end) {
      for (int i = start; i < end; i++) {
        Eigen::Vector3f newE = pointAt(Q, i + 1) - pointAt(Q, i);
        u.row(i) = parallelTransport(u.row(i), e.row(i), newE).normalized();
        v.row(i) = newE.cross(vec(u, i)).normalized();
        e.row(i) = newE.transpose();
      }
    }, 1, m - 1, step);
  EASY_END_BLOCK;

  EASY_BLOCK("Solve theta");
  bool shouldContinue = true;
  for (int iter = 0; shouldContinue && iter < 100; iter++) {
    shouldContinue = false;

    std::vector<float> thetaUpdate(m - 1, 0.0f);
    threading::runSequentialJob(thread_pool,
      [this, &thetaUpdate](int thread_id, int start, int end) {
        for (int i = start; i < end; i++) {
          float li = e.row(i).norm() + e.row(i-1).norm();
          thetaUpdate[i] = params.kTwist * 2 * (theta[i] - theta[i - 1] - thetaHat[i]) / li;
          if (i < m - 2) {
            float li_1 = e.row(i).norm() + e.row(i+1).norm();
            thetaUpdate[i] -= params.kTwist * 2 * (theta[i+1] - theta[i] - thetaHat[i]) / li_1;
          }
        }
      }, 1, m - 1, step);


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
  EASY_END_BLOCK;

  EASY_BLOCK("Update material frame");
  threading::runSequentialJob(thread_pool,
    [this](int thread_id, int start, int end) {
      for (int i = start; i < end; i++) {
        m1.row(i) = std::cos(theta[i]) * u.row(i) + std::sin(theta[i]) * v.row(i);
        m2.row(i) = std::sin(theta[i]) * u.row(i) + std::cos(theta[i]) * v.row(i);
      }
    }, 1, m - 1, step);
  EASY_END_BLOCK;

  EASY_BLOCK("updateThetaHat");
  threading::runSequentialJob(thread_pool,
    [this](int thread_id, int start, int end) {
      for (int i = start; i < end; i++) {
        float newValue = newThetaHat(e, u, i);
        if (newValue + pi * thetaHatOffset[i] - thetaHat[i] < - pi / 2) {
          thetaHatOffset[i] += 1;
        } else if (newValue + pi * thetaHatOffset[i] - thetaHat[i] > pi / 2) {
          thetaHatOffset[i] -= 1;
        }
        assert(std::abs(newValue + pi * thetaHatOffset[i] - thetaHat[i]) < pi / 2);
        thetaHat[i] = newValue + pi * thetaHatOffset[i];
      }
    }, 1, m - 1, step);
  EASY_END_BLOCK;
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
