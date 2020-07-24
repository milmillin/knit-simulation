#include "DiscreteSimulator.h"

#include "macros.h"
#include "Helper.h"
#include "./AABB.h"

#include <iostream>

#include<Eigen/SparseLU>
#include<glm/gtx/norm.hpp>

namespace simulator {

DiscreteSimulator::DiscreteSimulator() {
  yarns = file_format::YarnRepr();
  yarns.yarns.push_back(file_format::Yarn());
}

DiscreteSimulator::DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params) {
  this->yarns = yarns;
  this->params = params;

  auto &Q = this->yarns.yarns[0].points;

  // TODO: remove hard-coded control points
  pinControlPoints.push_back(0);

  // TODO: remove hard-coded slicing
  Eigen::MatrixXf newMatrix(10, 3);
  newMatrix = Q.block(0, 0, 20, 3);
  Q = newMatrix;

  // Initialize speed
  dQ = Eigen::MatrixXf(Q.rows(), 3);
  dQ.setZero();

  // Initialize force
  ddQ = Eigen::MatrixXf(Q.rows(), 3);

  // Initialize rest length
  for (int i = 0; i < Q.rows() - 1; i++) {
    glm::vec3 p1 = POINT_FROM_ROW(Q, i);
    glm::vec3 p2 = POINT_FROM_ROW(Q, i + 1);
    restLength.push_back(glm::length(p2 - p1));
  }

  // Initialize constrain
  constrain = Eigen::VectorXf(nConstrain(), 1);
  dConstrain = Eigen::SparseMatrix<float>(nConstrain(), 3 * Q.rows());

  // Initialize pin position
  for (int pinPoint : pinControlPoints) {
    pinControlPointsPosition.push_back(POINT_FROM_ROW(Q, pinPoint));
  }

  // Initialize AABB tree
  collisionTree = aabb::Tree(3, 0.05, Q.rows() - 4, true);
  std::vector<double> lowerBound;
  std::vector<double> upperBound;
  for (int i = 0; i < Q.rows() - 3; i++) {
    catmullRomBoundingBox(Q, i, &lowerBound, &upperBound, yarns.yarns[0].radius);
    collisionTree.insertParticle((unsigned)i, lowerBound, upperBound);
  }

  // Initialize frame
  initBendingForceMetadata();
}

void DiscreteSimulator::step() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < params.steps; i++) {
    // Calculate acceleration
    ddQ.setZero();
    applyGravity();
    applyContactForce();
    applyBendingForce();

    // Calculate velocity
    dQ += ddQ * params.h;
    applyGroundVelocityFilter();
    applyGlobalDamping();

    // Calculate position
    Eigen::MatrixXf originalQ = Q;
    Q += dQ * params.h;

    // Apply constrains to position
    fastProjection();
    // Also update the velocity after applying the constrains
    dQ = (Q - originalQ) / params.h;

    // Update AABB tree
    std::vector<double> lowerBound;
    std::vector<double> upperBound;
    for (int i = 0; i < Q.rows() - 3; i++) {
      catmullRomBoundingBox(Q, i, &lowerBound, &upperBound, yarns.yarns[0].radius);
      collisionTree.updateParticle((unsigned)i, lowerBound, upperBound);
    }

    // Update frame
    updateBendingForceMetadata();
  }
}

void DiscreteSimulator::applyGravity() {
  auto &Q = yarns.yarns[0].points;
	ddQ += Eigen::Vector3f(0, -params.gravity, 0).transpose().replicate(Q.rows(), 1);
}

void DiscreteSimulator::applyGroundVelocityFilter() {
  auto &Q = yarns.yarns[0].points;
  for (int i = 0; i < Q.rows(); i++) {
    if (Q(i, 1) < params.groundHeight && dQ(i, 1) < 0) {
      dQ(i, 0) *= params.groundFiction;
      dQ(i, 1) = 0;
      dQ(i, 2) *= params.groundFiction;
      Q(i, 1) = params.groundHeight;
    }
  }
}

#define SAMPLES 8

static glm::vec3 f(glm::vec3 gradient, float l, float r) {
  return -4 * r * r / l / l / l / l * gradient + 0.25f / r / r * gradient;
}

void DiscreteSimulator::applyContactForce() {
  auto &Q = yarns.yarns[0].points;
  float radius = yarns.yarns[0].radius;

  for (int i = 0; i < Q.rows() - 4; i++) {
    auto intersections = collisionTree.query(i);
    for (int j : intersections) {
      if (abs(i - j) <= 1) {
        continue;
      }

      float ds = 1.0f / SAMPLES;
      for (int m = 0; m < SAMPLES; m++) {
        float s1 = (float) m / SAMPLES;
        for (int n = 0; n < SAMPLES; n++) {
          float s2 = (float) n / SAMPLES;
          glm::vec3 p1 = catmullRomSample(Q, i, s1);
          glm::vec3 p2 = catmullRomSample(Q, j, s2);
          glm::vec3 dir = p2 - p1;
          float length = glm::length(dir);

          if (length < 2 * radius) {
            // std::cout << "Collision" <<i <<" " << j<< std::endl;
            glm::vec3 f0 = params.kContact * ds * ds * f(2 * simulator::b1(s1) * dir, length, radius);
            glm::vec3 f1 = params.kContact * ds * ds * f(2 * simulator::b2(s1) * dir, length, radius);
            glm::vec3 f2 = params.kContact * ds * ds * f(2 * simulator::b3(s1) * dir, length, radius);
            glm::vec3 f3 = params.kContact * ds * ds * f(2 * simulator::b4(s1) * dir, length, radius);

            ADD_TO_ROW(ddQ, i + 0, f0);
            ADD_TO_ROW(ddQ, i + 1, f1);
            ADD_TO_ROW(ddQ, i + 2, f2);
            ADD_TO_ROW(ddQ, i + 3, f3);
          }
        }
      }
    }
  }
  // std::cout << std::endl;
}

inline static void parallelTransport(Eigen::Vector3f tangent, Eigen::Vector3f oldUp,
    Eigen::Vector3f *newUp, Eigen::Vector3f *newCross) {
  *newCross = tangent.cross(oldUp).normalized();
  *newUp = newCross->cross(tangent).normalized();
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
  auto &Q = yarns.yarns[0].points;

  int nPoints = Q.rows();

  // No bending energy
  if (nPoints <= 2) {
    return;
  }

  // Allocate memory
  e.resize(nPoints - 1, 3);
  m1.resize(nPoints - 1, 3);
  m2.resize(nPoints - 1, 3);
  restOmega.resize(nPoints, 2);
  restOmega_1.resize(nPoints, 2);

  // Initialize tangent
  for (int i = 0; i < nPoints - 1; i++) {
    e.row(i) = Q.row(i + 1) - Q.row(i);
  }

  // Initialize direction 1 with arbitrary vector that's normal to u.row(0)
  if (std::abs(Q(0, 0) + Q(0, 1)) < std::abs(Q(0, 2))) {
    m1.row(0) = Eigen::Vector3f(1, 0, 0).cross(vec(e, 0));
  } else {
    m1.row(0) = Eigen::Vector3f(0, 0, 1).cross(vec(e, 0));
  }
  m1.row(0).normalize();

  // Initialize direction 2
  m2.row(0) = vec(e, 0).cross(vec(m1, 0)).normalized();

  // Fill in all frames
  for (int i = 1; i < nPoints - 1; i++) {
    Eigen::Vector3f newM1, newM2;
    parallelTransport(e.row(i), m1.row(i - 1), &newM1, &newM2);
    m1.row(i) = newM1;
    m2.row(i) = newM2;
  }

  // Calculate rest omega
  for (int i = 1; i < nPoints - 1; i++) {
    restOmega.row(i) = omega(e, m1, m2, i, i).transpose();
    restOmega_1.row(i) = omega(e, m1, m2, i, i - 1).transpose();
  }
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

void DiscreteSimulator::applyBendingForce() {
  auto &Q = yarns.yarns[0].points;

  // TODO: why starting from 1?
  for (int i = 1; i < Q.rows() - 1; i++) {
    Eigen::Vector3f force;
    force.setZero();
    for (int k = std::max(1, i - 1); k <= std::min((int)Q.rows() - 2, i + 1); k++) {
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
    ddQ.row(i) += params.kBend * force;
  }
}

void DiscreteSimulator::updateBendingForceMetadata() {
  auto &Q = yarns.yarns[0].points;

  int nPoints = Q.rows();

  if (nPoints <= 2) {
    return;
  }

  // Update tangent
  for (int i = 0; i < nPoints - 1; i++) {
    e.row(i) = Q.row(i + 1) - Q.row(i);
  }

  // Update frames
  for (int i = 1; i < nPoints - 1; i++) {
    Eigen::Vector3f newM1, newM2;
    parallelTransport(e.row(i), m1.row(i), &newM1, &newM2);
    m1.row(i) = newM1;
    m2.row(i) = newM2;
  }
}

void DiscreteSimulator::fastProjection() {
  auto &Q = yarns.yarns[0].points;

  // Current iteration
  int nIteration = 0;
  do {
    // Calculate constrain
    dConstrain.setZero();
    nextConstrainID = 0;
    applyLengthConstrain();
    applyPinConstrain();

    // Solve for lambda
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;

    solver.compute(dConstrain * dConstrain.transpose());
    if (solver.info() != Eigen::Success) {
      std::cerr << "Warning: failed to solve constrain at iteration "
        << nIteration << ". Stopping constrain fitting early." << std::endl;
      break;
    }

    Eigen::VectorXf lambda = solver.solve(constrain);
    if (solver.info() != Eigen::Success) {
      std::cerr << "Warning: failed to solve constrain at iteration "
        << nIteration << ". Stopping constrain fitting early." << std::endl;
      break;
    }

    // Evaluate position change
    Eigen::VectorXf deltaX = dConstrain.transpose() * lambda;
    Eigen::Map<Eigen::MatrixXf> mDeltaX(deltaX.data(), Q.rows(), 3);

    // Apply position change
    Q += mDeltaX;

    if (params.debug) {
      std::cout << "Iteration: " << nIteration << \
      " Constrain error: "<< constrain.norm() << std::endl;
    }
    // Check for termination
  } while (constrain.norm() / constrain.rows() > params.fastProjErrorCutoff
    && nIteration < params.fastProjMaxIter);
}

int DiscreteSimulator::pointIndex(int pointID, int axis) {
  return axis * yarns.yarns[0].points.rows() + pointID;
}

int DiscreteSimulator::nConstrain() {
  return pinControlPoints.size() + restLength.size();
}

void DiscreteSimulator::applyLengthConstrain() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < Q.rows() - 1; i++) {
    // Calculate constrain
    glm::vec3 p1 = POINT_FROM_ROW(Q, i);
    glm::vec3 p2 = POINT_FROM_ROW(Q, i + 1);
    glm::vec3 direction = p1 - p2;
    constrain(nextConstrainID) = glm::length2(direction) / restLength[i] - restLength[i];

    // Calculate constrain gradient
    glm::vec3 gradient = 2 / restLength[i] * direction;

    // Apply gradient
    dConstrain.coeffRef(nextConstrainID, pointIndex(i, 0)) -= gradient.x;
    dConstrain.coeffRef(nextConstrainID, pointIndex(i, 1)) -= gradient.y;
    dConstrain.coeffRef(nextConstrainID, pointIndex(i, 2)) -= gradient.z;

    dConstrain.coeffRef(nextConstrainID, pointIndex(i+1, 0)) += gradient.x;
    dConstrain.coeffRef(nextConstrainID, pointIndex(i+1, 1)) += gradient.y;
    dConstrain.coeffRef(nextConstrainID, pointIndex(i+1, 2)) += gradient.z;

    nextConstrainID++;
  }
}

void DiscreteSimulator::applyPinConstrain() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < pinControlPoints.size(); i++) {
    // Calculate constrain
    glm::vec3 p = POINT_FROM_ROW(Q, pinControlPoints[i]);
    constrain(nextConstrainID) = glm::length2(p - pinControlPointsPosition[i]);

    // Calculate constrain gradient
    glm::vec3 gradient = (p - pinControlPointsPosition[i]) * 2.0f;

    // Apply gradient
    dConstrain.coeffRef(nextConstrainID, pointIndex(pinControlPoints[i], 0)) -= gradient.x;
    dConstrain.coeffRef(nextConstrainID, pointIndex(pinControlPoints[i], 1)) -= gradient.y;
    dConstrain.coeffRef(nextConstrainID, pointIndex(pinControlPoints[i], 2)) -= gradient.z;

    nextConstrainID++;
  }
}

void DiscreteSimulator::applyGlobalDamping() {
  dQ *= exp(-params.kGlobal * params.h);
}

}  // namespace Simulator
