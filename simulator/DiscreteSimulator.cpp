#include "DiscreteSimulator.h"

#include "macros.h"

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
}

void DiscreteSimulator::step() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < params.steps; i++) {
    // Calculate acceleration
    ddQ.setZero();
    applyGravity();
    applyContactForce();

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

void DiscreteSimulator::applyContactForce() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < Q.rows(); i++) {
    for (int j = i + 1; j < Q.rows(); j++) {
      // Check contact
      glm::vec3 p1 = POINT_FROM_ROW(Q, i);
      glm::vec3 p2 = POINT_FROM_ROW(Q, j);
      glm::vec3 direction = p2 - p1;
      float distance = glm::length(direction);
      if (distance < 2 * yarns.yarns[0].radius) {
        // Relative distance
        distance = distance / 2 / yarns.yarns[0].radius;

        // Calculate force magnitude
        float force = 1 / distance / distance + distance * distance  - 2;
        force *= params.kContact;

        // Calculate force
        direction = glm::normalize(direction);
        glm::vec3 f = direction * force;

        // Apply force
        ADD_TO_ROW(ddQ, i + 1, f);
        SUBTRACT_FROM_ROW(ddQ, i, f);
      }
    }
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
