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
  for (int i = 0; i < Q.rows() - 4; i++) {
    catmullRomBoundingBox(Q, i, &lowerBound, &upperBound, yarns.yarns[0].radius);
    collisionTree.insertParticle((unsigned)i, lowerBound, upperBound);
  }
}

void DiscreteSimulator::step(const std::function<bool()>& cancelled) {
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

    // Update AABB tree
    std::vector<double> lowerBound;
    std::vector<double> upperBound;
    for (int i = 0; i < Q.rows() - 4; i++) {
      catmullRomBoundingBox(Q, i, &lowerBound, &upperBound, yarns.yarns[0].radius);
      collisionTree.updateParticle((unsigned)i, lowerBound, upperBound);
    }
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
