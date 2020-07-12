#include "DiscreteSimulator.h"

#include "macros.h"

namespace simulator {

DiscreteSimulator::DiscreteSimulator() {
  yarns = file_format::YarnRepr();
  yarns.yarns.push_back(file_format::Yarn());
}

DiscreteSimulator::DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params) {
  this->yarns = yarns;
  this->params = params;
  auto &Q = yarns.yarns[0].points;
  dQ = Eigen::MatrixXf(Q.rows(), 3);
  dQ.setZero();
  ddQ = Eigen::MatrixXf(Q.rows(), 3);

  for (int i = 0; i < Q.rows() - 1; i++) {
    glm::vec3 p1 = POINT_FROM_ROW(Q, i);
    glm::vec3 p2 = POINT_FROM_ROW(Q, i + 1);
    restLength.push_back(glm::length(p2 - p1));
  }
}

void DiscreteSimulator::step() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < params.steps; i++) {
    // Calculate acceleration
    ddQ.setZero();
    applyGravity();
    applyContactForce();
    applyLengthConstrain();

    // Calculate velocity
    dQ += ddQ * params.h;
    applyGroundVelocityFilter();

    // Calculate position
    Q += dQ * params.h;
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
    }
  }
}

void DiscreteSimulator::applyContactForce() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < Q.rows(); i++) {
    for (int j = i + 1; j < Q.rows(); j++) {
      glm::vec3 p1 = POINT_FROM_ROW(Q, i);
      glm::vec3 p2 = POINT_FROM_ROW(Q, j);
      glm::vec3 direction = p2 - p1;
      float distance = glm::length(direction);
      if (distance < 2 * yarns.yarns[0].radius) {
        distance = distance / 2 / yarns.yarns[0].radius;
        float force = 1 / distance / distance + distance * distance  - 2;
        force *= params.kContact;
        direction = glm::normalize(direction);
        glm::vec3 f2 = direction * force;
        glm::vec3 f1 = direction * force;
        ADD_TO_ROW(ddQ, i, f1);
        ADD_TO_ROW(ddQ, i, f2);
      }
    }
  }
}

void DiscreteSimulator::applyLengthConstrain() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < Q.rows() - 1; i++) {
    glm::vec3 p1 = POINT_FROM_ROW(Q, i);
    glm::vec3 p2 = POINT_FROM_ROW(Q, i + 1);
    glm::vec3 direction = p2 - p1;
    float offset = glm::length(direction) - restLength[i];
    direction = glm::normalize(direction);
    glm::vec3 force = 0.5f * params.kLen * offset * offset * direction;
    if (offset < restLength[i]) {
      SUBTRACT_FROM_ROW(ddQ, i, force);
      ADD_TO_ROW(ddQ, i + 1, force);
    } else {
      SUBTRACT_FROM_ROW(ddQ, i + 1, force);
      ADD_TO_ROW(ddQ, i, force);
    }
  }
}

}  // namespace Simulator