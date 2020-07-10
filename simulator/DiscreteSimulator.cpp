#include "DiscreteSimulator.h"

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
}

void DiscreteSimulator::step() {
  auto &Q = yarns.yarns[0].points;

  for (int i = 0; i < params.steps; i++) {
    // Calculate acceleration
    ddQ.setZero();
    applyGravity();

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

}  // namespace Simulator