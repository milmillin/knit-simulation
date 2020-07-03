#include "Simulator.h"

#include "SimulatorParams.h";

#include <Eigen/Core>

namespace simulator {

Simulator::Simulator(Eigen::MatrixXf q_, SimulatorParams params_) : q(q_), params(params_)  {
  // initial velocity is zero
  qD = Eigen::MatrixXf::Zero(q_.rows(), q_.cols());
}

int Simulator::getNumPoints() const {
  return q.rows();
}

const Eigen::MatrixXf& Simulator::getCurrentPoints() {
  return q;
}

void Simulator::step() {
  //TODO:
}

};  // namespace simulator
