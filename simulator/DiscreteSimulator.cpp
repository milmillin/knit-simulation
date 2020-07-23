#include "DiscreteSimulator.h"

#include "macros.h"
#include "Helper.h"
#include "./AABB.h"

#include <iostream>

#include<Eigen/SparseLU>
#include<glm/gtx/norm.hpp>

namespace simulator {

DiscreteSimulator::DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params):
  BaseSimulator(yarns, params) 
{
  initialize();
}

void DiscreteSimulator::stepImpl(const StateGetter& cancelled) {
  // Calculate acceleration
  F.setZero();
  applyGravity();
  applyContactForce(cancelled);

  Eigen::MatrixXf ddQ = F; // invM * F

  // Calculate velocity
  dQ += ddQ * params.h;
  applyGroundVelocityFilter();
  applyGlobalDamping();

  // Calculate position
  Eigen::MatrixXf originalQ = Q;
  Q += dQ * params.h;
}

void DiscreteSimulator::applyGravity() {
  auto &Q = yarns.yarns[0].points;
	F += Eigen::Vector3f(0, -params.gravity, 0).replicate(m, 1);
}

void DiscreteSimulator::applyGroundVelocityFilter() {
  for (int i = 0; i < m; i++) {
    if (coordAt(Q, i, 1) < params.groundHeight && coordAt(dQ, i, 1) < 0) {
      coordAt(dQ, i, 0) *= params.groundFiction;
      coordAt(dQ, i, 1) = 0;
      coordAt(dQ, i, 2) *= params.groundFiction;
      coordAt(Q, i, 1) = params.groundHeight;
    }
  }
}

void DiscreteSimulator::applyGlobalDamping() {
  dQ *= exp(-params.kGlobal * params.h);
}

void DiscreteSimulator::setUpConstraints() {
  // Length Constraints
  for (int i = 0; i < m - 1; i++) {
    addSegmentLengthConstraint(i);
  }

  // Add pin constraints
  addPinConstraint(0, pointAt(Q, 0));
}

}  // namespace Simulator
