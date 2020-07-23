#pragma once

#include "./BaseSimulator.h"

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"
#include "./AABB.h"

namespace simulator {

class DiscreteSimulator : public BaseSimulator {
public:

  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params);

  void applyGravity();
  void applyGroundVelocityFilter();
  void applyGlobalDamping();

protected:
  // Simulates next timestep.
  void stepImpl(const StateGetter& cancelled) override;

  // Set up constraints
  void setUpConstraints() override;
};

} // namespace simulatr 

