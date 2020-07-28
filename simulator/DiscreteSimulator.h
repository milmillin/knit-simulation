#pragma once

#include "./BaseSimulator.h"

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"

namespace simulator {

class DiscreteSimulator : public BaseSimulator {
public:

  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params);

  // Force
  void applyGravity();

  void initBendingForceMetadata();
  void applyBendingForce();
  void applyTwistingForce();
  void updateBendingForceMetadata();

  // Velocity filter
  void applyGroundVelocityFilter();
  void applyGlobalDamping();


 private:
  // === Constrains ===
  // Index of control points to be fixed
  std::vector<int> pinControlPoints;

  // === Bending Force ===
  Eigen::MatrixXf e;
  Eigen::MatrixXf m1;
  Eigen::MatrixXf m2;
  Eigen::MatrixXf u;
  Eigen::MatrixXf v;
  std::vector<float> theta;
  std::vector<float> thetaHat;
  std::vector<int> thetaHatOffset;
  Eigen::MatrixXf restOmega;
  Eigen::MatrixXf restOmega_1;

protected:
  // Simulates next timestep.
  void stepImpl(const StateGetter& cancelled) override;
  // Called after each step's calculation
  void postStep(const std::function<bool()>& cancelled) override;

  // Set up constraints
  void setUpConstraints() override;
};

} // namespace simulator 
