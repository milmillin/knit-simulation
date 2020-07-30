#pragma once

#include "./BaseSimulator.h"

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"

namespace simulator {

typedef
  Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>
  RowMatrixX3f;
typedef
  Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::RowMajor>
  RowMatrixX2f;

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
  RowMatrixX3f e;
  RowMatrixX3f m1;
  RowMatrixX3f m2;
  RowMatrixX3f u;
  RowMatrixX3f v;
  RowMatrixX3f curvatureBinormal;
  std::vector<float> theta;
  std::vector<float> thetaHat;
  std::vector<int> thetaHatOffset;
  RowMatrixX2f restOmega;
  RowMatrixX2f restOmega_1;

protected:
  // Simulates next timestep.
  void stepImpl(const StateGetter& cancelled) override;
  // Called after each step's calculation
  void postStep(const std::function<bool()>& cancelled) override;

  // Set up constraints
  void setUpConstraints() override;

private:
  // === Bending and twisting ===
  void curvatureBinormalTask(int thread_id, int start_index, int end_index);
  Eigen::Vector2f omega(int i, int j);
};

} // namespace simulator 
