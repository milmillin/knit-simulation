#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "./BaseSimulator.h"
#include "./SimulatorParams.h"
#include "file_format/yarnRepr.h"

namespace simulator {

typedef
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
  RowMatrixX3d;
typedef
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>
  RowMatrixX2d;

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
  void updateBendingForceMetadata();

  // Velocity filter
  void applyGroundVelocityFilter();
  void applyGlobalDamping();


 private:
  // === Bending Force ===
  RowMatrixX3d e;
  RowMatrixX3d m1;
  RowMatrixX3d m2;
  RowMatrixX3d u;
  RowMatrixX3d v;
  RowMatrixX3d curvatureBinormal;
  std::vector<std::vector<Eigen::Matrix3d>> gradCurvatureBinormal;
  std::vector<double> theta;
  std::vector<double> thetaHat;
  std::vector<int> thetaHatOffset;
  RowMatrixX2d restOmega;
  RowMatrixX2d restOmega_1;

protected:
  // Simulates next timestep.
  void stepImpl(const StateGetter& cancelled) override;
  // Called after each step's calculation
  void postStep(const StateGetter& cancelled) override;

  // Set up constraints
  void setUpConstraints() override;

private:
  // === Bending and twisting ===
  void curvatureBinormalTask(int thread_id, size_t index);
  void gradCurvatureBinormalTask(int thread_id, size_t index);
  void bendingForceTask(int thread_id, size_t index);
  void twistingForceTask(int thread_id, size_t index);

  Eigen::Vector2d omega(size_t i, size_t j);
};

} // namespace simulator 
