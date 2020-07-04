#pragma once

#include "SimulatorParams.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace simulator {

// Simulator contains all yarn simulation functionalities
//
// A yarn is modeled by a Catmull-Rom spline with #m control points
// denoted by q[i], i from 0 to m - 1. Each segment is governed by
// 4 control points, so there are #N = #m - 3 segments (0-indexed).
// segment[i] is governed by q[i], q[i + 1], q[i + 2], q[i + 3].
class Simulator {
private:
  size_t m;
  Eigen::MatrixXf q;
  SimulatorParams params;
  Eigen::SparseMatrix<float> M;
  Eigen::SparseMatrix<float> MInverse;

  // first-derivative of q
  Eigen::MatrixXf qD;

  // gradient of positional energy
  Eigen::MatrixXf gradE;

  // gradient of damping energy
  Eigen::MatrixXf gradD;

  // external force
  Eigen::MatrixXf f;

  void constructMassMatrix();
public:
  // Constructs a new Simulator
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  Simulator(Eigen::MatrixXf q_, SimulatorParams params_);

  // Returns #m x 3 matrix containing the current control points.
  Eigen::MatrixXf getControlPoints() const;

  // Simulates next timestep.
  void step();
};

}; // namespace simulator
