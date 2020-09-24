#pragma once

#include <vector>
#include <mutex>
#include <thread>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "./SimulatorParams.h"
#include "./Constraints.h"
#include "./macros.h"
#include "./BaseSimulator.h"
#include "file_format/yarnRepr.h"

namespace simulator{

// Simulator contains all yarn simulation functionalities
//
// A yarn is modeled by a Catmull-Rom spline with #m control points
// denoted by q[i], i from 0 to m - 1. Each segment is governed by
// 4 control points, so there are #N = #m - 3 segments (0-indexed).
// segment[i] is governed by q[i], q[i + 1], q[i + 2], q[i + 3].
class Simulator : public BaseSimulator {
private:
  void calculateBendingEnergy(int thread_id, size_t i);
  void calculateLengthEnergy(int thread_id, size_t i);
  void calculateGlobalDamping(int thread_id, size_t i);
  void calculateGravity();

public:
  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  Simulator(file_format::YarnRepr yarns, SimulatorParams params_);

protected:
  void constructMassMatrix() override;
  void stepImpl(const StateGetter& cancelled) override;
  void setUpConstraints() override;
};

} // namespace simulator
