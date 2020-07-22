#pragma once

#include "./SimulatorParams.h"
#include "./Constraints.h"
#include "../file_format/yarnRepr.h"
#include "./macros.h"
#include "./ParallelWorker.h"
#include "./BaseSimulator.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <mutex>
#include <thread>

namespace simulator{

// Simulator contains all yarn simulation functionalities
//
// A yarn is modeled by a Catmull-Rom spline with #m control points
// denoted by q[i], i from 0 to m - 1. Each segment is governed by
// 4 control points, so there are #N = #m - 3 segments (0-indexed).
// segment[i] is governed by q[i], q[i + 1], q[i + 2], q[i + 3].
class Simulator : public BaseSimulator {
private:
  // simulator attributes
  size_t m;
  int stepCnt;
  std::thread simulatorThread;
  mutable std::mutex dataLock;

  // position of control points
  Eigen::MatrixXf q;

  // segment length
  std::vector<float> segmentLength;

  // mass matrix
  Eigen::SparseMatrix<float> M;
  Eigen::SparseMatrix<float> MInverse;

  // first-derivative of q w.r.t. time
  Eigen::MatrixXf qD;

  // gradient of positional energy
  Eigen::MatrixXf gradE;

  // gradient of damping energy
  Eigen::MatrixXf gradD;

  // external force
  Eigen::MatrixXf f;

  mutable std::mutex gradLock;

  // constraints
  Constraints constraints;

  // debug
  void writeToFile() const;

  void calculateSegmentLength();
  void addSegmentLengthConstraint();
  void addControlPointLengthConstraint(int i, int j);

  void constructMassMatrix();

  void calculateGradient(const std::function<bool()> cancelled);
  void calculateBendingEnergyGradient(int i);
  void calculateLengthEnergyGradient(int i);
  void calculateCollisionGradient(int i, int j);

  void calculateGlobalDampingGradient(int i);

  void fastProjection();

public:
  // Empty constructor
  Simulator() : q(0, 1), constraints(0) {};

  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  Simulator(file_format::YarnRepr yarns, SimulatorParams params_);

  ~Simulator();

  // Returns current yarns
  virtual const file_format::YarnRepr& getYarns() override {
    return yarns;
  }

  void step(const std::function<bool()>& cancelled) override;

  // Gets a reference to constraint container
  Constraints& getConstraints() { return constraints; }
};

}; // namespace simulator
