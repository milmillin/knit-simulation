#pragma once


#include "./threading/ctpl_stl.h"
#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"
#include "./AABB.h"
#include "Constraints.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <functional>
#include <mutex>

namespace simulator
{

#define IFDEBUG if (params.debug)

class BaseSimulator {
public:
  using StateGetter = std::function<bool()>;

  virtual ~BaseSimulator() { }

  // Returns current yarns
  const file_format::YarnRepr& getYarns();

  file_format::YarnRepr getVelocityYarns();

  // Set velocity
  void setPosition(const file_format::YarnRepr& yarn);
  void setVelocity(const file_format::YarnRepr& yarn);

  // Simulates next timestep.
  void step(const StateGetter& cancelled);

  // TODO: Locks
  SimulatorParams getParams() const { return params; }
protected:
  ctpl::thread_pool thread_pool;

  // Position and meta-data
  file_format::YarnRepr yarns;

  // Simulation parameters
  SimulatorParams params;

  // Number of control points
  int m;

  // Flattened Coordinates [x0 y0 z0 x1 y1 x1 ...]
  // Position of control points
  Eigen::MatrixXf Q;
  // Velocity
  Eigen::MatrixXf dQ;
  // Force
  Eigen::MatrixXf F;
  mutable std::mutex lockF;

  // Acceleration can be derived from ddQ = invM * F
  // F = f - gradE - gradD

  // Mass Matrix
  Eigen::SparseMatrix<float> M;
  // Inverse
  Eigen::SparseMatrix<float> invM;

  // Constraints
  Constraints constraints;

  // Length for each segment
  std::vector<float> segmentLength;
  std::vector<float> catmullRomLength;

  // Collision Tree
  aabb::Tree collisionTree;

  ///////////////////////
  // Constructor

  BaseSimulator(file_format::YarnRepr _yarns, SimulatorParams _params);

  // Initialize mass matrix
  virtual void constructMassMatrix();

  // Add desired constraints
  virtual void setUpConstraints() = 0;

  // Call virtual initialization methods.
  // Must be called in derived class constructor.
  void initialize();

  ///////////////////////
  // Stepping

  // Update F, dQ; Do unconstrained step
  virtual void stepImpl(const StateGetter& cancelled) = 0;
  void fastProjection(const StateGetter& cancelled);
  void updateCollisionTree(const StateGetter& cancelled);
  // Called after all update is done
  virtual void postStep(const StateGetter& cancelled) {};


  ///////////////////////
  // Contact Force

  void applyContactForce(const StateGetter& cancelled);
  void contactForceBetweenSegments
      (int thread_id,
      std::vector<Eigen::MatrixXf> *forces,
      int ii, int jj);

  ///////////////////////
  // Constraints

  // Add length constraint of line segment defined by points i and i + 1.
  void addSegmentLengthConstraint(int i);
  // Add length constraint of Catmull-Rom segment defined by points i to i + 3.
  void addCatmullRomLengthConstraint(int i);
  // Add pin constraint of point i to a fixed position
  void addPinConstraint(int i, Eigen::Vector3f position);
};

} // namespace simulator
