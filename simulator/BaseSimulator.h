#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "spdlog/spdlog.h"

#include <functional>
#include <mutex>
#include <atomic>

#include "file_format/yarnRepr.h"
#include "./threading/ctpl_stl.h"
#include "./SimulatorParams.h"
#include "./AABB.h"
#include "./Constraints.h"
#include "./CacheArray2D.h"

namespace simulator
{

class BaseSimulator {
public:
  using StateGetter = std::function<bool()>;

  virtual ~BaseSimulator() { }

  // Returns current yarns
  virtual const file_format::YarnRepr& getYarns();

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
  Eigen::MatrixXd Q;
  // Velocity
  Eigen::MatrixXd dQ;
  // Force
  Eigen::MatrixXd F;
  mutable std::mutex lockF;

  // Acceleration can be derived from ddQ = invM * F
  // F = f - gradE - gradD

  // Mass Matrix
  Eigen::SparseMatrix<double> M;
  // Inverse
  Eigen::SparseMatrix<double> invM;

  // Constraints
  Constraints constraints;

  // Length for each segment
  std::vector<double> segmentLength;
  std::vector<double> catmullRomLength;

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

  typedef Eigen::Matrix<double, 12, 1> ControlPoints;
  typedef Eigen::Matrix<double, 12, 12> StiffnessMatrix;

  struct LinearizedContactModel {
    // Original segment location
    ControlPoints baseQI;
    ControlPoints baseQJ;

    // Original location relative to the center
    ControlPoints RelativeQI;
    ControlPoints RelativeQJ;

    // Original force
    ControlPoints baseForceI;
    ControlPoints baseForceJ;

    // Jecobian of force field.
    // `dForceIJ` means dF_i/dQ_j
    StiffnessMatrix dForceII;
    StiffnessMatrix dForceIJ;
    StiffnessMatrix dForceJI;
    StiffnessMatrix dForceJJ;

    // number of timesteps since last update
    int lastUpdate;

    // Is the model valid
    bool valid = false;
  };

  CacheArray2D<LinearizedContactModel> contactModelCache;

  // Calculate the contact force between segment `i` and `j`
  // Save the result in `forceI` and `forceJ`
  void contactForce(int i, int j, ControlPoints *forceI, ControlPoints *forceJ);
  // Build the linearized contact force model between segment `i` and `j`
  // Save the result in `model`
  void buildLinearModel(int i, int j, LinearizedContactModel *model);
  // Apply the force between segment `i` and `j` to the accumulator `forces`,
  // using the linearized model in `model`
  // Return `true` iff the approximation is valid
  bool applyApproxContactForce(int i, int j, Eigen::MatrixXd &forces, LinearizedContactModel *model);

  // catmullRomCoefficient(i, j) is the coefficient of the j^th control point
  // when the curve paramter s = (i + 0.5) / <number of samples>
  Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor>
    catmullRomCoefficient;
  // Initialize `catmullRomCoefficient`
  void initializeContactForceMetaData();
  void applyContactForce(const StateGetter& cancelled);
  void applyContactForceBetweenSegments
      (int thread_id,
      std::vector<Eigen::MatrixXd> *forces,
      int ii, int jj);

  ///////////////////////
  // Length spring

  // Use springs to constrain the segment length of each segments
  void applyLengthSpringForce();

  ///////////////////////
  // Constraints

  // Add length constraint of line segment defined by points i and i + 1.
  void addSegmentLengthConstraint(int i);
  // Add length constraint of Catmull-Rom segment defined by points i to i + 3.
  void addCatmullRomLengthConstraint(int i);
  // Add pin constraint of point i to a fixed position
  void addPinConstraint(int i, Eigen::Vector3d position);

  ///////////////////////
  // Statistics

  struct Statistics {
    std::atomic<int> totalContactCount;
    std::atomic<int> linearizedModelRebuildCount;
    std::atomic<int> approximationUsedCount;
    std::atomic<double> contactForceTotalError;
    std::atomic<int> contactForceErrorDivider;
    std::atomic<int> materialFrameUpdateCount;

   public:
    Statistics() {
      clear();
    }

    void clear() {
      totalContactCount = 0;
      linearizedModelRebuildCount = 0;
      approximationUsedCount = 0;
      contactForceTotalError = 0;
      contactForceErrorDivider = 0;
      materialFrameUpdateCount = 0;
    }
  } statistics;

  void printStatistics() {
    auto &stat = statistics;
    SPDLOG_INFO("* Total contacts: {}", stat.totalContactCount);
    SPDLOG_INFO("* Rebuild: {} ({}%)",
      stat.linearizedModelRebuildCount,
      100.0 * stat.linearizedModelRebuildCount / stat.totalContactCount);
    SPDLOG_INFO("* Approximation used: {} ({}%)",
      stat.approximationUsedCount,
      100.0 * stat.approximationUsedCount / stat.totalContactCount);
    SPDLOG_INFO("* Approximation error: {} ({})",
      stat.contactForceTotalError / stat.contactForceErrorDivider,
      stat.contactForceErrorDivider);
    SPDLOG_INFO("* Material frame update: {}", stat.materialFrameUpdateCount);
  }

  // === Debug helper ===

  // idea: only print the first occurance of NaN
  bool foundNaN = false;
  // If there's a NaN in `x`, print error message
  #define checkNaN(x) \
    if (x.hasNaN() && !foundNaN) { \
      SPDLOG_ERROR("NaN found!"); \
      foundNaN = true; \
    }
};

} // namespace simulator
