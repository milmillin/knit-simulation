#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"

namespace simulator {

class DiscreteSimulator {
 public:
  // Empty constructor
  DiscreteSimulator();

  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params);

  // Returns current yarns
  const file_format::YarnRepr &getYarns() const { return this->yarns; };

  // Simulates next timestep.
  void step();

  void applyGravity();
  void applyGroundVelocityFilter();
  void applyContactForce();
  void applyLengthConstrain();
  void applyPinConstrain();
  void fastProjection();
  void applyGlobalDamping();
  SimulatorParams params;
 private:
  // Velocity
  Eigen::MatrixXf dQ;
  // Acceleration
  Eigen::MatrixXf ddQ;
  // Position and meta-data
  file_format::YarnRepr yarns;
  // Rest Length for each segment
  std::vector<float> restLength;

  // === Constrains ===
  // Index of control points to be fixed
  std::vector<int> pinControlPoints;
  // Location of control points to be fixed
  std::vector<glm::vec3> pinControlPointsPosition;
  // Value of constrain function (one row per constrain)
  Eigen::VectorXf constrain;
  // Gradient of constrain function (one row per constrain and one colomn per coordiate)
  // The i^th axis of j^th point is stored in the colomn `i * nPoints + j`
  Eigen::SparseMatrix<float> dConstrain;
  // Next available constrain id
  int nextConstrainID = 0;

  inline int pointIndex(int pointID, int axis);
  inline int nConstrain();
};

} // namespace simulatr 

