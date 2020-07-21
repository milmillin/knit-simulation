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
  // Empty constructor
  DiscreteSimulator();

  // Constructs a new simulator with control points
  //
  // q_ : The #m x 3 matrix containing initial control points.
  // params_ : Simulation paramters
  DiscreteSimulator(file_format::YarnRepr yarns, SimulatorParams params);

  // Returns current yarns
  virtual const file_format::YarnRepr &getYarns() override { return this->yarns; };

  // Simulates next timestep.
  void step(const std::function<bool()>& cancelled) override;

  void applyGravity();
  void applyGroundVelocityFilter();
  void applyContactForce();
  void applyLengthConstrain();
  void applyPinConstrain();
  void fastProjection();
  void applyGlobalDamping();
 private:
  // Velocity
  Eigen::MatrixXf dQ;
  // Acceleration
  Eigen::MatrixXf ddQ;

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

  aabb::Tree collisionTree;
};

} // namespace simulatr 

