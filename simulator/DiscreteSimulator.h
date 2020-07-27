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
  void step() override;


  // Force
  void applyGravity();

  void applyContactForce();

  void initBendingForceMetadata();
  void applyBendingForce();
  void applyTwistingForce();
  void updateBendingForceMetadata();

  // Constrain
  void applyLengthConstrain();
  void applyPinConstrain();

  void fastProjection();

  // Velocity filter
  void applyGroundVelocityFilter();
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

  inline int pointIndex(int pointID, int axis);
  inline int nConstrain();

  aabb::Tree collisionTree;
};

} // namespace simulatr 

