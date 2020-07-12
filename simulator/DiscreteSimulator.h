#include <vector>
#include <Eigen/Core>

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
  SimulatorParams params;
 private:
  Eigen::MatrixXf dQ;
  Eigen::MatrixXf ddQ;
  file_format::YarnRepr yarns;
  std::vector<float> restLength;
};

} // namespace simulatr 

