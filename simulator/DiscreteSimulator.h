#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"

#include <Eigen/Core>

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
  SimulatorParams params;
 private:
  Eigen::MatrixXf dQ;
  Eigen::MatrixXf ddQ;
  file_format::YarnRepr yarns;
};

} // namespace simulatr 

