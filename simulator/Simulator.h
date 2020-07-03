#pragma once

#include "SimulatorParams.h"

#include <Eigen/Core>

namespace simulator {

// Simulator contains all yarn simulation functionalities
class Simulator {
  private:
    Eigen::MatrixXf q;
    SimulatorParams params;

    // first-derivative of q
    Eigen::MatrixXf qD;

    // gradient of positional energy
    Eigen::MatrixXf gradE;

    // gradient of damping energy
    Eigen::MatrixXf gradD;

    // external force
    Eigen::MatrixXf f;
  public:
    /**
     * Constructs a new Simulator
     *
     * @param q_ The #m x 3 matrix containing initial control points.
     * @param params_ Simulation paramters
     */
    Simulator(Eigen::MatrixXf q_, SimulatorParams params_);

    /**
     * @return the number of control points.
     */
    int getNumPoints() const;

    /**
     * @return #m x 3 matrix containing the current control points.
     */
    const Eigen::MatrixXf& getCurrentPoints();

    /**
     * Simulate the next timestep
     */
    void step();
};

}; // namespace simulator
