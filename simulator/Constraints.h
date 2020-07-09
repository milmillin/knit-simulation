#pragma once

#include <Eigen/Core>
#include <vector>
#include <functional>

namespace simulator {

class Constraints {
private:
  struct Entry {
    int index;
    std::function<float(const Eigen::MatrixXf&)> f;
  };

  // number of control points
  size_t m;
  
  // list of constraints
  std::vector<std::function<float(const Eigen::MatrixXf&)>> C;

  // list of jacobian derivatives
  std::vector<std::vector<Entry>> CD;

public:
  // Construct a constraint container with m control points.
  // To be called by Simulator only.
  Constraints(size_t m_) : m(m_) { }

  // Pins i-th control point to <x, y, z>.
  // i in [0, m - 1].
  // Adds 3 contraints.
  void addPinConstrain(int i, float x, float y, float z);

  // Glues i-th and j-th control points together
  // i, j in [0, m - 1]
  // Adds 3 constraints.
  void addGlueConstrain(int i, int j);

  // Gets Jacobian matrix of size c x (3 * m)
  // where c is the number of constraints
  // evaluated at q.
  Eigen::MatrixXf getJacobian(const Eigen::MatrixXf& q) const;

  // Evaluates the constraints at q
  // Returns the maximum of absolute of each constraint
  float calculate(const Eigen::MatrixXf& q) const;
};

} // namespace simulator