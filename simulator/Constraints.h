#pragma once

#include <Eigen/Core>
#include <vector>
#include <functional>

namespace simulator {

class Constraints {
public:
  struct Entry {
    int index;
    std::function<float(const Eigen::MatrixXf&)> f;
  };

  // Construct a constraint container with m control points.
  // To be called by Simulator only.
  Constraints(size_t m_) : m(m_) { }

  // Adds constraint
  // f: constrain function of q
  // fD: d(f)/d(q(index)) for every index with non-zero value
  void addConstraint(std::function<float(const Eigen::MatrixXf&)> f,
    std::vector<Entry> fD);

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
  float calculateMax(const Eigen::MatrixXf& q) const;

  // Evaluates the constraints at q
  // Returns the matrix of size c x 1
  Eigen::MatrixXf calculate(const Eigen::MatrixXf& q) const;

private:
  // number of control points
  size_t m;
  
  // list of constraints
  std::vector<std::function<float(const Eigen::MatrixXf&)>> C;

  // list of jacobian derivatives
  std::vector<std::vector<Entry>> CD;

};

} // namespace simulator