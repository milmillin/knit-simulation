#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Jacobi>
#include "spdlog/spdlog.h"

#include "./Helper.h"

namespace simulator {

void symmeticSchur(const Eigen::Matrix3d &A, int p, int q, double *c, double *s) {
  if (std::abs(A(p, q)) > 1e-6) {
    double tau = (A(q, q) - A(p, p)) / (2 * A(p, q));
    double t;
    if (tau >= 0.0) {
      t = 1.0 / (tau + std::sqrt(1 + tau * tau));
    } else {
      t = -1.0 / (-tau + std::sqrt(1 + tau * tau));
    }
    *c = 1.0 / std::sqrt(1 + t * t);
    *s = t * (*c);
  } else {
    *c = 1.0;
    *s = 0.0;
  }
}


// Run cyclic Jacobi method to diagonize matrix `A`.
// Suppose `A=X` at the beginning. When this function returns,
// `V * A * V.transpose() = X` where `A` is a nearly diagonized matrix.
//
// A: the matrix that will be diagonized
// V: the eigen values
// iterations: number of sweeps to run
void cyclicJacobi(Eigen::Matrix3d *A, Eigen::Matrix3d *V, int iterations) {
  for (int k = 0; k < iterations; k++) {
    for (int p = 0; p < A->rows(); p++) {
      for (int q = p + 1; q < A->rows(); q++) {
        double c, s;
        symmeticSchur(*A, p, q, &c, &s);
        Eigen::JacobiRotation<double> J(c, s);
        A->applyOnTheLeft(p, q, J.transpose());
        A->applyOnTheRight(p, q, J);
        V->applyOnTheRight(p, q, J);
      }
    }
  }
}

// Suppose that `V * A * V.transpose() = X` and `A` is nearly diagonized,
// this function sets $A = X^{-1/2}$
void inverseSquareRoot(Eigen::Matrix3d *A, const Eigen::Matrix3d &V) {
  for (int k = 0; k < A->rows(); k++) {
    (*A)(k, k) = 1.0 / std::sqrt(std::abs((*A)(k, k)));
  }
  *A = V * (*A) * V.transpose();
}

}  // namespace simulator
