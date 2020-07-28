#include "Constraints.h"

namespace simulator {

void Constraints::addConstraint(Func f, JacobianFunc fD) {
  constraints.push_back(Entry{ f, fD });
}

Eigen::SparseMatrix<float> Constraints::getJacobian(const Eigen::MatrixXf& q) const {
  int c = (int)constraints.size();
	Eigen::SparseMatrix<float> J = Eigen::SparseMatrix<float>(c, 3 * m);
  int i;
  Referrer ref = [&J, &i](int index, int ax)->float& {return J.coeffRef(i, index * 3ll + ax); };

  for (i = 0; i < c; i++) {
    constraints[i].fD(q, ref);
  }

  return J;
}

Eigen::MatrixXf Constraints::calculate(const Eigen::MatrixXf& q) const {
  Eigen::MatrixXf res(constraints.size(), 1);
  for (int i = 0; i < (int)constraints.size(); i++) {
    res(i) = constraints[i].f(q);
  }
  return res;
}

} // namespace simulator