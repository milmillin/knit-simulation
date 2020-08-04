#include "Constraints.h"

#include <mutex>

#include "easy_profiler_stub.h"
#include "./threading/threading.h"

namespace simulator {

void Constraints::addConstraint(Func f, JacobianFunc fD) {
  constraints.push_back(Entry{ f, fD });
}

Eigen::SparseMatrix<double> Constraints::getJacobian(const Eigen::MatrixXd& q) const {
  EASY_FUNCTION();
  int c = (int)constraints.size();
  std::vector<Eigen::SparseMatrix<double>> Js(thread_pool->size(), Eigen::SparseMatrix<double>(c, 3 * m));

  auto getJacobianTask = [this, &Js, &q](int thread_id, int start_index, int end_index) {
    int i;
    Referrer ref = [&Js, &i, thread_id](int index, int ax)->double& {
      return Js[thread_id].coeffRef(i, index * 3ll + ax); 
    };
    for (i = start_index; i < end_index; i++) {
      constraints[i].fD(q, ref);
    }
  };

  threading::runSequentialJob(*thread_pool, getJacobianTask, 0, (int)constraints.size());

  Eigen::SparseMatrix<double> J(c, 3 * m);
  for (const auto& j : Js) {
    J += j;
  }

  return J;
}

Eigen::MatrixXd Constraints::calculate(const Eigen::MatrixXd& q) const {
  EASY_FUNCTION();
  Eigen::MatrixXd res(constraints.size(), 1);

  auto calculateTask = [this, &res, &q](int thread_id, int start_index, int end_index) {
    for (int i = start_index; i < end_index; i++) {
      res(i) = constraints[i].f(q);
    }
  };

  threading::runSequentialJob(*thread_pool, calculateTask, 0, (int)constraints.size());

  return res;
}

} // namespace simulator