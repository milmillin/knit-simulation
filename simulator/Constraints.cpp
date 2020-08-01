#include "Constraints.h"

#include "../easy_profiler_stub.h"
#include "./threading/threading.h"

#include <mutex>

namespace simulator {

void Constraints::addConstraint(Func f, JacobianFunc fD) {
  constraints.push_back(Entry{ f, fD });
}

Eigen::SparseMatrix<float> Constraints::getJacobian(const Eigen::MatrixXf& q) const {
  EASY_FUNCTION();
  int c = (int)constraints.size();
	Eigen::SparseMatrix<float> J = Eigen::SparseMatrix<float>(c, 3 * m);
  std::mutex Jlock;

  auto getJacobianTask = [this, &J, &Jlock, &q](int thread_id, int start_index, int end_index) {
    int i;
    Referrer ref = [&J, &Jlock, &i](int index, int ax)->float& {
      std::lock_guard<std::mutex> lock(Jlock);
      return J.coeffRef(i, index * 3ll + ax); 
    };
    for (i = start_index; i < end_index; i++) {
      constraints[i].fD(q, ref);
    }
  };

  threading::runSequentialJob(*thread_pool, getJacobianTask, 0, (int)constraints.size());

  return J;
}

Eigen::MatrixXf Constraints::calculate(const Eigen::MatrixXf& q) const {
  EASY_FUNCTION();
  Eigen::MatrixXf res(constraints.size(), 1);

  auto calculateTask = [this, &res, &q](int thread_id, int start_index, int end_index) {
    for (int i = start_index; i < end_index; i++) {
      res(i) = constraints[i].f(q);
    }
  };

  threading::runSequentialJob(*thread_pool, calculateTask, 0, (int)constraints.size());

  return res;
}

} // namespace simulator