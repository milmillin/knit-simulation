#include "Constraints.h"

namespace simulator {

void Constraints::addConstraint(Func f, JacobianFunc fD) {
  constraints.push_back(Entry{ f, fD });
}

void Constraints::addPinConstrain(int i, Eigen::Vector3f point) {
  int index = i * 3;

  for (int kk = 0; kk < 3; kk++) {
    Func f = [=](const Eigen::MatrixXf& q)->float {
      return q(index + kk) - point(kk);
    };

    JacobianFunc fD = [=](const Eigen::MatrixXf& q, const Referrer& ref) {
      ref(index + kk) += 1;
    };

    constraints.push_back(Entry{ f, fD });
  }
}

void Constraints::addGlueConstrain(int i, int j) {
  int iIndex = i * 3;
  int jIndex = j * 3;

  for (int kk = 0; kk < 3; kk++) {
    Func f = [=](const Eigen::MatrixXf& q)->float {
      return q(iIndex + kk) - q(jIndex + kk);
    };

    JacobianFunc fD = [=](const Eigen::MatrixXf& q, const Referrer& ref) {
      ref(iIndex + kk) += 1;
      ref(jIndex + kk) += -1;
    };

    constraints.push_back(Entry{ f, fD });
  }
}

void Constraints::addLengthConstrain(int i, float length) {
  int index = i * 3;

  Func f = [=](const Eigen::MatrixXf& q)->float {
    DECLARE_POINTS2(p, q, index);
    float currentLength = integrate<float>([&](float s)->float {
      DECLARE_BASIS_D2(bD, s);
      return POINT_FROM_BASIS(p, bD).norm();
      }, 0, 1);
    return 1 - currentLength / length;
  };

  using Vec12 = Eigen::Matrix<float, 12, 1>;

  JacobianFunc fD = [=](const Eigen::MatrixXf& q, const Referrer& ref) {
    DECLARE_POINTS2(p, q, index);
    Vec12 res = integrate<Vec12>([&](float s)->Vec12 {
      Vec12 ans;
      DECLARE_BASIS_D2(bD, s);
      Eigen::Vector3f P = POINT_FROM_BASIS(p, bD);
      float norm = P.norm();

      for (int kk = 0; kk < 4; kk++) {
        ans.block<3, 1>(3 * kk, 0) = (bD[kk] / norm) * P;
      }
      return ans;
      }, 0, 1);

    res *= -1.f / length;

    for (int ii = 0; ii < 12; ii++) {
      ref(index + ii) += res(ii);
    }
  };

  constraints.push_back(Entry{ f, fD });
}

Eigen::SparseMatrix<float> Constraints::getJacobian(const Eigen::MatrixXf& q) const {
  int c = (int)constraints.size();
	Eigen::SparseMatrix<float> J = Eigen::SparseMatrix<float>(c, 3 * m);
  int i;
  Referrer ref = [&J, &i](int index)->float& {return J.coeffRef(i, index); };

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