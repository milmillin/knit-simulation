#include "Constraints.h"

namespace simulator {

void Constraints::addConstraint(std::function<float(const Eigen::MatrixXf&)> f,
  std::vector<Entry> fd) {
  C.push_back(f);
  CD.push_back(fd);
}

void Constraints::addPinConstrain(int i, float x, float y, float z) {
  // TODO: not sure if this really works

  int index = i * 3;

  // px - x
  C.push_back([=](const Eigen::MatrixXf& q)->float {
    return q(index) - x;
    });

  // d(px - x)/d(px) = 1
  CD.push_back(std::vector<Entry>());
  CD.back().push_back(Entry{ index, [=](const Eigen::MatrixXf& q)->float {
    return 1;
} });

  // py - y
  C.push_back([=](const Eigen::MatrixXf& q)->float {
    return q(index + 1) - y;
    });

  // d(py - y)/d(py) = 1
  CD.push_back(std::vector<Entry>());
  CD.back().push_back(Entry{ index + 1, [=](const Eigen::MatrixXf& q)->float {
    return 1;
} });

  // pz - z
  C.push_back([=](const Eigen::MatrixXf& q)->float {
    return q(index + 2) - z;
    });

  // d(pz - z)/d(py) = 1
  CD.push_back(std::vector<Entry>());
  CD.back().push_back(Entry{ index + 2, [=](const Eigen::MatrixXf& q)->float {
    return 1;
} });
}

void Constraints::addGlueConstrain(int i, int j) {
  int iIndex = i * 3;
  int jIndex = j * 3;

  for (int kk = 0; kk < 3; kk++) {
    // pix - pjx
    C.push_back([=](const Eigen::MatrixXf& q)->float {
      return q(iIndex + kk) - q(jIndex + kk);
      });

    CD.push_back(std::vector<Entry>());
    // d(pix - pjx)/d(pix) = 1
    CD.back().push_back(Entry{ iIndex + kk, [](const Eigen::MatrixXf& q)->float {
      return 1;
      } });
    // d(pix - pjx)/d(pjx) = -1
    CD.back().push_back(Entry{ jIndex + kk, [](const Eigen::MatrixXf& q)->float {
      return -1;
      } });
  }
}

Eigen::MatrixXf Constraints::getJacobian(const Eigen::MatrixXf& q) const {
  int c = (int)C.size();
	Eigen::MatrixXf J = Eigen::MatrixXf::Zero(c, 3 * m);

  for (int i = 0; i < c; i++) {
    for (int j = 0; j < CD[i].size(); j++) {
      J(i, CD[i][j].index) += CD[i][j].f(q);
    }
  }

  return J;
}

float Constraints::calculateMax(const Eigen::MatrixXf& q) const {
  float res = 0;
  for (auto& constraint : C) {
    res = std::max(res, fabs(constraint(q)));
  }
  return res;
}

Eigen::MatrixXf Constraints::calculate(const Eigen::MatrixXf& q) const {
  Eigen::MatrixXf res(C.size(), 1);
  for (int i = 0; i < (int)C.size(); i++) {
    res(i) = C[i](q);
  }
  return res;
}

} // namespace simulator