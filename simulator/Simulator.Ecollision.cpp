#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Collision Energy
//

static constexpr int SUBDIVIDE = 20;

static inline Eigen::Vector3f fff(Eigen::Vector3f gradient, float norm2, float thresh2) {
  return (-thresh2 / norm2 / norm2 + 1 / thresh2) * gradient;
}

void Simulator::calculateCollisionEnergyGradient(int ii, int jj) {
  int iIndex = ii * 3;
  int jIndex = jj * 3;

  attrLock.lock();
  DECLARE_POINTS2(pi, iIndex);
  DECLARE_POINTS2(pj, jIndex);

  float r = yarns.yarns[0].radius;
  float coefficient = params.kContact * segmentLength[ii] * segmentLength[jj];
  attrLock.unlock();

  float thresh2 = 4.0f * r * r;

  // calculate only for j > i
  coefficient *= 2;

  float step = 1.f / SUBDIVIDE;
  float halfStep = step / 2;

  Eigen::MatrixXf gradEi = Eigen::MatrixXf::Zero(12, 1);
  Eigen::MatrixXf gradEj = Eigen::MatrixXf::Zero(12, 1);

  for (int i = 0; i < SUBDIVIDE; i++) {
    float si = i * step + halfStep;
    DECLARE_BASIS(bi, si);
    Eigen::Vector3f pi = bi1 * pi1 + bi2 * pi2 + bi3 * pi3 + bi4 * pi4;
    for (int j = 0; j < SUBDIVIDE; j++) {
      float sj = j * step + halfStep;
      DECLARE_BASIS(bj, sj);
      Eigen::Vector3f pj = bj1 * pj1 + bj2 * pj2 + bj3 * pj3 + bj4 * pj4;

      Eigen::Vector3f diff = pj - pi;

      float norm2 = diff.squaredNorm();
      if (norm2 >= thresh2) continue;

      gradEi.block<3, 1>(0, 0) += fff(-2 * bi1 * diff, norm2, thresh2);
      gradEi.block<3, 1>(3, 0) += fff(-2 * bi2 * diff, norm2, thresh2);
      gradEi.block<3, 1>(6, 0) += fff(-2 * bi3 * diff, norm2, thresh2);
      gradEi.block<3, 1>(9, 0) += fff(-2 * bi4 * diff, norm2, thresh2);

      gradEj.block<3, 1>(0, 0) += fff(2 * bj1 * diff, norm2, thresh2);
      gradEj.block<3, 1>(3, 0) += fff(2 * bj2 * diff, norm2, thresh2);
      gradEj.block<3, 1>(6, 0) += fff(2 * bj3 * diff, norm2, thresh2);
      gradEj.block<3, 1>(9, 0) += fff(2 * bj4 * diff, norm2, thresh2);
    }
  }
  {
    std::lock_guard<std::mutex> lock(gradELock);
    gradE.block<12, 1>(iIndex, 0) += coefficient * gradEi * step * step;
    gradE.block<12, 1>(jIndex, 0) += coefficient * gradEj * step * step;
  }
}

}  // namespace simulator
