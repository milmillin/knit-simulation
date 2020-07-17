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

void Simulator::calculateCollisionGradient(int ii, int jj) {
  int iIndex = ii * 3;
  int jIndex = jj * 3;

  dataLock.lock();
  DECLARE_POINTS2(pi, q, iIndex);
  DECLARE_POINTS2(pj, q, jIndex);
  DECLARE_POINTS2(piD, qD, iIndex);
  DECLARE_POINTS2(pjD, qD, jIndex);

  float r = yarns.yarns[0].radius;
  float coeffE = params.kContact * segmentLength[ii] * segmentLength[jj];
  float coeffD = segmentLength[ii] * segmentLength[jj];
  float kDt = params.kDt;
  float kDn = params.kDn;
  dataLock.unlock();

  float thresh2 = 4.0f * r * r;

  // calculate only for j > i
  coeffE *= 2;
  coeffD *= 2;

  float step = 1.f / SUBDIVIDE;
  float halfStep = step / 2;

  Eigen::MatrixXf gradEi = Eigen::MatrixXf::Zero(12, 1);
  Eigen::MatrixXf gradEj = Eigen::MatrixXf::Zero(12, 1);

  Eigen::MatrixXf gradDi = Eigen::MatrixXf::Zero(12, 1);
  Eigen::MatrixXf gradDj = Eigen::MatrixXf::Zero(12, 1);

  for (int i = 0; i < SUBDIVIDE; i++) {
    float si = i * step + halfStep;
    DECLARE_BASIS2(bi, si);
    Eigen::Vector3f Pi = POINT_FROM_BASIS(pi, bi);
    Eigen::Vector3f PiD = POINT_FROM_BASIS(piD, bi);

    for (int j = 0; j < SUBDIVIDE; j++) {
      float sj = j * step + halfStep;
      DECLARE_BASIS2(bj, sj);

      Eigen::Vector3f Pj = POINT_FROM_BASIS(pj, bj);
      Eigen::Vector3f diff = Pj - Pi;

      float norm2 = diff.squaredNorm();
      if (norm2 >= thresh2) continue;

      Eigen::Vector3f PjD = POINT_FROM_BASIS(pjD, bj);
      Eigen::Vector3f diffD = PjD - PiD;

      float tmp = -2 * (kDt - kDn) / sqrt(norm2);

      for (int kk = 0; kk < 4; kk++) {
        // Contact Energy
        gradEi.block<3, 1>(3 * kk, 0) += fff(-2 * bi[kk] * diff, norm2, thresh2);
        gradEj.block<3, 1>(3 * kk, 0) += fff(2 * bj[kk] * diff, norm2, thresh2);

        // Contact Damping
        gradDi.block<3, 1>(3 * kk, 0) += kDt * (-2 * bi[kk] * diffD) + tmp * (-bi[kk] * diff);
        gradDj.block<3, 1>(3 * kk, 0) += kDt * (2 * bi[kk] * diffD) + tmp * (bi[kk] * diff);
      }
    }
  }
  {
    std::lock_guard<std::mutex> lock(gradLock);
    gradE.block<12, 1>(iIndex, 0) += coeffE * gradEi * step * step;
    gradE.block<12, 1>(jIndex, 0) += coeffE * gradEj * step * step;

    gradD.block<12, 1>(iIndex, 0) += coeffD * gradDi * step * step;
    gradD.block<12, 1>(jIndex, 0) += coeffD * gradDj * step * step;
  }
}

}  // namespace simulator
