#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Bending Energy
//

void Simulator::calculateBendingEnergy(int i) {
  int index = i * 3;
  DECLARE_POINTS2(p, Q, index);

  float coefficient = params.kBend * catmullRomLength[i];

  using Vec12 = Eigen::Matrix<float, 12, 1>;

  Vec12 res = integrate<Vec12>([&p](float s)->Vec12 {
    Vec12 ans;

    DECLARE_BASIS_D2(bD, s);
    DECLARE_BASIS_DD2(bDD, s);

    Eigen::Vector3f PD = POINT_FROM_BASIS(p, bD);
    Eigen::Vector3f PDD = POINT_FROM_BASIS(p, bDD);

    Eigen::Vector3f cross = PD.cross(PDD);

    float PDnorm2 = PD.squaredNorm();
    float PDnorm6 = PDnorm2 * PDnorm2 * PDnorm2;
    float PDnorm8 = PDnorm6 * PDnorm2;
    float crossNorm2 = cross.squaredNorm();

    for (int kk = 0; kk < 4; kk++) {
      ans(kk * 3ll) = 2.f * (cross(1) * (PD(2) * bDD[kk] - PDD(2) * bD[kk]) + cross(2) * (PDD(1) * bD[kk] - PD(1) * bDD[kk])) / PDnorm6;
      ans(kk * 3ll + 1) = 2.f * (cross(0) * (PDD(2) * bD[kk] - PD(2) * bDD[kk]) + cross(2) * (PD(0) * bDD[kk] - PDD(0) * bD[kk])) / PDnorm6;
      ans(kk * 3ll + 2) = 2.f * (cross(0) * (PD(1) * bDD[kk] - PDD(1) * bD[kk]) + cross(1) * (PDD(0) * bD[kk] - PD(0) * bDD[kk])) / PDnorm6;

      ans.block<3, 1>(kk * 3ll, 0) += (-6.f * bD[kk] * crossNorm2 / PDnorm8) * PD;
    }

    return ans;
    }, 0, 1);

  {
    std::lock_guard<std::mutex> lock(lockF);
    F.block<12, 1>(index, 0) -= coefficient * res;
  }
}

}  // namespace simulator
