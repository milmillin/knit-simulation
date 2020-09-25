#include "./Simulator.h"
#include "./Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Bending Energy
//

void Simulator::calculateBendingEnergy(int thread_id, size_t i) {
  EIGEN_UNUSED_VARIABLE(thread_id)

  int index = i * 3;
  DECLARE_POINTS2(p, Q, index);

  double coefficient = params.kBend * segmentLength;

  using Vec12 = Eigen::Matrix<double, 12, 1>;

  Vec12 res = integrate<Vec12>([&p](double s)->Vec12 {
    Vec12 ans;

    DECLARE_BASIS_D2(bD, s);
    DECLARE_BASIS_DD2(bDD, s);

    Eigen::Vector3d PD = POINT_FROM_BASIS(p, bD);
    Eigen::Vector3d PDD = POINT_FROM_BASIS(p, bDD);

    Eigen::Vector3d cross = PD.cross(PDD);

    double PDnorm2 = PD.squaredNorm();
    double PDnorm6 = PDnorm2 * PDnorm2 * PDnorm2;
    double PDnorm8 = PDnorm6 * PDnorm2;
    double crossNorm2 = cross.squaredNorm();

    for (int kk = 0; kk < 4; kk++) {
      ans(kk * 3ll) = 2.0 * (cross(1) * (PD(2) * bDD[kk] - PDD(2) * bD[kk]) + cross(2) * (PDD(1) * bD[kk] - PD(1) * bDD[kk])) / PDnorm6;
      ans(kk * 3ll + 1) = 2.0 * (cross(0) * (PDD(2) * bD[kk] - PD(2) * bDD[kk]) + cross(2) * (PD(0) * bDD[kk] - PDD(0) * bD[kk])) / PDnorm6;
      ans(kk * 3ll + 2) = 2.0 * (cross(0) * (PD(1) * bDD[kk] - PDD(1) * bD[kk]) + cross(1) * (PDD(0) * bD[kk] - PD(0) * bDD[kk])) / PDnorm6;

      ans.block<3, 1>(kk * 3ll, 0) += (-6.0 * bD[kk] * crossNorm2 / PDnorm8) * PD;
    }

    return ans;
    }, 0, 1);

  {
    std::lock_guard<std::mutex> lock(lockF);
    F.block<12, 1>(index, 0) -= coefficient * res;
  }
}

}  // namespace simulator
