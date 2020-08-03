#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Length Energy
//

void Simulator::calculateLengthEnergy(int i) {
  int index = i * 3;
  DECLARE_POINTS2(p, Q, index);

  double coefficient = params.kLen;
  double L = catmullRomLength[i];

  using Vec12 = Eigen::Matrix<double, 12, 1>;

  Vec12 res = integrate<Vec12>([L, &p](double s)->Vec12 {
    Vec12 ans;
    DECLARE_BASIS_D2(bD, s);
    Eigen::Vector3d P = POINT_FROM_BASIS(p, bD);
    double norm = P.norm();

    double tmp = 2.f * (norm / L - 1) / (L * norm);

    for (int kk = 0; kk < 4; kk++) {
      ans.block<3, 1>(kk * 3ll, 0) = (tmp * bD[kk]) * P;
    }
    return ans;
    }, 0, 1);

  {
    std::lock_guard<std::mutex> lock(lockF);
    F.block<12, 1>(index, 0) -= coefficient * res;
  }
}

}  // namespace simulator
