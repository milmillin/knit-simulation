#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Length Energy
//

void Simulator::calculateLengthEnergyGradient(int i) {
  int index = i * 3;
  DECLARE_POINTS2(p, q, index);

  float coefficient = params.kLen;
  float L = segmentLength[i];

  using Vec12 = Eigen::Matrix<float, 12, 1>;

  Vec12 res = integrate<Vec12>([L, &p](float s)->Vec12 {
    Vec12 ans;
    DECLARE_BASIS_D2(bD, s);
    Eigen::Vector3f P = POINT_FROM_BASIS(p, bD);
    float norm = P.norm();

    float tmp = 2.f * (norm / L - 1) / (L * norm);

    for (int kk = 0; kk < 4; kk++) {
      ans.block<3, 1>(kk * 3, 0) = (tmp * bD[kk]) * P;
    }
    return ans;
    }, 0, 1);

  gradE.block<12, 1>(index, 0) += coefficient * res;
}


}  // namespace simulator
