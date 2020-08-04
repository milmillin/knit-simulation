#include "./Simulator.h"
#include "./Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Global Damping
//

void Simulator::calculateGlobalDamping(int i) {
  int index = i * 3;

  DECLARE_POINTS2(pD, dQ, index);

  double coefficient = params.kGlobal;

  using Vec12 = Eigen::Matrix<double, 12, 1>;

  Vec12 res = integrate<Vec12>([&pD](double s)->Vec12 {
    Vec12 ans;
    DECLARE_BASIS2(b, s);

    Eigen::Vector3d P = POINT_FROM_BASIS(pD, b);

    for (int kk = 0; kk < 4; kk++) {
      ans.block<3, 1>(kk * 3ll, 0) = (2 * b[kk]) * P;
    }
    return ans;
    }, 0, 1);

  {
    std::lock_guard<std::mutex> lock(lockF);
    F.block<12, 1>(index, 0) -= coefficient * res;
  }
}


} // namespace simulator
