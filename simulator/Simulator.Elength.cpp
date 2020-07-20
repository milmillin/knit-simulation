#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Length Energy
//

void Simulator::calculateLengthEnergyGradient(int i) {
  int index = i * 3;
  DECLARE_POINTS(p, index);

  float coefficient = params.kLen;
  float L = segmentLength[i];

  // Bending Energy: px1
  gradE(index + 0) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD1 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: py1
  gradE(index + 1) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD1 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: pz1
  gradE(index + 2) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD1 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: px2
  gradE(index + 3) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD2 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: py2
  gradE(index + 4) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD2 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: pz2
  gradE(index + 5) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD2 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: px3
  gradE(index + 6) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD3 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: py3
  gradE(index + 7) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD3 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: pz3
  gradE(index + 8) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD3 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: px4
  gradE(index + 9) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD4 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: py4
  gradE(index + 10) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD4 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 2.0) / L;
    }, 0, 1);

  // Bending Energy: pz4
  gradE(index + 11) += coefficient * integrate<float>([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    return (bD4 * (sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) / L - 1.0) * 1.0 / sqrt(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0)) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 2.0) / L;
    }, 0, 1);
}

}  // namespace simulator
