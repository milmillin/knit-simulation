#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Bending Energy
//

void Simulator::calculateBendingEnergyGradient(int i) {
  int index = i * 3;
  DECLARE_POINTS(p, index);

  float coefficient = params.kBend * segmentLength[i];

  // Bending Energy: px1
  gradE(index + 0) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return -((bDD1 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD1 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 + (bDD1 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD1 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD1 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 6.0;
    }, 0, 1);

  // Bending Energy: py1
  gradE(index + 1) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD1 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD1 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 - (bDD1 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD1 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD1 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 6.0;
    }, 0, 1);

  // Bending Energy: pz1
  gradE(index + 2) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD1 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD1 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0 + (bDD1 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD1 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD1 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 6.0;
    }, 0, 1);

  // Bending Energy: px2
  gradE(index + 3) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return -((bDD2 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD2 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 + (bDD2 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD2 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD2 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 6.0;
    }, 0, 1);

  // Bending Energy: py2
  gradE(index + 4) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD2 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD2 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 - (bDD2 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD2 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD2 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 6.0;
    }, 0, 1);

  // Bending Energy: pz2
  gradE(index + 5) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD2 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD2 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0 + (bDD2 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD2 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD2 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 6.0;
    }, 0, 1);

  // Bending Energy: px3
  gradE(index + 6) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return -((bDD3 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD3 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 + (bDD3 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD3 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD3 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 6.0;
    }, 0, 1);

  // Bending Energy: py3
  gradE(index + 7) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD3 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD3 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 - (bDD3 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD3 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD3 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 6.0;
    }, 0, 1);

  // Bending Energy: pz3
  gradE(index + 8) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD3 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD3 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0 + (bDD3 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD3 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD3 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 6.0;
    }, 0, 1);

  // Bending Energy: px4
  gradE(index + 9) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return -((bDD4 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD4 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 + (bDD4 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD4 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD4 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * 6.0;
    }, 0, 1);

  // Bending Energy: py4
  gradE(index + 10) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD4 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD4 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4)) * 2.0 - (bDD4 * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) - bD4 * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD4 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * 6.0;
    }, 0, 1);

  // Bending Energy: pz4
  gradE(index + 11) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS_D(b, s);
    DECLARE_BASIS_DD(b, s);
    return ((bDD4 * (bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) - bD4 * (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4)) * ((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0 + (bDD4 * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) - bD4 * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4)) * ((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4)) * 2.0) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 3.0) - bD4 * (pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4), 2.0) + pow((bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * px1 + bDD2 * px2 + bDD3 * px3 + bDD4 * px4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0) + pow((bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4) * (bDD1 * pz1 + bDD2 * pz2 + bDD3 * pz3 + bDD4 * pz4) - (bDD1 * py1 + bDD2 * py2 + bDD3 * py3 + bDD4 * py4) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4), 2.0)) * 1.0 / pow(pow(bD1 * px1 + bD2 * px2 + bD3 * px3 + bD4 * px4, 2.0) + pow(bD1 * py1 + bD2 * py2 + bD3 * py3 + bD4 * py4, 2.0) + pow(bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4, 2.0), 4.0) * (bD1 * pz1 + bD2 * pz2 + bD3 * pz3 + bD4 * pz4) * 6.0;
    }, 0, 1);

}

}  // namespace simulator
