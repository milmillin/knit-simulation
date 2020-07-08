#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Global Damping
//

void Simulator::calculateGlobalDampingGradient(int i) {
  int index = i * 3;
  DECLARE_POINTS_D(index);

  float coefficient = params.kGlobal;

  // Global Damping: pxD1
  gradD(index + 0) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b1 * (b1 * pxD1 + b2 * pxD2 + b3 * pxD3 + b4 * pxD4) * 2.0;
    }, 0, 1);

  // Global Damping: pyD1
  gradD(index + 1) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b1 * (b1 * pyD1 + b2 * pyD2 + b3 * pyD3 + b4 * pyD4) * 2.0;
    }, 0, 1);

  // Global Damping: pzD1
  gradD(index + 2) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b1 * (b1 * pzD1 + b2 * pzD2 + b3 * pzD3 + b4 * pzD4) * 2.0;
    }, 0, 1);

  // Global Damping: pxD2
  gradD(index + 3) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b2 * (b1 * pxD1 + b2 * pxD2 + b3 * pxD3 + b4 * pxD4) * 2.0;
    }, 0, 1);

  // Global Damping: pyD2
  gradD(index + 4) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b2 * (b1 * pyD1 + b2 * pyD2 + b3 * pyD3 + b4 * pyD4) * 2.0;
    }, 0, 1);

  // Global Damping: pzD2
  gradD(index + 5) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b2 * (b1 * pzD1 + b2 * pzD2 + b3 * pzD3 + b4 * pzD4) * 2.0;
    }, 0, 1);

  // Global Damping: pxD3
  gradD(index + 6) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b3 * (b1 * pxD1 + b2 * pxD2 + b3 * pxD3 + b4 * pxD4) * 2.0;
    }, 0, 1);

  // Global Damping: pyD3
  gradD(index + 7) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b3 * (b1 * pyD1 + b2 * pyD2 + b3 * pyD3 + b4 * pyD4) * 2.0;
    }, 0, 1);

  // Global Damping: pzD3
  gradD(index + 8) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b3 * (b1 * pzD1 + b2 * pzD2 + b3 * pzD3 + b4 * pzD4) * 2.0;
    }, 0, 1);

  // Global Damping: pxD4
  gradD(index + 9) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b4 * (b1 * pxD1 + b2 * pxD2 + b3 * pxD3 + b4 * pxD4) * 2.0;
    }, 0, 1);

  // Global Damping: pyD4
  gradD(index + 10) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b4 * (b1 * pyD1 + b2 * pyD2 + b3 * pyD3 + b4 * pyD4) * 2.0;
    }, 0, 1);

  // Global Damping: pzD4
  gradD(index + 11) += coefficient * integrate([&](float s)->float {
    DECLARE_BASIS
    return b4 * (b1 * pzD1 + b2 * pzD2 + b3 * pzD3 + b4 * pzD4) * 2.0;
    }, 0, 1);
}


} // namespace simulator
