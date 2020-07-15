#include "Simulator.h"
#include "Helper.h"

namespace simulator {

//////////////////////////////////////////////
//
// Collision Energy
//

static constexpr int collisionSubdivide = 4;
static constexpr int collisionMaxDepth = 0;

void Simulator::calculateCollisionEnergyGradient(int i, int j) {
  int iIndex = i * 3;
  int jIndex = j * 3;

  float r = yarns.yarns[0].radius;
  float coefficient = params.kContact * segmentLength[i] * segmentLength[j];
  float constR = 4.0f * r * r;

  // calculate only for j > i
  coefficient *= 2;

  DECLARE_POINTS(pi, iIndex);
  DECLARE_POINTS(pj, jIndex);

  // Collision Energy: pix1
  gradE(iIndex + 0) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi1 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piy1
  gradE(iIndex + 1) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi1 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piz1
  gradE(iIndex + 2) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi1 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pix2
  gradE(iIndex + 3) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi2 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piy2
  gradE(iIndex + 4) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi2 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piz2
  gradE(iIndex + 5) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi2 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pix3
  gradE(iIndex + 6) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi3 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piy3
  gradE(iIndex + 7) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi3 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piz3
  gradE(iIndex + 8) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi3 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pix4
  gradE(iIndex + 9) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi4 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piy4
  gradE(iIndex + 10) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi4 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: piz4
  gradE(iIndex + 11) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bi4 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * 2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjx1
  gradE(jIndex + 0) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj1 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjy1
  gradE(jIndex + 1) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj1 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjz1
  gradE(jIndex + 2) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj1 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjx2
  gradE(jIndex + 3) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj2 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjy2
  gradE(jIndex + 4) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj2 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjz2
  gradE(jIndex + 5) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj2 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjx3
  gradE(jIndex + 6) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj3 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjy3
  gradE(jIndex + 7) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj3 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjz3
  gradE(jIndex + 8) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj3 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjx4
  gradE(jIndex + 9) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj4 * (bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjy4
  gradE(jIndex + 10) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj4 * (bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

  // Collision Energy: pjz4
  gradE(jIndex + 11) += coefficient * integrate([=](float si)->float {
    DECLARE_BASIS(bi, si);
    return integrate([=](float sj)->float {
      DECLARE_BASIS(bj, sj);
      float norm = pow(bi1 * pix1 + bi2 * pix2 + bi3 * pix3 + bi4 * pix4 - bj1 * pjx1 - bj2 * pjx2 - bj3 * pjx3 - bj4 * pjx4, 2.0) + pow(bi1 * piy1 + bi2 * piy2 + bi3 * piy3 + bi4 * piy4 - bj1 * pjy1 - bj2 * pjy2 - bj3 * pjy3 - bj4 * pjy4, 2.0) + pow(bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4, 2.0);
      if (norm >= constR) return 0;
      float normD = bj4 * (bi1 * piz1 + bi2 * piz2 + bi3 * piz3 + bi4 * piz4 - bj1 * pjz1 - bj2 * pjz2 - bj3 * pjz3 - bj4 * pjz4) * -2.0;
      return -constR * normD / (norm * norm) + normD / constR;
      }, 0, 1, collisionSubdivide, collisionMaxDepth);
    }, 0, 1, collisionSubdivide, collisionMaxDepth);

}

}  // namespace simulator
