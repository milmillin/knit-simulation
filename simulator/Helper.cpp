#include "Helper.h"

static float simpson(const std::function<float(float)>& f, float a, float b) {
  return (f(a) + 4 * f((a + b) / 2) + f(b))* (b - a) / 4;
}

static float integrateImpl(const std::function<float(float)>& f, float a, float b, float tot) {
  float m = (a + b) / 2;
  float l = simpson(f, a, m);
  float r = simpson(f, m, b);
  if (fabs(l + r - tot) < simulator::SIMPSON_EPS) {
    return l + r;
  }
  return integrateImpl(f, a, m, l) + integrateImpl(f, m, b, r);
}

float simulator::integrate(const std::function<float(float)>& f, float a, float b)
{
  return integrateImpl(f, a, b, simpson(f, a, b));
}

float simulator::catmullRomArcLength(const Eigen::MatrixXf& q, int index)
{
  const float& px1 = q(index);
  const float& py1 = q(index + 1);
  const float& pz1 = q(index + 2);
  const float& px2 = q(index + 3);
  const float& py2 = q(index + 4);
  const float& pz2 = q(index + 5);
  const float& px3 = q(index + 6);
  const float& py3 = q(index + 7);
  const float& pz3 = q(index + 8);
  const float& px4 = q(index + 9);
  const float& py4 = q(index + 10);
  const float& pz4 = q(index + 11);

  float k0 = (px1 * px1) / 4.0 + (px3 * px3) / 4.0 + (py1 * py1) / 4.0 + (py3 * py3) / 4.0 + (pz1 * pz1) / 4.0 + (pz3 * pz3) / 4.0;
  float k1 = -(px1 * px2 * -5.0 + px1 * px3 * 2.0 - px1 * px4 + px2 * px3 * 5.0 + px3 * px4 - py1 * py2 * 5.0 + py1 * py3 * 2.0 - py1 * py4 + py2 * py3 * 5.0 + py3 * py4 - pz1 * pz2 * 5.0 + pz1 * pz3 * 2.0 - pz1 * pz4 + pz2 * pz3 * 5.0 + pz3 * pz4 + (px1 * px1) * 2.0 - (px3 * px3) * 4.0 + (py1 * py1) * 2.0 - (py3 * py3) * 4.0 + (pz1 * pz1) * 2.0 - (pz3 * pz3) * 4.0);
  float k2 = (px1 * px2 * (-4.9E+1 / 2.0) + px1 * px3 * 1.9E+1 - px1 * px4 * (1.1E+1 / 2.0) - px2 * px3 * (7.1E+1 / 2.0) + px2 * px4 * 1.0E+1 - px3 * px4 * (1.3E+1 / 2.0) - py1 * py2 * (4.9E+1 / 2.0) + py1 * py3 * 1.9E+1 - py1 * py4 * (1.1E+1 / 2.0) - py2 * py3 * (7.1E+1 / 2.0) + py2 * py4 * 1.0E+1 - py3 * py4 * (1.3E+1 / 2.0) - pz1 * pz2 * (4.9E+1 / 2.0) + pz1 * pz3 * 1.9E+1 - pz1 * pz4 * (1.1E+1 / 2.0) - pz2 * pz3 * (7.1E+1 / 2.0) + pz2 * pz4 * 1.0E+1 - pz3 * pz4 * (1.3E+1 / 2.0) + (px1 * px1) * (1.1E+1 / 2.0) + (px2 * px2) * 2.5E+1 + (px3 * px3) * (2.3E+1 / 2.0) + px4 * px4 + (py1 * py1) * (1.1E+1 / 2.0) + (py2 * py2) * 2.5E+1 + (py3 * py3) * (2.3E+1 / 2.0) + py4 * py4 + (pz1 * pz1) * (1.1E+1 / 2.0) + (pz2 * pz2) * 2.5E+1 + (pz3 * pz3) * (2.3E+1 / 2.0) + pz4 * pz4) - (px1 * px3) / 2.0 - (py1 * py3) / 2.0 - (pz1 * pz3) / 2.0;
  float k3 = -(px1 * px2 * -3.3E+1 + px1 * px3 * 3.0E+1 - px1 * px4 * 9.0 - px2 * px3 * 8.1E+1 + px2 * px4 * 2.4E+1 - px3 * px4 * 2.1E+1 - py1 * py2 * 3.3E+1 + py1 * py3 * 3.0E+1 - py1 * py4 * 9.0 - py2 * py3 * 8.1E+1 + py2 * py4 * 2.4E+1 - py3 * py4 * 2.1E+1 - pz1 * pz2 * 3.3E+1 + pz1 * pz3 * 3.0E+1 - pz1 * pz4 * 9.0 - pz2 * pz3 * 8.1E+1 + pz2 * pz4 * 2.4E+1 - pz3 * pz4 * 2.1E+1 + (px1 * px1) * 6.0 + (px2 * px2) * 4.5E+1 + (px3 * px3) * 3.6E+1 + (px4 * px4) * 3.0 + (py1 * py1) * 6.0 + (py2 * py2) * 4.5E+1 + (py3 * py3) * 3.6E+1 + (py4 * py4) * 3.0 + (pz1 * pz1) * 6.0 + (pz2 * pz2) * 4.5E+1 + (pz3 * pz3) * 3.6E+1 + (pz4 * pz4) * 3.0);
  float k4 = (px1 * px2 * (-2.7E+1 / 2.0) + px1 * px3 * (2.7E+1 / 2.0) - px1 * px4 * (9.0 / 2.0) - px2 * px3 * (8.1E+1 / 2.0) + px2 * px4 * (2.7E+1 / 2.0) - px3 * px4 * (2.7E+1 / 2.0) - py1 * py2 * (2.7E+1 / 2.0) + py1 * py3 * (2.7E+1 / 2.0) - py1 * py4 * (9.0 / 2.0) - py2 * py3 * (8.1E+1 / 2.0) + py2 * py4 * (2.7E+1 / 2.0) - py3 * py4 * (2.7E+1 / 2.0) - pz1 * pz2 * (2.7E+1 / 2.0) + pz1 * pz3 * (2.7E+1 / 2.0) - pz1 * pz4 * (9.0 / 2.0) - pz2 * pz3 * (8.1E+1 / 2.0) + pz2 * pz4 * (2.7E+1 / 2.0) - pz3 * pz4 * (2.7E+1 / 2.0) + (px1 * px1) * (9.0 / 4.0) + (px2 * px2) * (8.1E+1 / 4.0) + (px3 * px3) * (8.1E+1 / 4.0) + (px4 * px4) * (9.0 / 4.0) + (py1 * py1) * (9.0 / 4.0) + (py2 * py2) * (8.1E+1 / 4.0) + (py3 * py3) * (8.1E+1 / 4.0) + (py4 * py4) * (9.0 / 4.0) + (pz1 * pz1) * (9.0 / 4.0) + (pz2 * pz2) * (8.1E+1 / 4.0) + (pz3 * pz3) * (8.1E+1 / 4.0) + (pz4 * pz4) * (9.0 / 4.0));

  return integrate([&](float s)->float {
    return k0 + k1 * s + k2 * (s * s) + k3 * (s * s * s) + k4 * (s * s * s * s);
    }, 0, 1);
}
