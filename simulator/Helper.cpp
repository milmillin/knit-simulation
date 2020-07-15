#include "Helper.h"

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include <vector>

#include "macros.h"

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


glm::vec3 simulator::catmullRomSample(const Eigen::MatrixXf &controlPoints, int index, float s) {
  glm::vec3 c0 = POINT_FROM_ROW(controlPoints, index);
  glm::vec3 c1 = POINT_FROM_ROW(controlPoints, index + 1);
  glm::vec3 c2 = POINT_FROM_ROW(controlPoints, index + 2);
  glm::vec3 c3 = POINT_FROM_ROW(controlPoints, index + 3);

  return simulator::b1(s) * c0 + simulator::b2(s) * c1 + simulator::b3(s) * c2 + simulator::b4(s) * c3;
}

// Generate samples on a catmull-rom curve
// Helper function for `catmullRomSequenceSample`
//
// controlPoints: one row per control points
// controlStartRow: the index range [controlStartRow, controlStartRow + 4) will be the control points
// nSamples: number of samples to generate
// samples: output variable to store the result
// samplesStartRow: store the samples in the index range [samplesStartRow, samplesStartRow + nSamples).
static inline void catmullRowSample
    (const Eigen::MatrixXf &controlPoints, int controlStartRow, int nSamples,
    Eigen::MatrixXf *samples, int samplesStartRow) {
  glm::vec3 c0 = POINT_FROM_ROW(controlPoints, controlStartRow);
  glm::vec3 c1 = POINT_FROM_ROW(controlPoints, controlStartRow + 1);
  glm::vec3 c2 = POINT_FROM_ROW(controlPoints, controlStartRow + 2);
  glm::vec3 c3 = POINT_FROM_ROW(controlPoints, controlStartRow + 3);

  for (int i = 0; i < nSamples; i++) {
    float s = (float) i / nSamples;
    auto c = simulator::catmullRomSample(controlPoints, controlStartRow, s);
    ROW_FROM_POINT(*samples, samplesStartRow, c);
    samplesStartRow++;
  }
}

Eigen::MatrixXf simulator::catmullRomSequenceSample(Eigen::MatrixXf points, int samplePerSegment) {
  int nPoints = points.rows();
  Eigen::MatrixXf result((nPoints - 3) * samplePerSegment + 1, 3);

  for (int i = 0; i < nPoints - 3; i++) {
    catmullRowSample(points, i, samplePerSegment, &result, i * samplePerSegment);
  }
  glm::vec3 lastPoint = POINT_FROM_ROW(points, nPoints - 2);
  ROW_FROM_POINT(result, result.rows() - 1, lastPoint);

  return result;
}

// See https://pomax.github.io/bezierinfo/#catmullconv
void simulator::catmullRomBoundingBox(const Eigen::MatrixXf &points, int index,
    std::vector<double> *lowerBound, std::vector<double> *upperBound, float radius) {
  Eigen::MatrixXf bezierControlPoints(4, 3);
  bezierControlPoints.row(0) = points.row(index + 1);
  bezierControlPoints.row(1) = points.row(index + 1) + (points.row(index + 2) - points.row(index + 0)) / 3;
  bezierControlPoints.row(2) = points.row(index + 2) - (points.row(index + 3) - points.row(index + 1)) / 3;
  bezierControlPoints.row(3) = points.row(index + 2);
  lowerBound->resize(3);
  upperBound->resize(3);
  for (int i = 0; i < 3; i++) {
    (*lowerBound)[i] = bezierControlPoints.col(i).minCoeff() - radius;
    (*upperBound)[i] = bezierControlPoints.col(i).maxCoeff() + radius;
  }
}
