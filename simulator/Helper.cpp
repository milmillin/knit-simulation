#include "Helper.h"

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "macros.h"

//static float simpson(const std::function<float(float)>& f, float a, float b) {
  //return (f(a) + 4 * f((a + b) / 2) + f(b))* (b - a) / 4;
//}

static inline float simpson(float fa, float fm, float fb, float a, float b) {
  return (fa + 4 * fm + fb) * (b - a) / 4;
}

static float integrateImpl(const std::function<float(float)>& f, float a, float b, float fa, float fm, float fb, int dep) {
  float m = (a + b) / 2;
  float am = (a + m) / 2;
  float mb = (m + b) / 2;

  float fam = f(am);
  float fmb = f(mb);

  float l = simpson(fa, fam, fm, a, m);
  float r = simpson(fm, fmb, fb, m, b);
  float tot = simpson(fa, fm, fb, a, b);

  if (fabs(l + r - tot) < simulator::SIMPSON_EPS || dep == 0) {
    return l + r;
  }
  return integrateImpl(f, a, m, fa, fam, fm, dep - 1) + integrateImpl(f, m, b, fm, fmb, fb, dep - 1);
}


// Helper function for `catmullRowSample`
static inline glm::vec3 bezierTerm(float t0, float t1, float t, glm::vec3& p0, glm::vec3& p1) {
  return (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1;
}

// Generate samples on a catmull-row curve
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
  glm::vec3 p0 = POINT_FROM_ROW(controlPoints, controlStartRow);
  glm::vec3 p1 = POINT_FROM_ROW(controlPoints, controlStartRow + 1);
  glm::vec3 p2 = POINT_FROM_ROW(controlPoints, controlStartRow + 2);
  glm::vec3 p3 = POINT_FROM_ROW(controlPoints, controlStartRow + 3);

  constexpr float t0 = 0;
  float t1 = t0 + glm::sqrt(glm::length(p1 - p0));
  float t2 = t1 + glm::sqrt(glm::length(p2 - p1));
  float t3 = t2 + glm::sqrt(glm::length(p3 - p2));

  for (int i = 0; i < nSamples; i++) {
    float t = t1 + (float)i / nSamples * (t2 - t1);

    auto a1 = bezierTerm(t0, t1, t, p0, p1);
    auto a2 = bezierTerm(t1, t2, t, p1, p2);
    auto a3 = bezierTerm(t2, t3, t, p2, p3);

    auto b1 = bezierTerm(t0, t2, t, a1, a2);
    auto b2 = bezierTerm(t1, t3, t, a2, a3);

    auto c = bezierTerm(t1, t2, t, b1, b2);
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

void simulator::writeMatrix(std::string filename, const Eigen::MatrixXf& q) {
  std::ofstream f;
	f.open(filename);
	f << q.format(CSVFormat) << "\n";
	f.close();
}

Eigen::MatrixXf simulator::flatten(const Eigen::MatrixXf& v) {
	Eigen::MatrixXf res(v.rows() * v.cols(), 1);
	int r = v.rows();
	int c = v.cols();
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			res(i * c + j, 0) = v(i, j);
		}
	}
	return res;
}

Eigen::MatrixXf simulator::inflate(const Eigen::MatrixXf& v, size_t col) {
	assert(v.cols() == 1);
	assert(v.rows() % col == 0);

	int r = v.rows() / col;

	Eigen::MatrixXf res(r, col);
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < col; j++) {
			res(i, j) = v(i * col + j, 0);
		}
	}
	return res;
}
