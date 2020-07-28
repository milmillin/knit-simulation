#include "Helper.h"

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include <vector>
#include <ctime>

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
void simulator::catmullRomBoundingBox(const Eigen::MatrixXf& points, int index,
    std::vector<double>& lowerBound, std::vector<double>& upperBound, float radius) {
  Eigen::MatrixXf controlPoints(3, 4);
  controlPoints.col(0) = pointAt(points, index + 1);
  controlPoints.col(1) = pointAt(points, index + 1) + (pointAt(points, index + 2) - pointAt(points, index + 0)) / 3;
  controlPoints.col(2) = pointAt(points, index + 2) - (pointAt(points, index + 3) - pointAt(points, index + 1)) / 3;
  controlPoints.col(3) = pointAt(points, index + 2);
  lowerBound.resize(3);
  upperBound.resize(3);
  for (int i = 0; i < 3; i++) {
    lowerBound[i] = controlPoints.row(i).minCoeff() - radius;
    upperBound[i] = controlPoints.row(i).maxCoeff() + radius;
  }
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

std::ostream& simulator::log()
{
	char buf[20];
	time_t rawTime;
	struct tm* timeInfo;

	time(&rawTime);
  timeInfo = localtime(&rawTime);

	strftime(buf, 20, "%D %T", timeInfo);
	return std::cout << "[" << buf << "] ";
}

Eigen::Block<Eigen::MatrixXf, 3, 1> simulator::pointAt(Eigen::MatrixXf& q, int index) {
  return q.block<3, 1>(index * 3ll, 0);
}

Eigen::Block<const Eigen::MatrixXf, 3, 1> simulator::pointAt(const Eigen::MatrixXf& q, int index) {
  return q.block<3, 1>(index * 3ll, 0);
}

float& simulator::coordAt(Eigen::MatrixXf& q, int index, int axis) {
  return q(index * 3ll + axis);
}

const float& simulator::coordAt(const Eigen::MatrixXf& q, int index, int axis) {
  return q(index * 3ll + axis);
}

float simulator::maxCoeff(const Eigen::MatrixXf& m) {
	if (m.rows() == 0 || m.cols() == 0) return 0;
	return std::max(fabs(m.maxCoeff()), fabs(m.minCoeff()));
}
