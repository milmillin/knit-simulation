#include "./Helper.h"

#include <vector>
#include <ctime>
#include <functional>

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "./macros.h"

//static double simpson(const std::function<double(double)>& f, double a, double b) {
  //return (f(a) + 4 * f((a + b) / 2) + f(b))* (b - a) / 4;
//}

static inline double simpson(double fa, double fm, double fb, double a, double b) {
  return (fa + 4 * fm + fb) * (b - a) / 4;
}

static double integrateImpl(const std::function<double(double)>& f, double a, double b, double fa, double fm, double fb, int dep) {
  double m = (a + b) / 2;
  double am = (a + m) / 2;
  double mb = (m + b) / 2;

  double fam = f(am);
  double fmb = f(mb);

  double l = simpson(fa, fam, fm, a, m);
  double r = simpson(fm, fmb, fb, m, b);
  double tot = simpson(fa, fm, fb, a, b);

  if (fabs(l + r - tot) < simulator::SIMPSON_EPS || dep == 0) {
    return l + r;
  }
  return integrateImpl(f, a, m, fa, fam, fm, dep - 1) + integrateImpl(f, m, b, fm, fmb, fb, dep - 1);
}



glm::dvec3 simulator::catmullRomSample(const Eigen::MatrixXd& controlPoints, int index, double s) {
  glm::dvec3 c0 = POINT_FROM_ROW(controlPoints, index);
  glm::dvec3 c1 = POINT_FROM_ROW(controlPoints, index + 1);
  glm::dvec3 c2 = POINT_FROM_ROW(controlPoints, index + 2);
  glm::dvec3 c3 = POINT_FROM_ROW(controlPoints, index + 3);

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
(const Eigen::MatrixXd& controlPoints, int controlStartRow, int nSamples,
  Eigen::MatrixXd* samples, int samplesStartRow) {
  for (int i = 0; i < nSamples; i++) {
    double s = (double)i / nSamples;
    auto c = simulator::catmullRomSample(controlPoints, controlStartRow, s);
    ROW_FROM_POINT(*samples, samplesStartRow, c);
    samplesStartRow++;
  }
}

Eigen::MatrixXd simulator::catmullRomSequenceSample(Eigen::MatrixXd points, int samplePerSegment) {
  int nPoints = points.rows();
  Eigen::MatrixXd result((nPoints - 3) * samplePerSegment + 1, 3);

  for (int i = 0; i < nPoints - 3; i++) {
    catmullRowSample(points, i, samplePerSegment, &result, i * samplePerSegment);
  }
  glm::dvec3 lastPoint = POINT_FROM_ROW(points, nPoints - 2);
  ROW_FROM_POINT(result, result.rows() - 1, lastPoint);

  return result;
}

// See https://pomax.github.io/bezierinfo/#catmullconv
void simulator::catmullRomBoundingBox(const Eigen::MatrixXd& points, int index,
  std::vector<double>& lowerBound, std::vector<double>& upperBound, double radius) {
  Eigen::MatrixXd controlPoints(3, 4);
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

void simulator::writeMatrix(std::string filename, const Eigen::MatrixXd& q) {
  std::ofstream f;
  f.open(filename);
  f << q.format(CSVFormat) << "\n";
  f.close();
}

Eigen::MatrixXd simulator::flatten(const Eigen::MatrixXd& v) {
  Eigen::MatrixXd res(v.rows() * v.cols(), 1);
  int r = v.rows();
  int c = v.cols();
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      res(i * c + j, 0) = v(i, j);
    }
  }
  return res;
}

Eigen::MatrixXd simulator::inflate(const Eigen::MatrixXd& v, size_t col) {
  assert(v.cols() == 1);
  assert(v.rows() % col == 0);

  int r = v.rows() / col;

  Eigen::MatrixXd res(r, col);
  for (int i = 0; i < r; i++) {
    for (size_t j = 0; j < col; j++) {
      res(i, j) = v(i * col + j, 0);
    }
  }
  return res;
}

double& simulator::coordAt(Eigen::MatrixXd& q, int index, int axis) {
  return q(index * 3ll + axis);
}

const double& simulator::coordAt(const Eigen::MatrixXd& q, int index, int axis) {
  return q(index * 3ll + axis);
}

double simulator::maxCoeff(const Eigen::MatrixXd& m) {
  if (m.rows() == 0 || m.cols() == 0) return 0;
  return std::max(fabs(m.maxCoeff()), fabs(m.minCoeff()));
}

std::string simulator::toString(Eigen::MatrixXd x) {
  std::string s = "[";
  for (int r = 0; r < x.rows(); r++) {
    s += "[";
    for (int c = 0; c < x.cols(); c++) {
      s += std::to_string(x(r, c)) + ", ";
    }
    s += "]";
    if (r != x.rows() - 1) {
      s += ", ";
    }
  }
  s += "]";
  return s;
}

Eigen::Vector3d simulator::parallelTransport(const Eigen::Vector3d& u, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2) {
  Eigen::Vector3d t1 = e1 / e1.norm();
  Eigen::Vector3d t2 = e2 / e2.norm();
  Eigen::Vector3d n = t1.cross(t2);
  if (n.norm() < 1e-10)
    return u;
  n /= n.norm();
  Eigen::Vector3d p1 = n.cross(t1);
  Eigen::Vector3d p2 = n.cross(t2);
  return u.dot(n) * n + u.dot(t1) * t2 + u.dot(p1) * p2;
}

static Eigen::RowVector3d getPointFromTime
  (const Eigen::Block<const Eigen::MatrixXd, 4, 3>& controlPoints, double s) {
  Eigen::RowVector3d vec = Eigen::RowVector3d::Zero();
  DECLARE_BASIS2(b, s);
  for (size_t i = 0; i < 4; i++) {
    vec += controlPoints.row(i) * b[i];
  }
  return vec;
}

static double calculateSegmentLength(const Eigen::Block<const Eigen::MatrixXd, 4, 3>& controlPoints,
  double a, double b) {
  Eigen::RowVector3d vec;
  return simulator::integrate<double>([&](double s)->double {
    DECLARE_BASIS_D2(bD, s);
    vec.setZero();
    for (size_t i = 0; i < 4; i++) {
      vec += controlPoints.row(i) * bD[i];
    }
    return vec.norm();
    }, a, b);
}

static const double TIME_EPS = 1e-5;
static double findTimeAtLength(const Eigen::Block<const Eigen::MatrixXd, 4, 3>& controlPoints, double length) {
  double lo = 0;
  double hi = 1;
  double mid;
  double toLoLength = 0;
  double curLength;
  while (hi - lo > TIME_EPS) {
    mid = (lo + hi) / 2;
    curLength = toLoLength + calculateSegmentLength(controlPoints, lo, mid);
    if (curLength <= length) {
      lo = mid;
      toLoLength = curLength;
    }
    else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

file_format::YarnRepr simulator::reparameterizeYarns(const file_format::YarnRepr& yarns, double avgLengthFactor, double* L)
{
  std::vector<double> yarnLength(yarns.yarns.size(), 0);
  std::vector<std::vector<double>> segmentLength(yarns.yarns.size());
  double totalLength = 0;
  size_t segmentCount = 0;

  for (size_t y = 0; y < yarns.yarns.size(); y++) {
    const auto& yarn = yarns.yarns[y];
    segmentLength[y].resize(std::max(yarn.size() - 3, 0ull));

    for (size_t i = yarn.begin; i < yarn.end - 3; i++) {
      double currentLength = calculateSegmentLength(yarns.vertices.block<4, 3>(i, 0), 0, 1);
      segmentLength[y][i - yarn.begin] = currentLength;
      yarnLength[y] += currentLength;
      totalLength += currentLength;
      segmentCount++;
    }
  }

  double avgLength = totalLength / segmentCount;
  *L = avgLengthFactor * avgLength;

  std::vector<Eigen::RowVector3d> newVertices;
  file_format::YarnRepr newYarns;
  newYarns.yarns = yarns.yarns;

  for (size_t y = 0; y < yarns.yarns.size(); y++) {
    const auto& yarn = yarns.yarns[y];
    auto& newYarn = newYarns.yarns[y];

    double offset = std::fmod(yarnLength[y], *L) / 2;
    size_t nNewPoints = (size_t)std::floor(yarnLength[y] / *L) + 1;
    // FIXME
    if (yarn.size() < 3 || nNewPoints == 0) {
      newYarn.begin = newYarn.end = newVertices.size();
      continue;
    }
    size_t nSegment = yarn.size() - 3;
    size_t i = 0;
    double curT = 0;
    double curLength = 0;
    auto advanceLength = [&](double length) {
      double toAdvance = length;

      // advance by segment
      while (i < nSegment && (segmentLength[y][i] - curLength) <= toAdvance) {
        toAdvance -= segmentLength[y][i] - curLength;
        curLength = 0;
        curT = 0;
        i++;
      }
      assert(i < nSegment || toAdvance <= 1e-5);

      // advance within segment
      if (toAdvance > 1e-5) {
        curT = findTimeAtLength(yarns.vertices.block<4, 3>(i + yarn.begin, 0), curLength + toAdvance);
        curLength += toAdvance;
      }
    };

    advanceLength(offset);
    newYarn.begin = newVertices.size();
    newYarn.end = newYarn.begin + nNewPoints;
    newVertices.push_back(getPointFromTime(yarns.vertices.block<4, 3>(i + yarn.begin, 0), curT));
    for (size_t ii = 1; ii < nNewPoints; ii++) {
      advanceLength(*L);
      newVertices.push_back(getPointFromTime(yarns.vertices.block<4, 3>(i + yarn.begin, 0), curT));
    }
  }

  newYarns.vertices.resize(newVertices.size(), 3);
  for (size_t i = 0; i < newVertices.size(); i++) {
    newYarns.vertices.row(i) = newVertices[i];
  }
  return newYarns;
}
