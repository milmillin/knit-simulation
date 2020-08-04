#pragma once

#include <functional>
#include <fstream>
#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <glm/glm.hpp>

namespace simulator {

constexpr double SIMPSON_EPS = 1e-5;

// Catmull-Rom coefficient
inline double b1(double s) {
  return -s / 2 + s * s - s * s * s / 2;
}
inline double b2(double s) {
  return 1 - s * s * 5 / 2 + s * s * s * 3 / 2;
}
inline double b3(double s) {
  return s / 2 + 2 * s * s - s * s * s * 3 / 2;
}
inline double b4(double s) {
  return -s * s / 2 + s * s * s / 2;
}

// Performs Simpson's integration of function f over [a, b]
// f will be called subdivide + 1 times
// 
// subdivide : number of subdivision
template<typename T>
T integrate(const std::function<T(double)>& f, double a, double b, int subdivide = 16)
{
  assert(subdivide % 2 == 0);
  double step = (b - a) / subdivide;

  T result = f(a) + f(b);

  for (int i = 1; i < subdivide; i += 2) {
    result += 4 * f(a + i * step);
  }
  
  for (int i = 2; i < subdivide; i += 2) {
    result += 2 * f(a + i * step);
  }

  return result * step / 3;
}

// Calculates Arc length of a Catmull-Rom Spline defined by 4 control points.
glm::dvec3 catmullRomSample(const Eigen::MatrixXd &controlPoints, int index, double s);

// Sample a Catmul-Rom curve
//
// points: one row for each point coordinate
// samplePerSegment: number of samples for each segment
// Return: samples (one row for each point coordinate)
Eigen::MatrixXd catmullRomSequenceSample(Eigen::MatrixXd points, int samplePerSegment);

void catmullRomBoundingBox(const Eigen::MatrixXd& controlPoints, int index,
  std::vector<double>& lowerBound, std::vector<double>& upperBound, double radius);

// Writes a matrix to a file in CSV format
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
void writeMatrix(std::string filename, const Eigen::MatrixXd& q);

// Reshapes #m x 3 matrix into a vector of (3 * #m) rows.
// Returns a new matrix.
Eigen::MatrixXd flatten(const Eigen::MatrixXd& v);

// Reshapes a vector of (3 * #m) rows to a #m x 3 matrix.
// Returns a new matrix.
Eigen::MatrixXd inflate(const Eigen::MatrixXd& v, size_t col = 3);

std::ostream& log();

Eigen::Block<Eigen::MatrixXd, 3, 1> pointAt(Eigen::MatrixXd& q, int index);
Eigen::Block<const Eigen::MatrixXd, 3, 1> pointAt(const Eigen::MatrixXd& q, int index);

double& coordAt(Eigen::MatrixXd& q, int index, int axis);
const double& coordAt(const Eigen::MatrixXd& q, int index, int axis);

double maxCoeff(const Eigen::MatrixXd& m);

} // namespace simulator
