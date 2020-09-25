#pragma once

#include <functional>
#include <fstream>
#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <glm/glm.hpp>

#include "../file_format/yarnRepr.h"

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

#define pointAt(q, index) (q).block<3,1>((index) * 3ull, 0)

double& coordAt(Eigen::MatrixXd& q, int index, int axis);
const double& coordAt(const Eigen::MatrixXd& q, int index, int axis);

double maxCoeff(const Eigen::MatrixXd& m);

std::string toString(Eigen::MatrixXd x);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d& u, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2);

// Reparameterize the Catmull-Rom curve by
// (1) finding points in the curve such that they are equal length L apart;
// (2) creating a new curve with these points; this will approximate the reparameterization.
// Length L is determined by "avgLengthFactor" * average length of each pair of adjacent control points.
// Since the length of the whole curve might not be divisible by L, the first and the last points will
// be at the same offset from their original endpoints (i.e., the first point will be located at 
// (curve length % L) / 2).
// Returns a new YarnRepr and the length L
file_format::YarnRepr reparameterizeYarns(const file_format::YarnRepr& yarn, double avgLengthFactor, double* L);

} // namespace simulator
