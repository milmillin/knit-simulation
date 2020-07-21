#pragma once

#include <functional>
#include <Eigen/Core>
#include <fstream>
#include <vector>
#include <iostream>

#include <glm/glm.hpp>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-5;

// Performs Simpson's Integration of function f over [a, b]
float integrate(const std::function<float(float)>& f, float a, float b);

// Catmull-Rom coefficient
inline float b1(float s) {
  return -s / 2 + s * s - s * s * s / 2;
}
inline float b2(float s) {
  return 1 - s * s * 5 / 2 + s * s * s * 3 / 2;
}
inline float b3(float s) {
  return s / 2 + 2 * s * s - s * s * s * 3 / 2;
}
inline float b4(float s) {
  return -s * s / 2 + s * s * s / 2;
}

// Performs Simpson's integration of function f over [a, b]
// f will be called subdivide + 1 times
// 
// subdivide : number of subdivision
template<typename T>
T integrate(const std::function<T(float)>& f, float a, float b, int subdivide = 64)
{
  assert(subdivide % 2 == 0);
  float step = (b - a) / subdivide;

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
glm::vec3 catmullRomSample(const Eigen::MatrixXf &controlPoints, int index, float s);

// Sample a Catmul-Rom curve
//
// points: one row for each point coordinate
// samplePerSegment: number of samples for each segment
// Return: samples (one row for each point coordinate)
Eigen::MatrixXf catmullRomSequenceSample(Eigen::MatrixXf points, int samplePerSegment);

void catmullRomBoundingBox(const Eigen::MatrixXf &controlPoints, int index,
  std::vector<double> *lowerBound, std::vector<double> *upperBound, float radius);

// Writes a matrix to a file in CSV format
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
void writeMatrix(std::string filename, const Eigen::MatrixXf& q);

// Reshapes #m x 3 matrix into a vector of (3 * #m) rows.
// Returns a new matrix.
Eigen::MatrixXf flatten(const Eigen::MatrixXf& v);

// Reshapes a vector of (3 * #m) rows to a #m x 3 matrix.
// Returns a new matrix.
Eigen::MatrixXf inflate(const Eigen::MatrixXf& v, size_t col = 3);

std::ostream& log();

} // namespace simulator
