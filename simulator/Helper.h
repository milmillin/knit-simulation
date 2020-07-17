#pragma once

#include <functional>
#include <Eigen/Core>
#include <fstream>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-5;

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

// Sample a Catmul-Rom curve
//
// points: one row for each point coordinate
// samplePerSegment: number of samples for each segment
// Return: samples (one row for each point coordinate)
Eigen::MatrixXf catmullRomSequenceSample(Eigen::MatrixXf points, int samplePerSegment);

// Writes a matrix to a file in CSV format
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
void writeMatrix(std::string filename, const Eigen::MatrixXf& q);

// Reshapes #m x 3 matrix into a vector of (3 * #m) rows.
// Returns a new matrix.
Eigen::MatrixXf flatten(const Eigen::MatrixXf& v);

// Reshapes a vector of (3 * #m) rows to a #m x 3 matrix.
// Returns a new matrix.
Eigen::MatrixXf inflate(const Eigen::MatrixXf& v, size_t col = 3);

}; // namespace simulator
