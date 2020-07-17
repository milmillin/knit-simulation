#pragma once

#include <functional>
#include <Eigen/Core>
#include <fstream>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-5;

// Performs adaptive Simpson's integration of function f over [a, b]
// 
// subdivide : minimum subdivision
float integrate(const std::function<float(float)>& f, float a, float b, int subdivide = 1, int maxDepth = 5);

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
