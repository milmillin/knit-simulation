#pragma once

#include <Eigen/Core>

void circleSweep(const Eigen::MatrixXf path, const float radius,
    Eigen::MatrixXf &vertex, Eigen::MatrixXi &triangles,
    const int stride);