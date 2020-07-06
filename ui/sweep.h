// Copyright 2020 Tom Lou

#pragma once

#include <Eigen/Core>

namespace UI {

// Construct a tube's mesh by sweeping a circle along a path
// path: a list of 3D points defining a curve
// radius: radius of tube
// vertex: vertices for the tube
// triangles: triangles for the tube
void circleSweep(const Eigen::MatrixXf path, const float radius, const int stride,
    Eigen::MatrixXf *vertex, Eigen::MatrixXi *triangles);

}  // namespace UI
