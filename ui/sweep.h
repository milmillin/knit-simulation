// Copyright 2020 Tom Lou

#pragma once

#include <Eigen/Core>

namespace UI {

// Construct a tube's mesh by sweeping a circle along a path
// path: a list of 3D points defining a curve
// radius: radius of tube
// vertex: vertices for the tube
// triangles: triangles for the tube
void circleSweep(const Eigen::MatrixXd path, const double radius, const int stride,
    Eigen::MatrixXd *vertex, Eigen::MatrixXi *triangles);

}  // namespace UI
