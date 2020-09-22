#pragma once

#include <vector>

#include <Eigen/Core>
#include <glm/common.hpp>

#include "./yarns.h"

namespace file_format {

// An internal representation of a yarn in the simulator
struct Yarn {
  // A list of 3D control points (row based)
  Eigen::MatrixXd points;
  // RGB color (in range 0-255)
  glm::ivec3 color = glm::ivec3(128, 128, 128);
  // Bishop frame directions (unit vector or zero)
  // The frame at index `i` is the frame for the segment between
  // `points.row(i)` and `points.row(i+1)`
  Eigen::MatrixX3d bishopFrameU;
  Eigen::MatrixX3d bishopFrameV;
  // Material frame directions (unit vector or zero)
  // The frame at index `i` is the frame for the segment between
  // `points.row(i)` and `points.row(i+1)`
  Eigen::MatrixX3d materialFrameU;
  Eigen::MatrixX3d materialFrameV;
  // Yarn radius
  double radius;
};

// An internal representation of yarns in the simulator 
class YarnRepr {
 public:
  // Empty constructor
  YarnRepr() {}
  // Convert `.yarns` file to internal representation
  YarnRepr(file_format::Yarns::Yarns &yarns);
  // Clone YarnRepr without cloning points
  YarnRepr createAlike() const;

  Yarns::Yarns toYarns() const;

  // Stores a list of yarns
  std::vector<Yarn> yarns;
};

}  // namespace file_format
