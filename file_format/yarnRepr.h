#pragma once

#include <vector>

#include <Eigen/Core>
#include <glm/common.hpp>

#include "./yarns.h"

namespace file_format {

// An internal representation of a yarn in the simulator
struct Yarn {
  // A list of 3D control points (row based)
  Eigen::Block<Eigen::MatrixXd, -1, -1, false> points;
  // RGB color (in range 0-1)
  Eigen::RowVector3d color = Eigen::RowVector3d(0.5, 0.5, 0.5);
  // Yarn radius
  double radius;
};

// An internal representation of yarns in the simulator 
class YarnRepr {
public:
  struct YarnInfo {
    size_t begin = 0;
    size_t end = 0; // exclusive
    Eigen::RowVector3d color = Eigen::RowVector3d(0.5, 0.5, 0.5);
    double radius = 0.1;

    inline size_t size() const { return end - begin; }
  };

  // Empty constructor
  YarnRepr() {}
  // Convert `.yarns` file to internal representation
  YarnRepr(file_format::Yarns::Yarns& yarns);
  // Clone YarnRepr without cloning points
  YarnRepr createAlike() const;

  Yarns::Yarns toYarns() const;

  // Vertices
  Eigen::MatrixXd vertices;
  
  Eigen::Block<Eigen::MatrixXd> getYarnPoints(size_t index);
  const Eigen::Block<const Eigen::MatrixXd> getYarnPoints(size_t index) const;

  size_t numYarns() const { return yarns.size(); }

  std::vector<YarnInfo> yarns;
};

}  // namespace file_format
