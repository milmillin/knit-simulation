#include "./yarnRepr.h"

#include <Eigen/Core>
#include <glm/common.hpp>

#include <vector>

#include "./yarns.h"

namespace file_format {

YarnRepr::YarnRepr(file_format::Yarns::Yarns yarns) {
  for (auto& yarn : yarns.yarns) {
    int nPoints = yarn.points.size();

    // Construct point list
    Eigen::MatrixXf P(nPoints, 3);
    for (int i = 0; i < nPoints; i++) {
      P(i, 0) = yarn.points[i][0];
      P(i, 1) = yarn.points[i][1];
      P(i, 2) = yarn.points[i][2];
    }

    // Other meta data
    Yarn internalYarn;
    internalYarn.points = P;
    internalYarn.color = yarn.color;
    internalYarn.radius = yarn.radius;

    // Save
    this->yarns.push_back(internalYarn);
  }
}

YarnRepr YarnRepr::createAlike() const {
  YarnRepr res;
  res.yarns.reserve(yarns.size());
  for (const Yarn& yarn : yarns) {
    Yarn internalYarn;
    internalYarn.color = yarn.color;
    internalYarn.radius = yarn.radius;
    res.yarns.push_back(internalYarn);
  }
  return res;
}

}  // namespace file_format
