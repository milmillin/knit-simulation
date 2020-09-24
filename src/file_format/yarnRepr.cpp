#include "./yarnRepr.h"

#include <vector>

#include <Eigen/Core>
#include <glm/common.hpp>

#include "./yarns.h"

namespace file_format {

YarnRepr::YarnRepr(file_format::Yarns::Yarns& yarns) {
  size_t n = 0;
  for (auto& yarn : yarns.yarns) {
    n += yarn.points.size();
  }
  vertices = Eigen::MatrixXd(n, 3);
  this->yarns.reserve(yarns.yarns.size());

  size_t idx = 0;
  for (auto& yarn : yarns.yarns) {
    size_t nPoints = yarn.points.size();

    this->yarns.emplace_back();
    YarnInfo& yarnInfo = this->yarns.back();
    yarnInfo.begin = idx;
    yarnInfo.end = idx + nPoints;
    yarnInfo.color = Eigen::RowVector3d(
      yarn.color[0] / 255., yarn.color[1] / 255., yarn.color[2] / 255.);
    yarnInfo.radius = yarn.radius;

    // Construct point list
    for (size_t i = 0; i < nPoints; i++) {
      vertices(idx, 0) = yarn.points[i][0];
      vertices(idx, 1) = yarn.points[i][1];
      vertices(idx, 2) = yarn.points[i][2];
      idx++;
    }
    // Initialize frames to zero
    this->bishopFrameU = Eigen::MatrixX3d::Zero(n, 3);
    this->bishopFrameU = Eigen::MatrixX3d::Zero(n, 3);
    this->materialFrameU = Eigen::MatrixX3d::Zero(n, 3);
    this->materialFrameV = Eigen::MatrixX3d::Zero(n, 3);
  }
}

Yarns::Yarns YarnRepr::toYarns() const {
  Yarns::Yarns res;
  res.yarns.reserve(yarns.size());
  for (const YarnInfo& yarn : yarns) {
    res.yarns.emplace_back();
    Yarns::Yarn& thatYarn = res.yarns.back();
    thatYarn.color = glm::u8vec4(
      yarn.color(0) * 255,
      yarn.color(1) * 255,
      yarn.color(2) * 255, 1);
    thatYarn.radius = yarn.radius;

    size_t n = yarn.end - yarn.begin;
    thatYarn.points.resize(yarn.end - yarn.begin);
    for (size_t i = 0; i < n; i++) {
      thatYarn.points[i] = glm::vec3(
        vertices(yarn.begin + i, 0),
        vertices(yarn.begin + i, 1),
        vertices(yarn.begin + i, 2));
    }
    thatYarn.sources.resize(n, 0);
  }
  return res;
}

YarnRepr YarnRepr::createAlike() const {
  YarnRepr res;
  res.yarns.reserve(yarns.size());
  for (const YarnInfo& yarn : yarns) {
    res.yarns.push_back(yarn);
  }
  return res;
}

Eigen::Block<Eigen::MatrixXd> YarnRepr::getYarnPoints(size_t index) {
  assert(index < yarns.size());
  const auto& yarn = yarns[index];
  return vertices.block(yarn.begin, 0, yarn.end - yarn.begin, 3);
}

const Eigen::Block<const Eigen::MatrixXd> YarnRepr::getYarnPoints(size_t index) const {
  assert(index < yarns.size());
  const auto& yarn = yarns[index];
  return vertices.block(yarn.begin, 0, yarn.end - yarn.begin, 3);
}

}  // namespace file_format
