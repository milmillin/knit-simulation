#include <glm/glm.hpp>
#include <glm/ext/scalar_constants.hpp>

#include "sweep.h"

#define POINT_FROM_ROW(matrix, index) glm::vec3(matrix((index), 0), matrix((index), 1), matrix((index), 2))
#define ROW_FROM_POINT(matrix, index, point) matrix((index), 0) = point.x; matrix((index), 1) = point.y; matrix((index), 2) = point.z;

static Eigen::MatrixXf circle(glm::vec3 offset, glm::vec3 d1, glm::vec3 d2, float radius, int stride) {
  Eigen::MatrixXf vertices(stride, 3);
  for (int i = 0; i < stride; i++) {
    glm::vec3 vertex = glm::sin((float)i * 2 * glm::pi<float>() / stride) * d1 
                     + glm::cos((float)i * 2 * glm::pi<float>() / stride) * d2
                     + offset;
    vertices(i, 0) = vertex.x;
    vertices(i, 1) = vertex.y;
    vertices(i, 2) = vertex.z;
  }
  return vertices;
}

static Eigen::MatrixXf circleEnd(glm::vec3 offset, glm::vec3 normal, float radius, int stride) {
  glm::vec3 d1;
  if (glm::abs(normal.x + normal.y) < glm::abs(normal.z)) {
    d1 = glm::cross(glm::vec3(1, 0, 0), normal);
  } else {
    d1 = glm::cross(glm::vec3(0, 0, 1), normal);
  }
  glm::vec3 d2 = glm::cross(normal, d1);

  d1 = glm::normalize(d1) * radius;
  d2 = glm::normalize(d2) * radius;

  return circle(offset, d1, d2, radius, stride);
}

static Eigen::MatrixXf circleJoint(glm::vec3 offset, glm::vec3 dir1, glm::vec3 dir2, float radius, int stride) {
  dir1 = glm::normalize(dir1);
  dir2 = glm::normalize(dir2);

  if (glm::length(dir1 + dir2) < 0.001) {
    return circleEnd(offset, dir1, radius, stride);
  }

  glm::vec3 d1 = glm::normalize(dir1 + dir2);
  glm::vec3 d2 = glm::cross(dir1, dir2);
  return circle(offset, d1, d2, radius, stride);
}

void circleSweep(const Eigen::MatrixXf path, const float radius,
    Eigen::MatrixXf &vertex, Eigen::MatrixXi &triangles,
    const int stride) {
  // const int nPoints = path.rows();
  const int nPoints = 4;

  if (nPoints <= 1) {
    vertex = Eigen::MatrixXf(0, 3);
    triangles = Eigen::MatrixXi(0, 3);
    return;
  }

  vertex = Eigen::MatrixXf(nPoints * stride + 2, 3);
  triangles = Eigen::MatrixXi((nPoints - 1) * stride * 2, 3);

  glm::vec3 startPoint = POINT_FROM_ROW(path, 0);
  glm::vec3 nextPoint = POINT_FROM_ROW(path, 1);
  vertex.block(0, 0, stride, 3) = \
    circleEnd(startPoint, nextPoint - startPoint, radius, stride);
  
  for (int i = 1; i < nPoints - 1; i++) {
    glm::vec3 prevPoint = POINT_FROM_ROW(path, i - 1);
    glm::vec3 currentPoint = POINT_FROM_ROW(path, i);
    glm::vec3 nextPoint = POINT_FROM_ROW(path, i + 1);

    vertex.block(i * stride, 0, stride, 3) = 
      circleJoint(currentPoint, prevPoint - currentPoint, nextPoint - currentPoint,
                  radius, stride);
  }
  glm::vec3 endPoint = POINT_FROM_ROW(path, nPoints - 1);
  glm::vec3 prevPoint = POINT_FROM_ROW(path, nPoints - 2);
  vertex.block((nPoints - 1) * stride, 0, stride, 3) = \
    circleEnd(endPoint, endPoint - prevPoint, radius, stride);
  
  ROW_FROM_POINT(vertex, vertex.rows() - 2, startPoint);
  ROW_FROM_POINT(vertex, vertex.rows() - 1, endPoint);

  int index = 0;
  for (int i = 1; i < nPoints; i++) {
    for (int j = 0; j < stride; j++) {
      triangles(index, 0) = (i-1)*stride + j;
      triangles(index, 1) = (i-1)*stride + (j + 1) % stride;
      triangles(index, 2) = i*stride + j;
      index++;
      triangles(index, 0) = i*stride + j;
      triangles(index, 1) = (i-1)*stride + (j + 1) % stride;
      triangles(index, 2) = i*stride + (j + 1) % stride;
      index++;
    }
  }

  // for (int j = 0; j < stride; j++) {
  //     triangles(index, 0) = j;
  //     triangles(index, 1) = (j + 1) % stride;
  //     triangles(index, 2) = nPoints*stride*2;
  //     index++;
  //     triangles(index, 0) = (nPoints - 1) * stride + j;
  //     triangles(index, 1) = (nPoints - 1) * stride + (j + 1) % stride;
  //     triangles(index, 2) = nPoints*stride*2 + 1;
  //     index++;
  // }
}