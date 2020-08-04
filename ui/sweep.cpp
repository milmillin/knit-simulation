#include "./sweep.h"

#include <glm/glm.hpp>
#include <glm/ext/scalar_constants.hpp>

// Get a point from a row in the matrix
#define POINT_FROM_ROW(matrix, index) \
  glm::dvec3(matrix((index), 0), matrix((index), 1), matrix((index), 2))
// Set a row in the matrix with a point
#define ROW_FROM_POINT(matrix, index, point) \
  (matrix)((index), 0) = point.x; \
  (matrix)((index), 1) = point.y; \
  (matrix)((index), 2) = point.z;

namespace UI {

// Return a list of 3D points sampled on a 2D circle
// offset: center of circle
// d1: one of two vectors defining the plane
// d2: one of two vectors defining the plane
// radius: circle radius
// stride: number of samples
inline static Eigen::MatrixXd circle(glm::dvec3 offset, glm::dvec3 d1, glm::dvec3 d2,
    double radius, int stride) {
  Eigen::MatrixXd vertices(stride, 3);
  for (int i = 0; i < stride; i++) {
    glm::dvec3 vertex = radius *
      (glm::sin((double)i * 2 * glm::pi<double>() / stride) * d1
        + glm::cos((double)i * 2 * glm::pi<double>() / stride) * d2)
      + offset;
    vertices(i, 0) = vertex.x;
    vertices(i, 1) = vertex.y;
    vertices(i, 2) = vertex.z;
  }
  return vertices;
}

// Calculate next parallel transform frame
// tangent: new tangent direction
// oldUp: old up direction
// newUp: new up direction
// newCross: new cross direction
inline static void parallelTransform(glm::dvec3 tangent, glm::dvec3 oldUp,
    glm::dvec3 *newUp, glm::dvec3 *newCross) {
  *newCross = glm::normalize(glm::cross(tangent, oldUp));
  *newUp = glm::normalize(glm::cross(*newCross, tangent));
}

void circleSweep(const Eigen::MatrixXd path, const double radius, const int stride,
    Eigen::MatrixXd *vertices, Eigen::MatrixXi *triangles) {
  // Number of points on the path
  const int nPoints = path.rows();

  // The path is empty
  if (nPoints < 2) {
    vertices->resize(0, 3);
    triangles->resize(0, 3);
    return;
  }

  // Allocate memory
  vertices->resize(nPoints * stride + 2, 3);
  triangles->resize(nPoints * stride * 2, 3);

  // Initial frame
  glm::dvec3 tangent;
  glm::dvec3 up;
  glm::dvec3 cross;

  tangent = POINT_FROM_ROW(path, 1) - POINT_FROM_ROW(path, 0);

  // Use an arbitrary vector that's normal to tangent line
  if (glm::abs(tangent.x + tangent.y) < glm::abs(tangent.z)) {
    up = glm::cross(glm::dvec3(1, 0, 0), tangent);
  } else {
    up = glm::cross(glm::dvec3(0, 0, 1), tangent);
  }
  up = glm::normalize(up);

  cross = glm::normalize(glm::cross(tangent, up));

  // Intial circle
  vertices->block(0, 0, stride, 3) = \
    circle(POINT_FROM_ROW(path, 0), up, cross, radius, stride);

  for (int i = 1; i < nPoints; i++) {
    // Update tangent
    glm::dvec3 prevPoint;
    glm::dvec3 currentPoint;
    glm::dvec3 nextPoint;

    prevPoint = POINT_FROM_ROW(path, i - 1);
    currentPoint = POINT_FROM_ROW(path, i);
    glm::dvec3 prevDir = currentPoint - prevPoint;
    if (i == nPoints - 1) {
      // This is the last point
      tangent = prevDir;
    } else {
      // Take the average of two segments
      nextPoint = POINT_FROM_ROW(path, i + 1);
      glm::dvec3 nextDir = nextPoint - currentPoint;
      tangent = glm::normalize(glm::normalize(prevDir)+glm::normalize(nextDir));
    }

    // Get next frame
    parallelTransform(tangent, up, &up, &cross);

    // Get circle
    vertices->block(i * stride, 0, stride, 3) =
      circle(currentPoint, up, cross, radius, stride);
  }
  glm::dvec3 endPoint = POINT_FROM_ROW(path, nPoints - 1);

  // The last two vertices is the start and end point
  ROW_FROM_POINT(*vertices, vertices->rows() - 2, POINT_FROM_ROW(path, 0));
  ROW_FROM_POINT(*vertices, vertices->rows() - 1, POINT_FROM_ROW(path, nPoints - 1));

  // Create faces for the side of tube
  int index = 0;
  for (int i = 1; i < nPoints; i++) {
    for (int j = 0; j < stride; j++) {
      (*triangles)(index, 0) = (i-1)*stride + j;
      (*triangles)(index, 1) = i*stride + j;
      (*triangles)(index, 2) = (i-1)*stride + (j + 1) % stride;
      index++;
      (*triangles)(index, 0) = i*stride + j;
      (*triangles)(index, 1) = i*stride + (j + 1) % stride;
      (*triangles)(index, 2) = (i-1)*stride + (j + 1) % stride;
      index++;
    }
  }

  // Create faces for the end of tube
  for (int j = 0; j < stride; j++) {
    // Starting end
    (*triangles)(index, 0) = j;
    (*triangles)(index, 1) = (j + 1) % stride;
    (*triangles)(index, 2) = nPoints*stride;
    index++;
    // Termnating end
    (*triangles)(index, 0) = (nPoints - 1) * stride + (j + 1) % stride;
    (*triangles)(index, 1) = (nPoints - 1) * stride + j;
    (*triangles)(index, 2) = nPoints*stride + 1;
    index++;
  }
}

}  // namespace UI
