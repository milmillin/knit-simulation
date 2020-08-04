#pragma once

#define DECLARE_POINTS2(p, q, index) \
  const Eigen::Vector3d p[4] = { \
    Eigen::Vector3d(q.block<3, 1>((index), 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 3, 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 6, 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 9, 0)) }

#define DECLARE_BASIS2(b, s) \
  const double b[4] = { \
    (s) * (-1.0 / 2.0) + (s) * (s) - ((s) * (s) * (s)) / 2.0, \
    ((s) * (s)) * (-5.0 / 2.0) + ((s) * (s) * (s)) * (3.0 / 2.0) + 1.0, \
    (s) / 2.0 + ((s) * (s)) * 2.0 - ((s) * (s) * (s)) * (3.0 / 2.0), \
    ((s) * (s)) * (-1.0 / 2.0) + ((s) * (s) * (s)) / 2.0 }

#define DECLARE_BASIS_D2(b, s) \
  const double b[4] = { \
    (s) * 2.0 - ((s) * (s)) * (3.0 / 2.0) - 1.0 / 2.0, \
    (s) * -5.0 + ((s) * (s)) * (9.0 / 2.0), \
    (s) * 4.0 - ((s) * (s)) * (9.0 / 2.0) + 1.0 / 2.0, \
    -(s) + ((s) * (s)) * (3.0 / 2.0) }

#define DECLARE_BASIS_DD2(b, s) \
  const double b[4] = { \
  (s) * -3.0 + 2.0, \
  (s) * 9.0 - 5.0, \
  (s) * -9.0 + 4.0, \
  (s) * 3.0 - 1.0 }


#define POINT_FROM_BASIS(points, basis) \
  (basis[0] * points[0] + basis[1] * points[1] + basis[2] * points[2] + basis[3] * points[3])

// Get a point from a row in the matrix
#define POINT_FROM_ROW(matrix, index) \
  glm::dvec3(matrix((index), 0), matrix((index), 1), matrix((index), 2))

// Set a row in the matrix with a point
#define ROW_FROM_POINT(matrix, index, point) \
  (matrix)((index), 0) = point.x; \
  (matrix)((index), 1) = point.y; \
  (matrix)((index), 2) = point.z;

// Add a row in the matrix by a vector
#define ADD_TO_ROW(matrix, index, point) \
  (matrix)((index), 0) += point.x; \
  (matrix)((index), 1) += point.y; \
  (matrix)((index), 2) += point.z;

// Subtract a row in the matrix by a vector
#define SUBTRACT_FROM_ROW(matrix, index, point) \
  (matrix)((index), 0) -= point.x; \
  (matrix)((index), 1) -= point.y; \
  (matrix)((index), 2) -= point.z;

