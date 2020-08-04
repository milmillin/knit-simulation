#pragma once

#define DECLARE_POINTS(p, index) \
  const double& p##x1 = q((index)); \
  const double& p##y1 = q((index) + 1); \
  const double& p##z1 = q((index) + 2); \
  const double& p##x2 = q((index) + 3); \
  const double& p##y2 = q((index) + 4); \
  const double& p##z2 = q((index) + 5); \
  const double& p##x3 = q((index) + 6); \
  const double& p##y3 = q((index) + 7); \
  const double& p##z3 = q((index) + 8); \
  const double& p##x4 = q((index) + 9); \
  const double& p##y4 = q((index) + 10); \
  const double& p##z4 = q((index) + 11);

#define DECLARE_POINTS_D(p, index) \
  const double& p##xD1 = qD((index)); \
  const double& p##yD1 = qD((index) + 1); \
  const double& p##zD1 = qD((index) + 2); \
  const double& p##xD2 = qD((index) + 3); \
  const double& p##yD2 = qD((index) + 4); \
  const double& p##zD2 = qD((index) + 5); \
  const double& p##xD3 = qD((index) + 6); \
  const double& p##yD3 = qD((index) + 7); \
  const double& p##zD3 = qD((index) + 8); \
  const double& p##xD4 = qD((index) + 9); \
  const double& p##yD4 = qD((index) + 10); \
  const double& p##zD4 = qD((index) + 11);

#define DECLARE_POINTS2(p, q, index) \
  const Eigen::Vector3d p[4] = { \
    Eigen::Vector3d(q.block<3, 1>((index), 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 3, 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 6, 0)), \
    Eigen::Vector3d(q.block<3, 1>((index) + 9, 0)) }

#define DECLARE_BASIS(b, s) \
  double b##1 = s * (-1.0 / 2.0) + s * s - (s * s * s) / 2.0; \
  double b##2 = (s * s) * (-5.0 / 2.0) + (s * s * s) * (3.0 / 2.0) + 1.0; \
  double b##3 = s / 2.0 + (s * s) * 2.0 - (s * s * s) * (3.0 / 2.0); \
  double b##4 = (s * s) * (-1.0 / 2.0) + (s * s * s) / 2.0;

#define DECLARE_BASIS_D(b, s) \
  double b##D1 = s * 2.0 - (s * s) * (3.0 / 2.0) - 1.0 / 2.0; \
  double b##D2 = s * -5.0 + (s * s) * (9.0 / 2.0); \
  double b##D3 = s * 4.0 - (s * s) * (9.0 / 2.0) + 1.0 / 2.0; \
  double b##D4 = -s + (s * s) * (3.0 / 2.0);

#define DECLARE_BASIS_DD(b, s) \
  double b##DD1 = s * -3.0 + 2.0; \
  double b##DD2 = s * 9.0 - 5.0; \
  double b##DD3 = s * -9.0 + 4.0; \
  double b##DD4 = s * 3.0 - 1.0;

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

