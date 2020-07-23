#pragma once

#define DECLARE_POINTS(p, index) \
  const float& p##x1 = q((index)); \
  const float& p##y1 = q((index) + 1); \
  const float& p##z1 = q((index) + 2); \
  const float& p##x2 = q((index) + 3); \
  const float& p##y2 = q((index) + 4); \
  const float& p##z2 = q((index) + 5); \
  const float& p##x3 = q((index) + 6); \
  const float& p##y3 = q((index) + 7); \
  const float& p##z3 = q((index) + 8); \
  const float& p##x4 = q((index) + 9); \
  const float& p##y4 = q((index) + 10); \
  const float& p##z4 = q((index) + 11);

#define DECLARE_POINTS_D(p, index) \
  const float& p##xD1 = qD((index)); \
  const float& p##yD1 = qD((index) + 1); \
  const float& p##zD1 = qD((index) + 2); \
  const float& p##xD2 = qD((index) + 3); \
  const float& p##yD2 = qD((index) + 4); \
  const float& p##zD2 = qD((index) + 5); \
  const float& p##xD3 = qD((index) + 6); \
  const float& p##yD3 = qD((index) + 7); \
  const float& p##zD3 = qD((index) + 8); \
  const float& p##xD4 = qD((index) + 9); \
  const float& p##yD4 = qD((index) + 10); \
  const float& p##zD4 = qD((index) + 11);

#define DECLARE_POINTS2(p, q, index) \
  const Eigen::Vector3f p[4] = { \
    Eigen::Vector3f(q.block<3, 1>((index), 0)), \
    Eigen::Vector3f(q.block<3, 1>((index) + 3, 0)), \
    Eigen::Vector3f(q.block<3, 1>((index) + 6, 0)), \
    Eigen::Vector3f(q.block<3, 1>((index) + 9, 0)) }

#define DECLARE_BASIS(b, s) \
  float b##1 = s * (-1.0f / 2.0f) + s * s - (s * s * s) / 2.0f; \
  float b##2 = (s * s) * (-5.0f / 2.0f) + (s * s * s) * (3.0f / 2.0f) + 1.0f; \
  float b##3 = s / 2.0f + (s * s) * 2.0f - (s * s * s) * (3.0f / 2.0f); \
  float b##4 = (s * s) * (-1.0f / 2.0f) + (s * s * s) / 2.0f;

#define DECLARE_BASIS_D(b, s) \
  float b##D1 = s * 2.0f - (s * s) * (3.0f / 2.0f) - 1.0f / 2.0f; \
  float b##D2 = s * -5.0f + (s * s) * (9.0f / 2.0f); \
  float b##D3 = s * 4.0f - (s * s) * (9.0f / 2.0f) + 1.0f / 2.0f; \
  float b##D4 = -s + (s * s) * (3.0f / 2.0f);

#define DECLARE_BASIS_DD(b, s) \
  float b##DD1 = s * -3.0f + 2.0f; \
  float b##DD2 = s * 9.0f - 5.0f; \
  float b##DD3 = s * -9.0f + 4.0f; \
  float b##DD4 = s * 3.0f - 1.0f;

#define DECLARE_BASIS2(b, s) \
  const float b[4] = { \
    (s) * (-1.0f / 2.0f) + (s) * (s) - ((s) * (s) * (s)) / 2.0f, \
    ((s) * (s)) * (-5.0f / 2.0f) + ((s) * (s) * (s)) * (3.0f / 2.0f) + 1.0f, \
    (s) / 2.0f + ((s) * (s)) * 2.0f - ((s) * (s) * (s)) * (3.0f / 2.0f), \
    ((s) * (s)) * (-1.0f / 2.0f) + ((s) * (s) * (s)) / 2.0f }

#define DECLARE_BASIS_D2(b, s) \
  const float b[4] = { \
    (s) * 2.0f - ((s) * (s)) * (3.0f / 2.0f) - 1.0f / 2.0f, \
    (s) * -5.0f + ((s) * (s)) * (9.0f / 2.0f), \
    (s) * 4.0f - ((s) * (s)) * (9.0f / 2.0f) + 1.0f / 2.0f, \
    -(s) + ((s) * (s)) * (3.0f / 2.0f) }

#define DECLARE_BASIS_DD2(b, s) \
  const float b[4] = { \
  (s) * -3.0f + 2.0f, \
  (s) * 9.0f - 5.0f, \
  (s) * -9.0f + 4.0f, \
  (s) * 3.0f - 1.0f }


#define POINT_FROM_BASIS(points, basis) \
  (basis[0] * points[0] + basis[1] * points[1] + basis[2] * points[2] + basis[3] * points[3])

// Get a point from a row in the matrix
#define POINT_FROM_ROW(matrix, index) \
  glm::vec3(matrix((index), 0), matrix((index), 1), matrix((index), 2))

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

