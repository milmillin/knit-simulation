#pragma once

#define DECLARE_POINTS(index) \
  const float& px1 = q((index)); \
  const float& py1 = q((index) + 1); \
  const float& pz1 = q((index) + 2); \
  const float& px2 = q((index) + 3); \
  const float& py2 = q((index) + 4); \
  const float& pz2 = q((index) + 5); \
  const float& px3 = q((index) + 6); \
  const float& py3 = q((index) + 7); \
  const float& pz3 = q((index) + 8); \
  const float& px4 = q((index) + 9); \
  const float& py4 = q((index) + 10); \
  const float& pz4 = q((index) + 11);

#define DECLARE_POINTS_D(index) \
  const float& pxD1 = qD((index)); \
  const float& pyD1 = qD((index) + 1); \
  const float& pzD1 = qD((index) + 2); \
  const float& pxD2 = qD((index) + 3); \
  const float& pyD2 = qD((index) + 4); \
  const float& pzD2 = qD((index) + 5); \
  const float& pxD3 = qD((index) + 6); \
  const float& pyD3 = qD((index) + 7); \
  const float& pzD3 = qD((index) + 8); \
  const float& pxD4 = qD((index) + 9); \
  const float& pyD4 = qD((index) + 10); \
  const float& pzD4 = qD((index) + 11);

#define DECLARE_BASIS \
  float b1 = s * (-1.0f / 2.0f) + s * s - (s * s * s) / 2.0f; \
  float b2 = (s * s) * (-5.0f / 2.0f) + (s * s * s) * (3.0f / 2.0f) + 1.0f; \
  float b3 = s / 2.0f + (s * s) * 2.0f - (s * s * s) * (3.0f / 2.0f); \
  float b4 = (s * s) * (-1.0f / 2.0f) + (s * s * s) / 2.0f;

#define DECLARE_BASIS_D \
  float bD1 = s * 2.0f - (s * s) * (3.0f / 2.0f) - 1.0f / 2.0f; \
  float bD2 = s * -5.0f + (s * s) * (9.0f / 2.0f); \
  float bD3 = s * 4.0f - (s * s) * (9.0f / 2.0f) + 1.0f / 2.0f; \
  float bD4 = -s + (s * s) * (3.0f / 2.0f);

#define DECLARE_BASIS_DD \
  float bDD1 = s * -3.0f + 2.0f; \
  float bDD2 = s * 9.0f - 5.0f; \
  float bDD3 = s * -9.0f + 4.0f; \
  float bDD4 = s * 3.0f - 1.0f;

// Get a point from a row in the matrix
#define POINT_FROM_ROW(matrix, index) \
  glm::vec3(matrix((index), 0), matrix((index), 1), matrix((index), 2))

// Set a row in the matrix with a point
#define ROW_FROM_POINT(matrix, index, point) \
  (matrix)((index), 0) = point.x; \
  (matrix)((index), 1) = point.y; \
  (matrix)((index), 2) = point.z;
