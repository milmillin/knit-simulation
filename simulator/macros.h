#pragma once

#define DECLARE_POINTS(q, index) \
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

#define DECLARE_BASIS \
  float b1 = s * (-1.0 / 2.0) + s * s - (s * s * s) / 2.0; \
  float b2 = (s * s) * (-5.0 / 2.0) + (s * s * s) * (3.0 / 2.0) + 1.0; \
  float b3 = s / 2.0 + (s * s) * 2.0 - (s * s * s) * (3.0 / 2.0); \
  float b4 (s * s) * (-1.0 / 2.0) + (s * s * s) / 2.0;

#define DECLARE_BASIS_D \
  float bD1 = s * 2.0 - (s * s) * (3.0 / 2.0) - 1.0 / 2.0; \
  float bD2 = s * -5.0 + (s * s) * (9.0 / 2.0); \
  float bD3 = s * 4.0 - (s * s) * (9.0 / 2.0) + 1.0 / 2.0; \
  float bD4 = -s + (s * s) * (3.0 / 2.0);

#define DECLARE_BASIS_DD \
  float bDD1 = s * -3.0 + 2.0; \
  float bDD2 = s * 9.0 - 5.0; \
  float bDD3 = s * -9.0 + 4.0; \
  float bDD4 = s * 3.0 - 1.0;


