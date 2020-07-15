
#include <functional>
#include <Eigen/Core>

#include <vector>

#include <glm/glm.hpp>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-6;

// Performs Simpson's Integration of function f over [a, b]
float integrate(const std::function<float(float)>& f, float a, float b);



// Catmull-Rom coefficient
inline float b1(float s) {
  return -s / 2 + s * s - s * s * s / 2;
}
inline float b2(float s) {
  return 1 - s * s * 5 / 2 + s * s * s * 3 / 2;
}
inline float b3(float s) {
  return s / 2 + 2 * s * s - s * s * s * 3 / 2;
}
inline float b4(float s) {
  return -s * s / 2 + s * s * s / 2;
}

// Calculates Arc length of a Catmull-Rom Spline defined by 4 control points.
// 
// q : flattened coordinates #m x 1
// index : index of x-coord of the first control point in q
float catmullRomArcLength(const Eigen::MatrixXf& q, int index);


glm::vec3 catmullRomSample(const Eigen::MatrixXf &controlPoints, int index, float s);

// Sample a Catmul-Rom curve
//
// points: one row for each point coordinate
// samplePerSegment: number of samples for each segment
// Return: samples (one row for each point coordinate)
Eigen::MatrixXf catmullRomSequenceSample(Eigen::MatrixXf points, int samplePerSegment);

void catmullRomBoundingBox(const Eigen::MatrixXf &controlPoints, int index,
  std::vector<double> *lowerBound, std::vector<double> *upperBound, float radius);

}; // namespace simulator
