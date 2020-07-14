
#include <functional>
#include <Eigen/Core>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-5;

// Performs adaptive Simpson's integration of function f over [a, b]
// 
// subdivide : minimum subdivision
float integrate(const std::function<float(float)>& f, float a, float b, int subdivide = 1, int maxDepth = 16);

// Calculates Arc length of a Catmull-Rom Spline defined by 4 control points.
// 
// q : flattened coordinates #m x 1
// index : index of x-coord of the first control point in q
float catmullRomArcLength(const Eigen::MatrixXf& q, int index);

// Sample a Catmul-Rom curve
//
// points: one row for each point coordinate
// samplePerSegment: number of samples for each segment
// Return: samples (one row for each point coordinate)
Eigen::MatrixXf catmullRomSequenceSample(Eigen::MatrixXf points, int samplePerSegment);

}; // namespace simulator
