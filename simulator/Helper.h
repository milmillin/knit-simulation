
#include <functional>
#include <Eigen/Core>

namespace simulator {

constexpr float SIMPSON_EPS = 1e-6;

// Performs Simpson's Integration of function f over [a, b]
float integrate(const std::function<float(float)>& f, float a, float b);

// Calculates Arc length of a Catmull-Rom Spline defined by 4 control points.
// 
// q : flattened coordinates #m x 1
// index : index of x-coord of the first control point in q
float catmullRomArcLength(const Eigen::MatrixXf& q, int index);

}; // namespace simulator