#include <igl/opengl/glfw/Viewer.h>

#include "sm.hpp"
#include "Simulator.h"
#include "SimulatorParams.h"

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  std::cout << "Hello, world!" << std::endl;

  // Load .yarns file
  sm::Yarns yarn;
  try {
    yarn = sm::Yarns::load("../helloworld.yarns");
  } catch (const std::runtime_error& e) {
    std::cout << e.what() << std::endl;
    return 0;
  }

  // Construct polyline
  int n = yarn.yarns[0].points.size();
  Eigen::MatrixXf P(n, 3);
  for (int i = 0; i < n; i++) {
    P(i, 0) = yarn.yarns[0].points[i][0];
    P(i, 1) = yarn.yarns[0].points[i][1];
    P(i, 2) = yarn.yarns[0].points[i][2];
  }

  Eigen::MatrixXi E(n - 1, 2);
  for (int i = 0; i < n - 1; i++) {
    E(i, 0) = i;
    E(i, 1) = i + 1;
  }

  // Feed to the simulator
  Simulator simulator(P, SimulatorParams::Default());
  simulator.step();

  // TODO: Construct tube mesh from polyline

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_edges(simulator.getCurrentPoints().cast<double>(),
      E, Eigen::RowVector3d(1, 1, 1));
  viewer.launch();
}
