#include <Eigen/Core>

#include "file_format/yarns.h"
#include "simulator/Simulator.h"
#include "simulator/SimulatorParams.h"
#include "ui/viewer.h"

int main(int argc, char *argv[]) {

  // Load .yarns file
  std::cout << "Loading model" << std::endl;
  file_format::Yarns yarn;
  try {
    yarn = file_format::Yarns::load("../helloworld.yarns");
  } catch (const std::runtime_error& e) {
    std::cout << e.what() << std::endl;
    return 0;
  }

  // Construct polyling
  std::cout << "Constructing line" << std::endl;
  int n = yarn.yarns[0].points.size();
  Eigen::MatrixXf P(n, 3);
  for (int i = 0; i < n; i++) {
    P(i, 0) = yarn.yarns[0].points[i][0];
    P(i, 1) = yarn.yarns[0].points[i][1];
    P(i, 2) = yarn.yarns[0].points[i][2];
  }

  // Feed to the simulator
  std::cout << "Simulating" << std::endl;
  simulator::Simulator simulator(P, simulator::SimulatorParams::Default());
  simulator.step();

  // Launch viewer
  std::cout << "Launching viewer" << std::endl;
  UI::Viewer viewer;
  viewer.plot(simulator.getControlPoints());
  viewer.launch();
}
