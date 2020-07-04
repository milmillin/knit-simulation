#include <Eigen/Core>

#include "file_format/yarns.h"
#include "simulator/Simulator.h"
#include "simulator/SimulatorParams.h"
#include "ui/viewer.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " yarns-filename" << std::endl;
    return -1;
  }

  // Load .yarns file
  std::cout << "Loading model: " << argv[1] << std::endl;
  file_format::Yarns yarn;
  try {
    yarn = file_format::Yarns::load(argv[1]);
  } catch (const std::runtime_error& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }

  // Construct polyline
  std::cout << "Constructing line" << std::endl;
  int n = yarn.yarns[0].points.size();
  Eigen::MatrixXf P(n, 3);
  for (int i = 0; i < n; i++) {
    P(i, 0) = yarn.yarns[0].points[i][0];
    P(i, 1) = yarn.yarns[0].points[i][1];
    P(i, 2) = yarn.yarns[0].points[i][2];
  }

  std::cout << "Found: " << std::to_string(n) << " control points." << std::endl;

  // Feed to the simulator
  std::cout << "Initializing Simulator" << std::endl;
  simulator::Simulator simulator(P, simulator::SimulatorParams::Default());
  simulator.step();

  // Launch viewer
  std::cout << "Launching viewer" << std::endl;
  UI::Viewer viewer;
  viewer.plot(simulator.getControlPoints());
  viewer.launch();
}
