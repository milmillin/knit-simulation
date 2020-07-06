#include "viewer.h"
#include "sweep.h"
#include "../simulator/SimulatorParams.h"
#include "../file_format/yarns.h"

namespace UI {

int Viewer::launch(bool resizable, bool fullscreen, const std::string &name, int width, int height) {
  // Add menu
  _menu.reset(new Menu());
  this->plugins.push_back(_menu.get());

  // Start trackball mode (allow panning)
  this->core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

  // Disable wireframe display
  this->data().show_lines = 0;

  // Load yarns
  refresh();

  // Launch
  return igl::opengl::glfw::Viewer::launch(resizable, fullscreen, name, width, height);
}

void Viewer::refresh() {
  // Clear old mesh
  this->data().clear();

  // Get yarn shape
  Eigen::MatrixXf points = _simulator.getControlPoints();

  // Create mesh
  // TODO: don't hard-code radius
  Eigen::MatrixXf vertices;
  Eigen::MatrixXi triangles;
  circleSweep(points, 0.1, vertices, triangles, 8);
  this->data().set_mesh(vertices.cast<double>(), triangles);

  // Draw line
  if (points.rows() >= 2) {
    Eigen::MatrixXi E(points.rows() - 1, 2);
    for (int i = 0; i < points.rows() - 1; i++) {
      E(i, 0) = i;
      E(i, 1) = i + 1;
    }
    this->data().set_edges(points.cast<double>(),
        E, Eigen::RowVector3d(1, 1, 1));
  }
}

void Viewer::loadYarn(std::string filename) {
  // Load .yarns file
  std::cout << "Loading model: " << filename << std::endl;
  file_format::Yarns yarn;
  try {
    yarn = file_format::Yarns::load(filename);
  } catch (const std::runtime_error& e) {
    std::cout << "Failed to load " << filename << std::endl;
    std::cout << e.what() << std::endl;
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

  // Update simulator
  _simulator = simulator::Simulator(P, simulator::SimulatorParams::Default());
  this->refresh();
}

} // namespace UI