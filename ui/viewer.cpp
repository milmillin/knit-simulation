#include "./viewer.h"

#include "./sweep.h"
#include "../simulator/SimulatorParams.h"
#include "../file_format/yarns.h"
#include "../file_format/yarnRepr.h"
#include "../simulator/Helper.h"

namespace UI {

int Viewer::launch(bool resizable, bool fullscreen, const std::string &name, int width, int height) {
  // Add menu
  _menu.reset(new Menu());
  this->plugins.push_back(_menu.get());

  // Start trackball mode (allow panning)
  this->core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

  // Disable wireframe display
  this->data().show_lines = 0;

  // Set background color
  this->core().background_color.data()[0] = 0;
  this->core().background_color.data()[1] = 0;
  this->core().background_color.data()[2] = 0;

  // Load yarns
  refresh();

  // Launch
  return igl::opengl::glfw::Viewer::launch(resizable, fullscreen, name, width, height);
}

void Viewer::refresh() {

  // Get yarn shape
  const file_format::YarnRepr yarns = simulator.getYarns();

  this->selected_data_index = 0;
  Eigen::MatrixXf groundPoints(4, 3);
  groundPoints << 10, simulator.params.groundHeight, 10,
    10, simulator.params.groundHeight, -10,
    -10, simulator.params.groundHeight, -10,
    -10, simulator.params.groundHeight, 10;
  Eigen::MatrixXi groundTrianges(2, 3);
  groundTrianges << 2, 0, 1,
    3, 0, 2;
  this->data().set_mesh(groundPoints.cast<double>(), groundTrianges);
  this->data_list.push_back(igl::opengl::ViewerData());

  for (int i = 0; i < yarns.yarns.size(); i++) {
    // Clear old mesh
    this->selected_data_index = i + 1;
    this->data().clear();

    // Get curve
    auto &yarn = yarns.yarns[i];

    // Use Catmull-Rom to smooth the curve
    Eigen::MatrixXf points = simulator::catmullRomSequenceSample(yarn.points, curveSamples);

    // Create mesh for tube
    Eigen::MatrixXf vertices;
    Eigen::MatrixXi triangles;
    circleSweep(points, yarn.radius, circleSamples, &vertices, &triangles);
    this->data().set_mesh(vertices.cast<double>(), triangles);

    // Set color
    Eigen::MatrixXd color(1, 3);
    color(0, 0) = yarn.color.r;
    color(0, 1) = yarn.color.g;
    color(0, 2) = yarn.color.b;
    this->data().set_colormap(color);

    // Draw center line
    if (yarn.points.rows() >= 2) {
      Eigen::MatrixXi E(yarn.points.rows() - 1, 2);
      for (int i = 0; i < yarn.points.rows() - 1; i++) {
        E(i, 0) = i;
        E(i, 1) = i + 1;
      }
      this->data().set_edges(yarn.points.cast<double>(),
          E, Eigen::RowVector3d(1, 1, 1));
    }
  }
}

void Viewer::loadYarn(std::string filename) {
  // Load .yarns file
  std::cout << "Loading model: " << filename << std::endl;
  file_format::Yarns::Yarns yarns;
  try {
    yarns = file_format::Yarns::Yarns::load(filename);
  } catch (const std::runtime_error& e) {
    std::cout << "Failed to load " << filename << std::endl;
    std::cout << e.what() << std::endl;
  }

  // Update simulator
  simulator = simulator::DiscreteSimulator(file_format::YarnRepr(yarns),
    simulator::SimulatorParams::Default());
  this->refresh();
}

} // namespace UI
