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

void Viewer::clearCache() {
  _cache.clear();
}

void Viewer::refresh() {
  auto it = _cache.find(currentStep);

  data().clear();

  if (it != _cache.end()) {
    data().set_mesh(it->second.vertices, it->second.triangles);
  } else {
    // Get yarn shape
    const file_format::YarnRepr& yarns = _simulator.getYarns(currentStep);

    for (int i = 0; i < yarns.yarns.size(); i++) {
      // Clear old mesh
      this->selected_data_index = i;
      this->data().clear();

      // Get curve
      auto& yarn = yarns.yarns[i];

      // Use Catmull-Rom to smooth the curve
      Eigen::MatrixXf points = simulator::catmullRomSequenceSample(yarn.points, curveSamples);

      // Create mesh for tube
      Eigen::MatrixXf vertices;
      Eigen::MatrixXi triangles;
      circleSweep(points, yarn.radius, circleSamples, &vertices, &triangles);

      _cache[currentStep] = Geometry{ vertices.cast<double>(), triangles };

      this->data().set_mesh(vertices.cast<double>(), triangles);


    }
  }

  const file_format::Yarn& yarn = _simulator.getYarns(currentStep).yarns[0];

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

  // Set color
  Eigen::MatrixXd color(1, 3);
  color(0, 0) = yarn.color.r;
  color(0, 1) = yarn.color.g;
  color(0, 2) = yarn.color.b;
  this->data().set_colormap(color);
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
  _simulator = simulator::Simulator(file_format::YarnRepr(yarns),
    simulator::SimulatorParams::Default());
  clearCache();
  this->refresh();
}

bool Viewer::viewNext() {
  if (currentStep + 1 < _simulator.numStep()) {
    currentStep++;
    return true;
  }
  return false;
}

bool Viewer::viewPrev() {
  if (currentStep > 0) {
    currentStep--;
    return true;
  }
  return false;
}

} // namespace UI
