#include "./viewer.h"

#include "./sweep.h"
#include "../simulator/SimulatorParams.h"
#include "../file_format/yarns.h"
#include "../file_format/yarnRepr.h"
#include "../simulator/Helper.h"
#include "../simulator/Simulator.h"
#include "../simulator/DiscreteSimulator.h"

namespace UI {

Viewer::Viewer()
    : _simulator(new simulator::Simulator()),
      _history(new HistoryManager(this, file_format::YarnRepr())),
      _animationManager(new AnimationManager(this)) {
}

int Viewer::launch(bool resizable, bool fullscreen, const std::string &name, int width, int height) {
  // Add menu
  _menu.reset(new Menu());
  this->plugins.push_back(_menu.get());

  // Start trackball mode (allow panning)
  this->core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);

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
  std::lock_guard<std::mutex> guard(_refreshLock);

  // Get yarn shape
  const file_format::YarnRepr &yarns = _history.get()->curentFrame();
  const simulator::SimulatorParams &params = _simulator.get()->params;

  // Draw ground
  this->selected_data_index = 0;
  Eigen::MatrixXf groundPoints(4, 3);
  groundPoints << 10, params.groundHeight, 10,
    10, params.groundHeight, -10,
    -10, params.groundHeight, -10,
    -10, params.groundHeight, 10;
  Eigen::MatrixXi groundTrianges(2, 3);
  groundTrianges << 2, 0, 1,
    3, 0, 2;
  this->data().clear();
  this->data().set_mesh(groundPoints.cast<double>(), groundTrianges);

  // Create new mesh
  while (this->data_list.size() <= yarns.yarns.size()) {
    this->data_list.push_back(igl::opengl::ViewerData());
  }

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

    // Set vertex labels
    std::vector<std::string> labels;
    for (int i = 0; i < yarn.points.rows(); i++) {
      labels.push_back(std::to_string(i));
    }
    this->data().set_labels(yarn.points.cast<double>(), labels);
    this->data().label_color = Eigen::Vector4f(1, 1, 1, 1);
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

  _animationManager.reset(new AnimationManager(this));

  file_format::YarnRepr yarnsRepr(yarns);
  simulator::SimulatorParams params;
  // Update simulator
  switch (simulatorClass)
  {
  case SimulatorClass::Contenious:
    std::cout << "Using contenious simulator" << std::endl;
    _simulator.reset(new simulator::Simulator(yarnsRepr, params));
    break;

  case SimulatorClass::Discrete:
    std::cout << "Using discrete simulator" << std::endl;
    _simulator.reset(new simulator::DiscreteSimulator(yarnsRepr, params));
    break;

  default:
    assert(false && "Invalid simulator class");
    break;
  }

  // Reset history
  _history.reset(new HistoryManager(this, yarnsRepr));

  this->refresh();
}

void Viewer::nextFrame() {
  HistoryManager &history = *_history.get(); 
  history.next();
  refresh();
}

void Viewer::prevFrame() {
  HistoryManager &history = *_history.get(); 
  history.prev();
  refresh();
}

void Viewer::setAnimationMode(bool animating) {
  core().is_animating = animating;
}

} // namespace UI
