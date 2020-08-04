#include "./viewer.h"

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include "./sweep.h"
#include "../simulator/SimulatorParams.h"
#include "../file_format/yarns.h"
#include "../file_format/yarnRepr.h"
#include "../simulator/Helper.h"
#include "../simulator/Simulator.h"
#include "../simulator/DiscreteSimulator.h"

namespace UI {

static constexpr const char* VIEWER_STATE_NAME = "viewer-state.txt";

Viewer::Viewer(std::string outputDirectory, bool reload) : _outputDirectory(outputDirectory), _reload(reload) {
  callback_pre_draw = [&](igl::opengl::glfw::Viewer&)-> bool {
    //_refreshLock.lock();
    if (_needRefresh) {
      refresh();
    }
    return false;
  };

  callback_post_draw = [&](igl::opengl::glfw::Viewer&)-> bool {
    //_refreshLock.unlock();
    return true;
  };

  //loadYarn(filename);
  if (_outputDirectory.back() != '/') {
    _outputDirectory.push_back('/');
  }
}

int Viewer::launch(bool resizable, bool fullscreen, const std::string &name, int width, int height) {
  assert(_loaded);
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

void Viewer::launchNoGUI() {
  assert(_loaded);
  this->animationManager()->startSimulation();
  this->animationManager()->join();   // wait indefinitely
}

void Viewer::refresh() {
  assert(std::this_thread::get_id() == _threadId
    && "refresh() must be called from the thread the viewer has been constructed");
  std::lock_guard<std::recursive_mutex> lock(_mutex);
  // std::lock_guard<std::recursive_mutex> guard(_refreshLock);

  // Get yarn shape
  const file_format::YarnRepr yarns = _history->getFrame(_currentFrame);
  const simulator::SimulatorParams& params = _simulator->getParams();

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
    color(0, 0) = yarn.color.r / 255.f;
    color(0, 1) = yarn.color.g / 255.f;
    color(0, 2) = yarn.color.b / 255.f;
    this->data().set_colors(color);

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

  _needRefresh = false;
}

void Viewer::createSimulator() {
  delete (_animationManager.release());
  delete (_history.release());
  delete (_simulator.release());

  // Update simulator
  switch (simulatorType)
  {
  case simulator::SimulatorType::Continuous:
    std::cout << "Using continuous simulator" << std::endl;
    _simulator.reset(new simulator::Simulator(_yarnsRepr, params));
    break;

  case simulator::SimulatorType::Discrete:
    std::cout << "Using discrete simulator" << std::endl;
    _simulator.reset(new simulator::DiscreteSimulator(_yarnsRepr, params));
    break;

  default:
    assert(false && "Invalid simulator class");
    break;
  }

  // Reset manager
  _animationManager.reset(new AnimationManager(this));
  _history.reset(new HistoryManager(this, _yarnsRepr));
  
  _currentFrame = 0;

  this->refresh();
}

void Viewer::loadYarn(const std::string& filename) {
  // Load .yarns file
  simulator::log() << "Loading model: " << filename << std::endl;
  file_format::Yarns::Yarns yarns;
  try {
    yarns = file_format::Yarns::Yarns::load(filename);
  }
  catch (const std::runtime_error& e) {
    std::cout << "Failed to load " << filename << std::endl;
    std::cout << e.what() << std::endl;
  }

  _yarnsRepr = file_format::YarnRepr(yarns);

  if (_reload) {
    // Restoring State
    simulator::log() << "Try restoring state and history from " << _outputDirectory << std::endl;
    fs::create_directory(_outputDirectory);

    file_format::ViewerState state(_outputDirectory + VIEWER_STATE_NAME);
    simulatorType = state.getType();
    params = state.getParams();
    int numSteps = state.getNumSteps();

    createSimulator();

    char positionName[200];
    char velocityName[200];
    for (int i = 2; i <= numSteps; i++) {
      snprintf(positionName, 200, "%sposition-%05d.yarns", _outputDirectory.c_str(), i);
      snprintf(velocityName, 200, "%svelocity-%05d.yarns", _outputDirectory.c_str(), i);
      if (!fs::exists(positionName)
        || !fs::exists(velocityName)) {
        std::cout << "WARNING: " << i - 1 << "frames out of " << numSteps << " restored." << std::endl;
        numSteps = i - 1;
        break;
      }

      file_format::Yarns::Yarns position = file_format::Yarns::Yarns::load(positionName);
      _history->addFrame(file_format::YarnRepr(position));
    }

    // Load last position and velocity
    if (numSteps > 1) {
      snprintf(velocityName, 200, "%svelocity-%05d.yarns", _outputDirectory.c_str(), numSteps);
      file_format::Yarns::Yarns velocity = file_format::Yarns::Yarns::load(velocityName);
      _simulator->setVelocity(file_format::YarnRepr(velocity));
      _simulator->setPosition(_history->getFrame(_history->totalFrameNumber() - 1));
    }

    simulator::log() << "> Frame 1 to " << numSteps << " restored." << std::endl;
  }
  else {
    createSimulator();
  }

  _loaded = true;
}

void Viewer::saveYarn(const std::string& filename) {
  _history->getFrame(_currentFrame).toYarns().save(filename);
}

void Viewer::nextFrame() {
  std::lock_guard<std::recursive_mutex> lock(_mutex);
  if (_currentFrame + 1 < _history->totalFrameNumber()) {
    _currentFrame++;
    invalidate();
  }
}

void Viewer::prevFrame() {
  std::lock_guard<std::recursive_mutex> lock(_mutex);
  if (_currentFrame > 0) {
    _currentFrame--;
    invalidate();
  }
}

void Viewer::setAnimationMode(bool animating) {
  std::lock_guard<std::recursive_mutex> lock(_mutex);
  core().is_animating = animating;
}


void Viewer::saveState() const {
  file_format::ViewerState state(simulatorType,
    _simulator->getParams(),
    _history->totalFrameNumber());

  state.save(_outputDirectory + VIEWER_STATE_NAME);
}

} // namespace UI
