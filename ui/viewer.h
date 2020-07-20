#pragma once

#include <Eigen/Core>
#include <imgui/imgui.h>

#include <unordered_set>
#include <string>
#include <memory>
#include <map>

#include "./menu.h"
#include "../simulator/BaseSimulator.h"

namespace UI {

class Menu;

class Viewer : igl::opengl::glfw::Viewer {
 public:
  Viewer() {}
  int launch(bool resizable = true, bool fullscreen = false,
    const std::string &name = "GRAIL Knit Simulator", int width = 0, int height = 0);
  void refresh();
  void loadYarn(std::string filename);
  inline void step() {
    _simulator->step();
  }
  simulator::BaseSimulator* simulator() const { return _simulator.get(); }
  inline simulator::SimulatorParams &getParameters() {
    return _simulator.get()->params;
  }

  // Number of samples for yarn cross-section (circle)
  int circleSamples = 8;
  // Number of samples for Catmull-Rom curve
  int curveSamples = 1;
  // Type of simulator to use
  enum SimulatorClass {
    Contenious = 0,
    Discrete = 1
  } simulatorClass = Contenious;

 private:
  std::unique_ptr<Menu> _menu;
  std::unique_ptr<simulator::BaseSimulator> _simulator;
};

}  // namespace UI
