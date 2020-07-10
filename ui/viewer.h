#pragma once

#include <Eigen/Core>
#include <imgui/imgui.h>

#include <unordered_set>
#include <string>
#include <memory>

#include "./menu.h"
#include "../simulator/DiscreteSimulator.h"

namespace UI {

class Menu;

class Viewer : igl::opengl::glfw::Viewer {
 public:
  Viewer() {}
  int launch(bool resizable = true, bool fullscreen = false,
    const std::string &name = "GRAIL Knit Simulator", int width = 0, int height = 0);
  void refresh();
  void loadYarn(std::string filename);
  inline void step() { simulator.step(); }

  // Number of samples for yarn cross-section (circle)
  int circleSamples = 8;
  // Number of samples for Catmull-Rom curve
  int curveSamples = 1;
  simulator::DiscreteSimulator simulator;

 private:
  std::unique_ptr<Menu> _menu;
};

}  // namespace UI
