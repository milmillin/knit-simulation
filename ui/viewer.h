#pragma once

#include <unordered_set>
#include <string>
#include <memory>

#include <Eigen/Core>
#include <imgui/imgui.h>

#include "menu.h"
#include "../simulator/Simulator.h"

namespace UI {

class Menu;

class Viewer : igl::opengl::glfw::Viewer {
 public:
  Viewer() {};
  int launch(bool resizable = true, bool fullscreen = false, const std::string &name = "libigl viewer", int width = 0, int height = 0);
  void refresh();
  void loadYarn(std::string filename);
  void step() { _simulator.step(); }

 private:
  std::unique_ptr<Menu> _menu;
  simulator::Simulator _simulator;
};

} // namespace UI