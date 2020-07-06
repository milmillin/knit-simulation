#pragma once

#include <unordered_set>
#include <string>
#include <memory>

#include <Eigen/Core>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include "menu.h"

namespace UI {

class Menu;

class Viewer : igl::opengl::glfw::Viewer {
 public:
  int launch(bool resizable = true, bool fullscreen = false, const std::string &name = "libigl viewer", int width = 0, int height = 0);
  void plot(Eigen::MatrixXf points, std::unordered_set<int> checkpoints);
  void loadYarn(std::string filename);

 private:
  std::unique_ptr<Menu> menu;
};

} // namespace UI