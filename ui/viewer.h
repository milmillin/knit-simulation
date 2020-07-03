#pragma once

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

namespace UI {

class Viewer {
 public:
  void launch();
  void plot(Eigen::MatrixXf points);

 private:
  igl::opengl::glfw::Viewer viewer;
};

} // namespace UI