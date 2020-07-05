#pragma once

#include <unordered_set>

#include <Eigen/Core>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

namespace UI {

class Viewer {
 public:
  void launch();
  void plot(Eigen::MatrixXf points, std::unordered_set<int> checkpoints);

 private:
  igl::opengl::glfw::Viewer viewer;
};

} // namespace UI