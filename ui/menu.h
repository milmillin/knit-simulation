// Copyright 2020 Tom Lou

#pragma once

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <string>

#include "./viewer.h"

namespace UI {

class Menu : public igl::opengl::glfw::imgui::ImGuiMenu {
 public:
  void init(igl::opengl::glfw::Viewer *_viewer) override;

 private:
  bool load(std::string filename) override;
  bool save(std::string filename) override;
};

}  // namespace UI
