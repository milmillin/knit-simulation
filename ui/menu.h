#pragma once

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include "viewer.h"

namespace UI {
  class Menu : public igl::opengl::glfw::imgui::ImGuiMenu {
   public:
    virtual void init(igl::opengl::glfw::Viewer *_viewer) override;

   private:
    virtual bool load(std::string filename) override;
  };
} // namespace UI