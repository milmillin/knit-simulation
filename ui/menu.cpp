#include "menu.h"

#include <iostream>

namespace UI {
  bool Menu::load(std::string filename) {
    reinterpret_cast<Viewer*>(viewer)->loadYarn(filename);
    return true;
  }

  void Menu::init(igl::opengl::glfw::Viewer *_viewer) {
    igl::opengl::glfw::imgui::ImGuiMenu::init(_viewer);

    this->callback_draw_viewer_menu = [&]() {
      // Draw parent menu content
      draw_viewer_menu();
      
      // Add simulator menu
      if (ImGui::CollapsingHeader("Simulator", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("Step", ImVec2(-1,0))) {
          std::cout << "Step button clicked" << std::endl;
        }
      }
    };
  }
} //namespace UI