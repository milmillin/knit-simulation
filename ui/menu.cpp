// Copyright 2020 Tom Lou

#include "./menu.h"

#include <iostream>
#include <string>

namespace UI {
  bool Menu::load(std::string filename) {
    reinterpret_cast<Viewer*>(viewer)->loadYarn(filename);
    return true;
  }

  void Menu::init(igl::opengl::glfw::Viewer *_viewer) {
    igl::opengl::glfw::imgui::ImGuiMenu::init(_viewer);

    this->callback_draw_viewer_menu = [&]() {
      // Based on 'ImGuiMenu::draw_viewer_menu'

      Viewer* myviewer = reinterpret_cast<Viewer*>(viewer);
      bool needRefresh = false;

      // Mesh
      if (ImGui::CollapsingHeader("Yarn", ImGuiTreeNodeFlags_DefaultOpen)) {
        float w = ImGui::GetContentRegionAvailWidth();
        float p = ImGui::GetStyle().FramePadding.x;
        if (ImGui::Button("Load##Mesh", ImVec2((w-p)/2.f, 0))) {
          viewer->open_dialog_load_mesh();
        }
        ImGui::SameLine(0, p);

        // TODO: The save button doesn't work
        ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        if (ImGui::Button("Save##Mesh", ImVec2((w-p)/2.f, 0))) {
          viewer->open_dialog_save_mesh();
        }
        ImGui::PopStyleVar();
      }

      // Viewing options
      if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("Center object", ImVec2(-1, 0))) {
          viewer->core().align_camera_center(viewer->data().V, viewer->data().F);
        }

        // Zoom
        ImGui::PushItemWidth(80 * menu_scaling());
        ImGui::DragFloat("Zoom", &(viewer->core().camera_zoom), 0.05f, 0.1f, 20.0f);

        // Select rotation type
        int rotation_type = static_cast<int>(viewer->core().rotation_type);
        static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
        static bool orthographic = true;
        if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0")) {
          using RT = igl::opengl::ViewerCore::RotationType;
          auto new_type = static_cast<RT>(rotation_type);
          if (new_type != viewer->core().rotation_type) {
            if (new_type == RT::ROTATION_TYPE_NO_ROTATION) {
              trackball_angle = viewer->core().trackball_angle;
              orthographic = viewer->core().orthographic;
              viewer->core().trackball_angle = Eigen::Quaternionf::Identity();
              viewer->core().orthographic = true;
            } else if (viewer->core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION) {
              viewer->core().trackball_angle = trackball_angle;
              viewer->core().orthographic = orthographic;
            }
            viewer->core().set_rotation_type(new_type);
          }
        }

        // Orthographic view
        ImGui::Checkbox("Orthographic view", &(viewer->core().orthographic));
        ImGui::PopItemWidth();
      }

      // Helper for setting viewport specific mesh options
      auto make_checkbox = [&](const char *label, unsigned int &option) {
        return ImGui::Checkbox(label,
          [&]() { return viewer->core().is_set(option); },
          [&](bool value) { return viewer->core().set(option, value); });
      };

      // Draw options
      if (ImGui::CollapsingHeader("Rendering Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Checkbox("Face-based shading", &(viewer->data().face_based))) {
          viewer->data().dirty = igl::opengl::MeshGL::DIRTY_ALL;
        }
        make_checkbox("Show texture", viewer->data().show_texture);
        if (ImGui::Checkbox("Invert normals", &(viewer->data().invert_normals))) {
          viewer->data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
        }
        make_checkbox("Show overlay", viewer->data().show_overlay);
        make_checkbox("Show overlay depth", viewer->data().show_overlay_depth);
        ImGui::ColorEdit4("Background color", viewer->core().background_color.data(),
            ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
        ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
        ImGui::DragFloat("Shininess", &(viewer->data().shininess), 0.05f, 0.0f, 100.0f);
        if (ImGui::DragInt("Cross-section samples",
            &(reinterpret_cast<Viewer*>(viewer)->circleSamples),
            1, 4, 10)) {
          needRefresh = true;
        }
        if (ImGui::DragInt("Curve samples",
            &(reinterpret_cast<Viewer*>(viewer)->curveSamples),
            1, 1, 10)) {
          needRefresh = true;
        }
        ImGui::PopItemWidth();
      }

      // Overlays
      if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
        make_checkbox("Wireframe", viewer->data().show_lines);
        make_checkbox("Fill", viewer->data().show_faces);
      }


      // Add simulator menu
      if (ImGui::CollapsingHeader("Simulator", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("Step", ImVec2(-1, 0))) {
          myviewer->step();
          myviewer->viewNext();
          needRefresh = true;
        }
        ImGui::Separator();
        ImGui::Text("Viewing %d of %d\n", myviewer->currentStep + 1, myviewer->numStep());

        if (ImGui::Button("View Prev", ImVec2(-1, 0))) {
          needRefresh = needRefresh || myviewer->viewPrev();
        }
        if (ImGui::Button("View Next", ImVec2(-1, 0))) {
          needRefresh = needRefresh || myviewer->viewNext();
        }
        if (ImGui::Button("View First", ImVec2(-1, 0))) {
          myviewer->currentStep = 0;
          needRefresh = true;
        }
      }

      // Refresh the mesh when config changes
      if (needRefresh) {
        reinterpret_cast<Viewer*>(viewer)->refresh();
      }
    };
  }
}  // namespace UI
