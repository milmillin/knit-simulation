#include "./menu.h"

#include <spdlog/spdlog.h>

#include <iostream>
#include <string>

constexpr int inputBoxWidth = 200;

namespace UI {

bool Menu::load(std::string filename) {
  reinterpret_cast<Viewer*>(viewer)->loadYarn(filename);
  return true;
}

bool Menu::save(std::string filename) {
  reinterpret_cast<Viewer*>(viewer)->saveYarn(filename);
  return true;
}

void Menu::init(igl::opengl::glfw::Viewer* _viewer) {
  igl::opengl::glfw::imgui::ImGuiMenu::init(_viewer);

  this->callback_draw_viewer_menu = [&]() {
    // Based on 'ImGuiMenu::draw_viewer_menu'
    bool needRefresh = false;
    Viewer* yarnViewer = reinterpret_cast<Viewer*>(viewer);
    simulator::SimulatorParams& params = yarnViewer->params;

    // Mesh
    if (ImGui::CollapsingHeader("Yarn", ImGuiTreeNodeFlags_DefaultOpen)) {
      float w = ImGui::GetContentRegionAvailWidth();
      float p = ImGui::GetStyle().FramePadding.x;
      if (ImGui::Button("Load##Mesh", ImVec2((w - p) / 2.f, 0))) {
        viewer->open_dialog_load_mesh();
      }
      ImGui::SameLine(0, p);

      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
      if (ImGui::Button("Save##Mesh", ImVec2((w - p) / 2.f, 0))) {
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
          }
          else if (viewer->core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION) {
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
    auto make_checkbox = [&](const char* label, unsigned int& option) {
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
      ImGui::Checkbox("Show control point labels", &(viewer->data().show_labels));
      if (ImGui::Checkbox("Show material frames", &(yarnViewer->showMaterialFrames)))
        needRefresh = true;
      if (ImGui::Checkbox("Show bishop frames", &(yarnViewer->showBishopFrame)))
        needRefresh = true;
      ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
      ImGui::DragFloat("Shininess", &(viewer->data().shininess), 0.05f, 0.0f, 100.0f);
      if (ImGui::DragInt("Cross-section samples",
        &(yarnViewer->circleSamples),
        1, 4, 10)) {
        needRefresh = true;
      }
      if (ImGui::DragInt("Curve samples",
        &(yarnViewer->curveSamples),
        1, 1, 10)) {
        needRefresh = true;
      }
      ImGui::PopItemWidth();
    }

    if (ImGui::CollapsingHeader("Colors")) {
      ImGui::ColorEdit4("Background color", viewer->core().background_color.data(),
        ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
      if (ImGui::ColorEdit3("Material Frame U", yarnViewer->materialFrameUColor.data()))
        needRefresh = true;
      if (ImGui::ColorEdit3("Material Frame V", yarnViewer->materialFrameVColor.data()))
        needRefresh = true;
      if (ImGui::ColorEdit3("Bishop Frame U", yarnViewer->bishopFrameUColor.data()))
        needRefresh = true;
      if (ImGui::ColorEdit3("Bishop Frame V", yarnViewer->bishopFrameVColor.data()))
        needRefresh = true;
    }

    // Overlays
    if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
      make_checkbox("Wireframe", viewer->data().show_lines);
      make_checkbox("Fill", viewer->data().show_faces);
    }

    // === Simulator menu ===
    ImGui::Begin("Simulator", NULL, ImGuiWindowFlags_AlwaysAutoResize);
    if (ImGui::CollapsingHeader("Simulation and Playback", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::PushID("Animation");
      ImGui::PushItemWidth(inputBoxWidth);
      if (yarnViewer->animationManager()->isSimulationRunning()) {
        if (ImGui::Button("Stop Simulation")) {
          yarnViewer->animationManager()->stopSimulation();
        }
      }
      else {
        if (ImGui::Button("Start Simulation")) {
          yarnViewer->animationManager()->startSimulation();
        }
      }
      ImGui::PushID("Player");
      ImGui::Text("Frame %d of %d",
        yarnViewer->currentFrame() + 1,
        yarnViewer->history()->totalFrameNumber());
      if (ImGui::Button("First")) {
        yarnViewer->currentFrame() = 0;
        needRefresh = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Prev")) {
        yarnViewer->prevFrame();
      }
      ImGui::SameLine();
      if (yarnViewer->animationManager()->isAnimationRunning()) {
        if (ImGui::Button("Pause")) {
          yarnViewer->animationManager()->stopAnimation();
        }
      }
      else {
        if (ImGui::Button("Play")) {
          yarnViewer->animationManager()->startAnimation();
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Next")) {
        yarnViewer->nextFrame();
      }
      ImGui::SameLine();
      if (ImGui::Button("Last")) {
        yarnViewer->currentFrame() = yarnViewer->history()->totalFrameNumber() - 1;
        needRefresh = true;
      }
      // TODO: concurrency issue here
      ImGui::InputInt("replay frame interval(ms)", &(yarnViewer->animationPlaybackInterval),
        10, 100);
      ImGui::PopItemWidth();
      ImGui::PopID();;
      ImGui::PopID();;
    }

    if (ImGui::CollapsingHeader("Create Simulator", 0)) {
      ImGui::PushItemWidth(inputBoxWidth);
      if (ImGui::Combo("Simulator class",
        reinterpret_cast<int*>(&yarnViewer->simulatorType),
        "Simulator\0"
        "Discrete Simulator\0\0")) {
        SPDLOG_INFO("Simulator Changed");
      }
      ImGui::Checkbox("Debug mode", &(params.debug));
      ImGui::Checkbox("Print statistics", &(params.statistics));
      ImGui::Checkbox("Enable ground", &(params.enableGround));
      ImGui::InputDouble("Time resolution", &(params.h));
      ImGui::InputInt("steps per frame", &(params.steps),
        10, 100);

      ImGui::Separator();

      ImGui::Text("Constraint and Fast Projection");
      ImGui::InputDouble("Scale fator", &(params.cInit));
      ImGui::Checkbox("Enable length constrain", &(params.enableLengthConstrain));
      ImGui::InputDouble("Target error", &(params.fastProjErrorCutoff));
      ImGui::InputInt("Max iterations", &(params.fastProjMaxIter), 1, 5);

      ImGui::Separator();

      ImGui::Text("Parameters");
      ImGui::InputDouble("Gravity", &(params.gravity));
      ImGui::InputDouble("Ground height", &(params.groundHeight));
      ImGui::InputDouble("Ground fiction", &(params.groundFriction));
      ImGui::InputInt("Contact force samples", &(params.contactForceSamples),
        1, 5);
      ImGui::InputDouble("Contact force model tolerance", &(params.contactModelTolerance));
      ImGui::InputInt("Max model update interval", &(params.maxContactModelUpdateInterval));
      ImGui::InputDouble("Contact force", &(params.kContact));
      ImGui::InputDouble("Global damping", &(params.kGlobalDamping));
      ImGui::InputDouble("Length Constant", &(params.kLen));
      ImGui::InputDouble("Bending Constant", &(params.kBend));
      ImGui::InputDouble("Material frame tolerance", &(params.materialFrameTolerance));
      ImGui::InputInt("Material frame max iterations", &(params.materialFrameMaxUpdate), 1, 5);
      ImGui::InputDouble("Twisting Constant", &(params.kTwist));
      ImGui::PopItemWidth();
      if (ImGui::Button("Create")) {
        yarnViewer->createSimulator();
      }
    }
    ImGui::End();

    if (needRefresh) {
      yarnViewer->invalidate();
    }
  };
}

}  // namespace UI
