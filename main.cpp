#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include "file_format/yarns.h"
#include "simulator/Simulator.h"
#include "simulator/SimulatorParams.h"

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  std::cout << "Hello, world!" << std::endl;

  // Load .yarns file
  file_format::Yarns yarn;
  try {
    yarn = file_format::Yarns::load("../../../helloworld.yarns");
  } catch (const std::runtime_error& e) {
    std::cout << e.what() << std::endl;
    return 0;
  }

  // Construct polyline
  int n = yarn.yarns[0].points.size();
  Eigen::MatrixXf P(n, 3);
  for (int i = 0; i < n; i++) {
    P(i, 0) = yarn.yarns[0].points[i][0];
    P(i, 1) = yarn.yarns[0].points[i][1];
    P(i, 2) = yarn.yarns[0].points[i][2];
  }

  Eigen::MatrixXi E(n - 1, 2);
  for (int i = 0; i < n - 1; i++) {
    E(i, 0) = i;
    E(i, 1) = i + 1;
  }

  // Feed to the simulator
  simulator::Simulator simulator(P, simulator::SimulatorParams::Default());
  simulator.step();

  // TODO: Construct tube mesh from polyline

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  double doubleVariable = 0.1;

  menu.callback_draw_viewer_menu = [&]() {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Expose variable directly ...
      ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

      // ... or using a custom callback
      static bool boolVariable = true;
      if (ImGui::Checkbox("bool", &boolVariable))
      {
        // do something
        std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
      }

      // Expose an enumeration type
      enum Orientation { Up=0, Down, Left, Right };
      static Orientation dir = Up;
      ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

      // We can also use a std::vector<std::string> defined dynamically
      static int num_choices = 3;
      static std::vector<std::string> choices;
      static int idx_choice = 0;
      if (ImGui::InputInt("Num letters", &num_choices))
      {
        num_choices = std::max(1, std::min(26, num_choices));
      }
      if (num_choices != (int) choices.size())
      {
        choices.resize(num_choices);
        for (int i = 0; i < num_choices; ++i)
          choices[i] = std::string(1, 'A' + i);
        if (idx_choice >= num_choices)
          idx_choice = num_choices - 1;
      }
      ImGui::Combo("Letter", &idx_choice, choices);

      // Add a button
      if (ImGui::Button("Print Hello", ImVec2(-1,0)))
      {
        std::cout << "Hello\n";
      }
    }
  };

  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.data().set_edges(simulator.getControlPoints().cast<double>(),
      E, Eigen::RowVector3d(1, 1, 1));
  viewer.launch();
}
