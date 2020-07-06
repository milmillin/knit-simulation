#include "./viewer.h"
#include "sweep.h"

namespace UI {

int Viewer::launch(bool resizable, bool fullscreen, const std::string &name, int width, int height) {
  // Add menu
  menu.reset(new Menu());
  this->plugins.push_back(menu.get());

  // Start trackball mode (allow panning)
  this->core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

  // Launch
  return igl::opengl::glfw::Viewer::launch(resizable, fullscreen, name, width, height);
}

void Viewer::plot(Eigen::MatrixXf points, std::unordered_set<int> checkpoints) {
  Eigen::MatrixXf vertices;
  Eigen::MatrixXi triangles;

  // TODO: don't hard-code radius
  circleSweep(points, 0.1, vertices, triangles, 8);

  this->data().set_mesh(vertices.cast<double>(), triangles);

  Eigen::MatrixXi E(points.rows() - 1, 2);
  for (int i = 0; i < points.rows() - 1; i++) {
    E(i, 0) = i;
    E(i, 1) = i + 1;
  }
  this->data().set_edges(points.cast<double>(),
      E, Eigen::RowVector3d(1, 1, 1));
  
  for (int checkpoint : checkpoints) {
    this->data().add_label(points.row(checkpoint).cast<double>(), "checkpoint");
  }
}

void Viewer::loadYarn(std::string filename) {
  std::cout << "Load yarn" << std::endl;
}

} // namespace UI