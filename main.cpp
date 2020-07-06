#include "ui/viewer.h"

int main(int argc, char *argv[]) {
  // Launch viewer
  UI::Viewer viewer;
  viewer.loadYarn("../helloworld.yarns");
  viewer.launch();
}
