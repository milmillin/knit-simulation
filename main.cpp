#include "ui/viewer.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " yarns-filename" << std::endl;
    return -1;
  }

  // Launch viewer
  UI::Viewer viewer("output/");
  viewer.loadYarn(argv[1]);
  viewer.launch();
}
