#include "./easy_profiler_stub.h"

#include "ui/viewer.h"


int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " yarns-filename" << std::endl;
    return -1;
  }

  EASY_PROFILER_ENABLE;

  // Launch viewer
  UI::Viewer viewer(argv[1]);
  viewer.launch();

  #ifdef ENABLE_EASY_PROFILER
    profiler::dumpBlocksToFile("test_profile.prof");
  #endif
}
