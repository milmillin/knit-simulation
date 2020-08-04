#include "./easy_profiler_stub.h"

#include "ui/viewer.h"


int main(int argc, char *argv[]) {

  std::vector<std::string> args(argc);
  for (int i = 0; i < argc; i++) {
    args[i] = argv[i];
  }

  if (argc < 2) {
    std::cout << "Usage: " << args[0] << " yarns-filename [--no-restore] [--no-gui] [--output output-directory]" << std::endl;
    return -1;
  }

  // Restore /output directory
  bool restore = true;
  // Launch GUI
  bool gui = true;
  // Output directory
  std::string outputDir = "output/";

  for (int i = 2; i < argc; i++) {
    const std::string& arg = args[i];
    if (arg == "--no-restore") {
      restore = false;
    }
    else if (arg == "--no-gui") {
      gui = false;
    }
    else if (arg == "--output") {
      if (i + 1 < argc) {
        outputDir = args[i + 1];
        if (outputDir.back() != '/' && outputDir.back() != '\\') {
          outputDir.push_back('/');
        }
        i++;
      }
      else {
        std::cout << "ERROR: cannot read --output option parameter" << std::endl;
        return -1;
      }
    }
    else {
      std::cout << "WARNING: " << arg << " does not exist" << std::endl;
    }
  }


  EASY_PROFILER_ENABLE;

  // Launch viewer
  UI::Viewer viewer(outputDir, restore);
  viewer.loadYarn(args[1]);
  if (gui) viewer.launch();
  else viewer.launchNoGUI();

  #ifdef USE_EASY_PROFILER
    profiler::dumpBlocksToFile("test_profile.prof");
  #endif
}
