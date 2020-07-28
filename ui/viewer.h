#pragma once

#include <Eigen/Core>
#include <imgui/imgui.h>

#include <unordered_set>
#include <string>
#include <memory>
#include <map>
#include <mutex>

#include "./menu.h"
#include "./HistoryManager.h"
#include "./AnimationManager.h"
#include "../simulator/BaseSimulator.h"
#include "../simulator/Enum.h"
#include "../file_format/ViewerState.h"

namespace UI {

class Menu;
class HistoryManager;
class AnimationManager;

class Viewer : igl::opengl::glfw::Viewer {
 public:
  Viewer(std::string filename);
  int launch(bool resizable = true, bool fullscreen = false,
    const std::string &name = "GRAIL Knit Simulator", int width = 0, int height = 0);
  void refresh();
  void loadYarn(const std::string& filename);
  void saveYarn(const std::string& filename);
  void createSimulator();

  void setAnimationMode(bool animating);
  void nextFrame();
  void prevFrame();

  void saveState() const;

  simulator::BaseSimulator* simulator() const { return _simulator.get(); }

  // Number of samples for yarn cross-section (circle)
  int circleSamples = 8;
  // Number of samples for Catmull-Rom curve
  int curveSamples = 4;
  // Animation playback interval (millisecond)
  int animationPlaybackInterval = 1000/10;
  // Type of simulator to use
  simulator::SimulatorType simulatorType = simulator::SimulatorType::Discrete;

  simulator::SimulatorParams params;

  HistoryManager* history() {return _history.get();}
  AnimationManager* animationManager() {return _animationManager.get();}
  int& currentFrame() { return _currentFrame; }

 private:
  std::unique_ptr<Menu> _menu;
  std::unique_ptr<simulator::BaseSimulator> _simulator;
  std::unique_ptr<HistoryManager> _history;
  std::unique_ptr<AnimationManager> _animationManager;
  mutable std::recursive_mutex _refreshLock;
  int _currentFrame = 0;
  file_format::YarnRepr _yarnsRepr;
};

}  // namespace UI
