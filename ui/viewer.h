#pragma once

#include <unordered_set>
#include <string>
#include <memory>
#include <map>
#include <mutex>
#include <thread>

#include <Eigen/Core>
#include <imgui/imgui.h>

#include "simulator/BaseSimulator.h"
#include "simulator/Enum.h"
#include "file_format/ViewerState.h"
#include "./menu.h"
#include "./HistoryManager.h"
#include "./AnimationManager.h"

namespace UI {

class Menu;
class HistoryManager;
class AnimationManager;

enum ViewerLayerID {
  GROUND = 0,
  MATERIAL_FRAMES,
  YARNS
};

class Viewer : igl::opengl::glfw::Viewer {
 public:
  Viewer(std::string outputDirectory, bool reload = true);
  int launch(bool resizable = true, bool fullscreen = false,
    const std::string &name = "GRAIL Knit Simulator", int width = 0, int height = 0);
  void launchNoGUI();
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
  const std::string& outputDirectory() const { return _outputDirectory; }
  void invalidate() {
    std::lock_guard<std::recursive_mutex> lock(_mutex);
    _needRefresh = true;
  }

  bool showMaterialFrames = false;
  bool showBishopFrame = false;
  Eigen::Vector3f materialFrameUColor = Eigen::Vector3f(211, 67, 67) / 255;
  Eigen::Vector3f materialFrameVColor = Eigen::Vector3f(243, 207, 68) / 255;
  Eigen::Vector3f bishopFrameUColor = Eigen::Vector3f(167, 104, 246) / 255;
  Eigen::Vector3f bishopFrameVColor = Eigen::Vector3f(99, 240, 153) / 255;
 private:
  std::unique_ptr<Menu> _menu;
  std::unique_ptr<simulator::BaseSimulator> _simulator;
  std::unique_ptr<HistoryManager> _history;
  std::unique_ptr<AnimationManager> _animationManager;
  int _currentFrame = 0;
  file_format::YarnRepr _yarnsRepr;
  std::string _outputDirectory;
  bool _loaded = false;
  bool _reload;
  bool _needRefresh = false;
  std::thread::id _threadId = std::this_thread::get_id();

  mutable std::recursive_mutex _mutex;

  void refresh();

  // Generate lines that visualize the material and bishop frames.
  // yarnRepr: the yarn that need to be visualized
  // V: vertices
  // E: lines
  // C: colors
  void visualizeMaterialAndBishopFrames(const file_format::YarnRepr &yarnRepr,
    Eigen::MatrixX3d *V, Eigen::MatrixX2i *E, Eigen::MatrixX3f *C);
};

}  // namespace UI
