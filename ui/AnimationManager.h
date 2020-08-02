#pragma once

#include <mutex>
#include <thread>
#include <condition_variable>

#include "./viewer.h"
#include "./HistoryManager.h"

namespace UI {

class Viewer;
class HistoryManager;

class AnimationManager {
 public:
  AnimationManager(Viewer *parent);
  ~AnimationManager();

  bool isSimulationRunning() const;
  void startSimulation();
  void stopSimulation();

  bool isAnimationRunning() const;
  void startAnimation();
  void stopAnimation();

  bool isTerminated() const;

  void join();

 private:
  Viewer *_parent;

  bool _terminate;

  bool _simulationRunning;
  std::thread _simulationThread;

  bool _animationRunning;
  std::thread _animationThread;

  mutable std::mutex _statusLock;
  std::condition_variable _statusChanged;

  void runSimulation();
  void runAnimation();
};

}  // namspace UI