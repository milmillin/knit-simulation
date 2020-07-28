#include "./AnimationManager.h"
#include <chrono>

namespace UI {

AnimationManager::AnimationManager(Viewer *parent) : 
    _parent(parent), _terminate(false),
    _animationRunning(false), _simulationRunning(false) {
  _simulationThread = std::thread(&AnimationManager::runSimulation, this);
  _animationThread = std::thread(&AnimationManager::runAnimation, this);
}

AnimationManager::~AnimationManager() {
  _terminate = true;
  _statusChanged.notify_all();

  _simulationThread.join();
  _animationThread.join();
}

bool AnimationManager::isSimulationRunning() const {
  std::lock_guard<std::mutex> lock(_statusLock);
  return _simulationRunning;
}

void AnimationManager::startSimulation() {
  {
    std::lock_guard<std::mutex> lock(_statusLock);
    _simulationRunning = true;
  }
  _statusChanged.notify_all();
}

void AnimationManager::stopSimulation() {
  {
    std::lock_guard<std::mutex> lock(_statusLock);
    _simulationRunning = false;
  }
  _statusChanged.notify_all();
}

bool AnimationManager::isTerminated() const {
  std::lock_guard<std::mutex> lock(_statusLock);
  return _terminate;
}

bool AnimationManager::isAnimationRunning() const {
  std::lock_guard<std::mutex> lock(_statusLock);
  return _animationRunning;
}

void AnimationManager::startAnimation() {
  {
    std::lock_guard<std::mutex> lock(_statusLock);
    _animationRunning = true;
    _parent->setAnimationMode(true);
  }
  _statusChanged.notify_all();
}

void AnimationManager::stopAnimation() {
  {
    std::lock_guard<std::mutex> lock(_statusLock);
    _animationRunning = false;
    _parent->setAnimationMode(false);
  }
  _statusChanged.notify_all();
}

void AnimationManager::runSimulation() {
  std::function<bool()> terminated = [this]()->bool { return isTerminated(); };
  while (true) {
    std::unique_lock<std::mutex> lock(_statusLock);
    _statusChanged.wait(lock, [this]{return _terminate || _simulationRunning;});
    if (_terminate) {
      return;
    }
    lock.unlock();

    _parent->simulator()->step(terminated);
    _parent->history()->addFrame(_parent->simulator()->getYarns());
    _parent->saveState();
  }
}

void AnimationManager::runAnimation() {
  while (true) {
    std::unique_lock<std::mutex> lock(_statusLock);
    _statusChanged.wait(lock, [this]{return _terminate || _animationRunning;});
    if (_terminate) {
      return;
    }
    lock.unlock();

    std::this_thread::sleep_for(
      std::chrono::milliseconds(_parent->animationPlaybackInterval));
    _parent->nextFrame();
  }
}

}  // namespace UI
