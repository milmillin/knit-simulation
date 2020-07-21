#include "./HistoryManager.h"

namespace UI {

HistoryManager::HistoryManager(Viewer *parent, const file_format::YarnRepr &firstFrame)
    : _parent(parent) {
  _history.push_back(firstFrame);
}

bool HistoryManager::hasNext() const {
  return _currentFrame < _history.size() - 1;
}

bool HistoryManager::hasPrev() const {
  return _currentFrame > 0;
}

void HistoryManager::next() {
  std::lock_guard<std::mutex> guard(_lock);
  if (hasNext()) {
    _currentFrame++;
  }
}

void HistoryManager::prev() {
  std::lock_guard<std::mutex> guard(_lock);
  if (hasPrev()) {
    _currentFrame--;
  }
}

void HistoryManager::goToStart() {
  std::lock_guard<std::mutex> guard(_lock);
  _currentFrame = 0;
}

void HistoryManager::goToEnd() {
  std::lock_guard<std::mutex> guard(_lock);
  _currentFrame = _history.size() - 1;
}

int HistoryManager::currentFrameNumber() const {
  return _currentFrame;
}

int HistoryManager::totalFrameNumber() const {
  return _history.size();
}

void HistoryManager::addFrame(const file_format::YarnRepr &yarn) {
  std::lock_guard<std::mutex> guard(_lock);
  _history.push_back(yarn);
  _currentFrame++;
}

const file_format::YarnRepr HistoryManager::curentFrame() {
  std::lock_guard<std::mutex> guard(_lock);
  return _history.at(_currentFrame);
}

}  // namespace UI