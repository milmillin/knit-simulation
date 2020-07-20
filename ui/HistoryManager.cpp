#include "./HistoryManager.h"

namespace UI {

bool HistoryManager::hasNext() const {
  return _currentFrame < _history.size() - 1;
}

bool HistoryManager::hasPrev() const {
  return _currentFrame > 0;
}

void HistoryManager::next() {
  if (hasNext()) {
    _currentFrame++;
  }
}

void HistoryManager::prev() {
  if (hasPrev()) {
    _currentFrame--;
  }
}
void HistoryManager::end() {
  _currentFrame = _history.size() - 1;
}

int HistoryManager::currentFrameNumber() const {
  return _currentFrame;
}

int HistoryManager::totalFrameNumber() const {
  return _history.size();
}

void HistoryManager::addFrame(const file_format::YarnRepr &yarn) {
  _history.push_back(yarn);
  next();
}

const file_format::YarnRepr& HistoryManager::curentFrame() {
  return _history.at(_currentFrame);
}

}  // namespace UI