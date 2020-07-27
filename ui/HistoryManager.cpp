#include "./HistoryManager.h"

namespace UI {

HistoryManager::HistoryManager(Viewer *parent, const file_format::YarnRepr &firstFrame)
    : _parent(parent) {
  _history.push_back(firstFrame);
}

int HistoryManager::totalFrameNumber() const {
  return _history.size();
}

void HistoryManager::addFrame(const file_format::YarnRepr &yarn) {
  std::lock_guard<std::mutex> guard(_lock);
  _history.push_back(yarn);
}

file_format::YarnRepr HistoryManager::getFrame(int index) const {
  std::lock_guard<std::mutex> guard(_lock);
  return _history[index];
}

}  // namespace UI