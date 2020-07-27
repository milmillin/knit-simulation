#pragma once

#include <vector>
#include <mutex>

#include "./viewer.h"
#include "../file_format/yarnRepr.h"

namespace UI {

class Viewer;

class HistoryManager {
 public:
  HistoryManager(Viewer *parent, const file_format::YarnRepr &firstFrame);
  int totalFrameNumber() const;
  void addFrame(const file_format::YarnRepr &yarn);
  file_format::YarnRepr getFrame(int index) const;

 private:
  Viewer *_parent;
  std::vector<file_format::YarnRepr> _history;
  int _currentFrame = 0;
  mutable std::mutex _lock;
};

}  // namespace UI
