#pragma once

#include <vector>

#include "./viewer.h"
#include "../file_format/yarnRepr.h"

namespace UI {

class Viewer;

class HistoryManager {
 public:
  HistoryManager(Viewer *parent) : _parent(parent) {}
  bool hasNext() const;
  bool hasPrev() const;
  void next();
  void prev();
  void end();
  int currentFrameNumber() const;
  int totalFrameNumber() const;
  void addFrame(const file_format::YarnRepr &yarn);
  const file_format::YarnRepr& curentFrame();

 private:
  Viewer *_parent;
  std::vector<file_format::YarnRepr> _history;
  int _currentFrame = 0;
};

}  // namespace UI
