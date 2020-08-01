#include "ViewerState.h"

#include <iomanip>
#include <fstream>
#include <cassert>

namespace file_format {

#define ADD_ENTRY(label, value) \
  _entries.push_back(std::unique_ptr<IEntry>(new Entry<decltype(value)>(label, &value)))


ViewerState::ViewerState(simulator::SimulatorType type, simulator::SimulatorParams params, int numSteps)
  : _type(type), _params(params), _numSteps(numSteps) {
  ADD_ENTRY("numSteps", _numSteps);
  ADD_ENTRY("simulatorType", _type);
  ADD_ENTRY("debug", _params.debug);
  ADD_ENTRY("m", _params.m);
  ADD_ENTRY("kLen", _params.kLen);
  ADD_ENTRY("kBend", _params.kBend);
  ADD_ENTRY("kGlobal", _params.kGlobal);
  ADD_ENTRY("kContact", _params.kContact);
  ADD_ENTRY("kDt", _params.kDt);
  ADD_ENTRY("kDn", _params.kDn);
  ADD_ENTRY("aSmall", _params.aSmall);
  ADD_ENTRY("aLarge", _params.aLarge);
  ADD_ENTRY("h", _params.h);
  ADD_ENTRY("steps", _params.steps);
  ADD_ENTRY("gravity", _params.gravity);
  ADD_ENTRY("groundHeight", _params.groundHeight);
  ADD_ENTRY("groundFriction", _params.groundFriction);
  ADD_ENTRY("cInit", _params.cInit);
  ADD_ENTRY("fastProjMaxIter", _params.fastProjMaxIter);
  ADD_ENTRY("fastProjErrorCutoff", _params.fastProjErrorCutoff);
}

ViewerState::ViewerState(const std::string& filename) : ViewerState() {
  std::ifstream f(filename);
  std::string fieldName;
  bool fieldFound;
  while (f >> fieldName) {
    fieldFound = false;
    for (auto& entry : _entries) {
      if (fieldName == entry->_fieldName) {
        assert(entry->read(f));
        fieldFound = true;
        break;
      }
    }
    assert(fieldFound);
  }
}

void ViewerState::save(const std::string& filename) const {
  std::ofstream f(filename);
  for (auto& entry : _entries) {
    entry->write(f);
  }
}

}