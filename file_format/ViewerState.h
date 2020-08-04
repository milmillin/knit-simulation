#pragma once

#include <string>
#include <vector>
#include <memory>

#include "simulator/Enum.h"
#include "simulator/SimulatorParams.h"
#include "./Entry.h"

namespace file_format {

class ViewerState {
public:
  ViewerState(simulator::SimulatorType type = simulator::SimulatorType::Continuous,
    simulator::SimulatorParams params = simulator::SimulatorParams(),
    int numSteps = 1);
  ViewerState(const std::string& filename);
  void save(const std::string& filename) const;

  const simulator::SimulatorType& getType() const { return _type; }
  const simulator::SimulatorParams& getParams() const { return _params; }
  int getNumSteps() const { return _numSteps; }
private:
  simulator::SimulatorType _type;
  simulator::SimulatorParams _params;
  int _numSteps;

  std::vector<std::unique_ptr<IEntry>> _entries;
};

} // namespace file_format


