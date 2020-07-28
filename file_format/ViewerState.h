#pragma once

#include <string>
#include <vector>
#include <memory>

#include "./Entry.h"
#include "../simulator/Enum.h"
#include "../simulator/SimulatorParams.h"

namespace file_format {

class ViewerState {
public:
  ViewerState(simulator::SimulatorType type = simulator::SimulatorType::Continuous,
    simulator::SimulatorParams params = simulator::SimulatorParams(),
    int numSteps = 0);
  ViewerState(const std::string& filename);
  void save(const std::string& filename) const;
private:
  simulator::SimulatorType _type;
  simulator::SimulatorParams _params;
  int _numSteps;

  std::vector<std::unique_ptr<IEntry>> _entries;
};

} // namespace file_format


