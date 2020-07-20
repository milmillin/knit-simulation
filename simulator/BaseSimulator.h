#pragma once

#include "../file_format/yarnRepr.h"
#include "./SimulatorParams.h"
#include "./AABB.h"

namespace simulator
{

class BaseSimulator {
 public:
  // Returns current yarns
  virtual const file_format::YarnRepr &getYarns() const { return this->yarns; };

  // Simulates next timestep.
  virtual void step() = 0;

  // Simulation parameters
  SimulatorParams params;

 protected:
  // Position and meta-data
  file_format::YarnRepr yarns;
};

} // namespace simulator
