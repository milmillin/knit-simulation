#pragma once

namespace simulator {

// SimulatorParams contains all parameters for Simulator
struct SimulatorParams {
  // Debug mode
  bool debug = false;

  double m = 0.006f;
  double kLen = 10000.f;
  // Bending coefficient
  double kBend = 0.05f;
  // Twisting coefficient
  double kTwist = 0.05f;
  // Global damping
  double kGlobal = 1.5f;
  // Contact force coefficient
  double kContact = 3250.f;
  // Contact force damping
  double kDt = 0.003f;
  double kDn = 0.03f;
  // Contact force samples per segment
  int contactForceSamples = 11;
  double aSmall = 0.3f;
  double aLarge = 0.3f;
  // Time delta for each step (time resolution)
  double h = 0.001;
  // Number of steps to run for each button click
  int steps = 100;
  // Gravity acceleration
  double gravity = 9.8f;
  // y coordinate of the ground
  double groundHeight = -5.f;
  // Ground fiction
  double groundFriction = 0.5f;

  double cInit = 1.0f;

  // === Fast projection ===
  // maximum iteration
  int fastProjMaxIter = 20;
  // Early termination when error is small
  double fastProjErrorCutoff = 1e-5;

  static SimulatorParams Default() {
    SimulatorParams param;
    return param;
  }
};

}; // namespace simulator
