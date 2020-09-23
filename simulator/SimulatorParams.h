#pragma once

namespace simulator {

// SimulatorParams contains all parameters for Simulator
struct SimulatorParams {
  // Debug mode
  bool debug = false;
  // Show statistics
  bool statistics = true;
  // Enable ground collision
  bool enableGround = false;
  // Enable length constrain;
  bool enableLengthConstrain = true;

  double m = 0.006;
  double kLen = 10;
  // Bending coefficient
  double kBend = 0.0;
  // Allowed tolerance when updating the material frame
  double materialFrameTolerance = 0.01;
  // Maximum update iteration then updating the material frame
  int materialFrameMaxUpdate = 20;
  // Twisting coefficient
  double kTwist = 0.05;
  // Global damping
  double kGlobalDamping = 0.00001;
  // Contact force coefficient
  double kContact = 3250;
  // Contact force damping
  double kDt = 0.003;
  double kDn = 0.03;
  // Contact force samples per segment
  int contactForceSamples = 11;
  double aSmall = 0.3;
  double aLarge = 0.3;
  // Allowed deviation for linearized contact force approximation
  double contactModelTolerance = 0.1;
  // Allowed update interval (steps) between linearized model update
  int maxContactModelUpdateInterval = 100;
  // Time delta for each step (time resolution)
  double h = 0.001;
  // Number of steps to run for each button click
  int steps = 100;
  // Gravity acceleration
  double gravity = 0.0;
  // y coordinate of the ground
  double groundHeight = -5.0;
  // Ground fiction
  double groundFriction = 0.5;

  double cInit = 0.5;

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
