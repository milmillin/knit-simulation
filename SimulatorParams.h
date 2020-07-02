#pragma once

// SimulatorParams contains all parameters for Simulator
struct SimulatorParams {
  float r;
  float m;
  float kLen;
  float kBend;
  float kGlobal;
  float kContact;
  float kDt;
  float kDn;
  float aSmall;
  float aLarge;
  float h;

  static SimulatorParams Default() {
    return SimulatorParams
    {
      0.125f,
      0.006f,
      10000.f,
      0.005f,
      1.5f,
      3250.f,
      0.003f,
      0.03f,
      0.3f,
      0.3f,
      1.f / 11800
    };
  }
};

