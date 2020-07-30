#pragma once

#include <fstream>
#include <iomanip>

#include "../simulator/Enum.h"

namespace file_format {

class IEntry {
public:
  const char* _fieldName;
  virtual void write(std::ofstream& f) = 0;
  virtual bool read(std::ifstream& f) = 0;

protected:
  IEntry(const char* fieldName) : _fieldName(fieldName) { }
};

template<typename T>
class Entry : public IEntry {
public:
  Entry(const char* fieldName, T* value) : IEntry(fieldName), _value(value) { }

  void write(std::ofstream& f) override {
    f << _fieldName << " " << *_value << "\n";
  }

  bool read(std::ifstream& f) override {
    f >> *_value;
    return true;
  }
private:
  T* _value;
};

/*
template<>
void Entry<bool>::write(std::ofstream& f);

template<>
bool Entry<bool>::read(std::ifstream& f);

template<>
void Entry<float>::write(std::ofstream& f);

template<>
void Entry<simulator::SimulatorType>::write(std::ofstream& f);

template<>
bool Entry<simulator::SimulatorType>::read(std::ifstream& f);
*/

template<>
void Entry<bool>::write(std::ofstream& f) {
  f << _fieldName << " " << ((*_value) ? "True" : "False") << "\n";
}

template<>
bool Entry<bool>::read(std::ifstream& f) {
  std::string a;
  f >> a;
  if (a == "True") *_value = true;
  else if (a == "False") *_value = false;
  else return false;
  return true;
}

template<>
void Entry<float>::write(std::ofstream& f) {
  f << _fieldName << " " << std::setprecision(9) << *_value << "\n";
}

template<>
void Entry<simulator::SimulatorType>::write(std::ofstream& f) {
  f << _fieldName << " ";
  switch (*_value)
  {
  case simulator::SimulatorType::Continuous:
    f << "Continuous";
    break;
  case simulator::SimulatorType::Discrete:
    f << "Discrete";
    break;
  default:
    f << "__unsupported";
    break;
  }
  f << "\n";
}

template<>
bool Entry<simulator::SimulatorType>::read(std::ifstream& f) {
  std::string a;
  f >> a;
  if (a == "Continuous") *_value = simulator::SimulatorType::Continuous;
  else if (a == "Discrete") *_value = simulator::SimulatorType::Discrete;
  else return false;
  return true;
}

}  // namespace file_format
