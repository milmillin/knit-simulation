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

template<>
void Entry<bool>::write(std::ofstream& f);

template<>
bool Entry<bool>::read(std::ifstream& f);

template<>
void Entry<double>::write(std::ofstream& f);

template<>
void Entry<simulator::SimulatorType>::write(std::ofstream& f);

template<>
bool Entry<simulator::SimulatorType>::read(std::ifstream& f);

}  // namespace file_format
