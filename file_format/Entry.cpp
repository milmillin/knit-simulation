#include "./Entry.h"

namespace file_format {

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
void Entry<double>::write(std::ofstream& f) {
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

template class Entry<double>;
template class Entry<bool>;
template class Entry<simulator::SimulatorType>;

}  // namepsace file_format
