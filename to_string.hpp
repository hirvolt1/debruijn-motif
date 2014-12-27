#ifndef TO_STRING_H_
#define TO_STRING_H_
#include <sstream>
#include <string>

template<typename T>
std::string to_string(const T &t){
  std::stringstream ss;
  ss<<t;
  return ss.str();
}
#endif  // TO_STRING_H
