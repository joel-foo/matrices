#ifndef TWO_D_VECTOR_H
#define TWO_D_VECTOR_H

#include <iostream>
#include <vector>

#include "FloatingPoint.h"

// typedefs do not work with template paramteters, but using does!
template <typename T> 
using TwoDVector = std::vector<std::vector<T>>;

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
  out << "(";
  int len = vec.size();
  for (int i = 0; i < len; ++i) {
    out << vec[i];
    if (i != len - 1) {
      out << ", ";
    }
  }
  out << ")";
  return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const TwoDVector<T>& vec) {
  int rows = vec.size();
  out << "(";
  for (int i = 0; i < rows; ++i) {
    out << vec[i];
    if (i != rows - 1) {
      out << ", ";
    }
  }
  out << ")";
  return out;
}

bool operator==(const std::vector<double>& v1, const std::vector<double>& v2);
bool operator==(const TwoDVector<double> &v1, const TwoDVector<double> &v2);

#endif
