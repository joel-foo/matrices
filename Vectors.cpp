#include "Vectors.h"

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
  out << "(";
  std::size_t len = vec.size();
  for(std::size_t i = 0; i < len; ++i) {
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
  std::size_t rows = vec.size();
  out << "(";
  for(std::size_t i = 0; i < rows; ++i) {
    out << vec[i];
    if (i != rows - 1) {
      out << ", ";
    }
  }
  out << ")";
  return out;
}

bool operator==(const std::vector<double>& v1, const std::vector<double>& v2){
    return std::equal(v1.begin(), v1.end(), v2.begin(),
                      [](double value1, double value2)
                      {
                          constexpr double epsilon = THRESHOLD; // Choose whatever you need here
                          return std::fabs(value1 - value2) < epsilon;
                      });
}


bool operator==(const TwoDVector<double> &v1, const TwoDVector<double> &v2) {
  return std::equal(v1.begin(), v1.end(), v2.begin(), [](auto p1, auto p2) {return p1 == p2;});
}
