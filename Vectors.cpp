#include "Vectors.h"

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
