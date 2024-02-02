#ifndef FLOATING_POINT_H
#define FLOATING_POINT_H

#include <cmath>

inline constexpr double THRESHOLD = 1e-10;

inline void round_if_below_threshold(double& v) {
  v = std::abs(v) < THRESHOLD ? 0 : v;
}

#endif 
