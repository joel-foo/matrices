#include <cmath>
#include <iomanip>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "Solution.h"

std::ostream& operator<< (std::ostream& out, Solution s) {
  //to 3sf, if necessary.
  out << std::setprecision(3);
  if (s.freeVariables.empty()) {
    int variableCount = 1;
    for (std::size_t i = 0; i < s.m_solutions.size(); ++i) {
      out << 'x' + std::to_string(variableCount++) << ": " << s.m_solutions[i].first;
      if (i != s.m_solutions.size() - 1) {
        out << ", ";
      }
    }
    return out;
  }
  std::unordered_map<int, std::string> freeVariableToValMap;
  int freeVariableCount = 1;
  // actual parameterizations of the free variables i.e. a1, a2,...
  std::vector<std::string> freeVariableVec;
  for (const auto& var: s.freeVariables) {
    freeVariableVec.emplace_back('a' + std::to_string(freeVariableCount));
    freeVariableToValMap[var] = freeVariableVec.back();
    ++freeVariableCount;
  }
  int variableCount = 1;
  for (std::size_t i = 0; i < s.m_solutions.size(); ++i) {
    out << 'x' + std::to_string(variableCount++) << ":";
    if (s.freeVariables.contains(i)) {
      out << freeVariableToValMap[i];
    } else {
      auto& solution = s.m_solutions[i];
      if (solution.first != 0) {
        out << solution.first;
      }
      if (!solution.second.has_value()) continue;
      auto map = solution.second.value();
      for (const auto& var: s.freeVariables) {
        if (map[var] == 0) continue;
        double val = map[var];
        out << (val > 0 ? "+": "-");
        val = std::abs(val);
        if (val != 1) {
          out << val;
        }
        out << freeVariableToValMap[var];
      }
    }
    if (i != s.m_solutions.size() - 1) {
      out << ", ";
    }
  }
  out << " where ";
  for (std::size_t i = 0; i < freeVariableVec.size(); ++i) {
    out << freeVariableVec[i];
    if (i != freeVariableVec.size() - 1) {
      out << ", ";
    }
  }
  out << " are arbitrary parameters";
  return out;
}
