#include <cmath>
#include <iomanip>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "Solution.h"

std::ostream& operator<< (std::ostream& out, Solution::SystemSolution s) {
  //to 3sf, if necessary.
  out << std::setprecision(3);
  if (!s.contains_free_variables()) {
    int variable_count = 1;
    for (std::size_t i = 0; i < s.m_solutions.size(); ++i) {
      out << 'x' + std::to_string(variable_count++) << ": " << s.m_solutions[i].val;
      if (i != s.m_solutions.size() - 1) {
        out << ", ";
      }
    }
    return out;
  }
  std::unordered_map<int, std::string> free_variable_map;
  int free_variable_count = 1;
  // actual parameterizations of the free variables i.e. a1, a2,...
  std::vector<std::string> free_variables;
  for (const auto& var: s.free_variables) {
    free_variables.emplace_back('a' + std::to_string(free_variable_count++));
    free_variable_map[var] = free_variables.back();
  }
  int variable_count = 1;
  for (std::size_t i = 0; i < s.m_solutions.size(); ++i) {
    out << 'x' + std::to_string(variable_count++) << ":";
    if (s.free_variables.contains(i)) {
      out << free_variable_map[i];
    } else {
      auto& solution = s.m_solutions[i];
      if (solution.val != 0) {
        out << solution.val;
      }
      if (!solution.contains_free_variables()) continue;
      auto& map = solution.variable_count_map;
      for (const auto& var: s.free_variables) {
        if (map[var] == 0) continue;
        double val = map[var];
        out << (val > 0 ? "+": "-");
        val = std::abs(val);
        if (val != 1) {
          out << val;
        }
        out << free_variable_map[var];
      }
    }
    if (i != s.m_solutions.size() - 1) {
      out << ", ";
    }
  }
  out << " where ";
  for (std::size_t i = 0; i < free_variables.size(); ++i) {
    out << free_variables[i];
    if (i != free_variables.size() - 1) {
      out << ", ";
    }
  }
  out << " are arbitrary parameters";
  return out;
}
