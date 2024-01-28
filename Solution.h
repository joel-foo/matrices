#ifndef SOLUTION_H
#define SOLUTION_H 

#include <optional>
#include <set>
#include <sstream> 
#include <string>
#include <unordered_map>
#include <utility>

#include "Vectors.h"

enum class SolutionType {
  NO_SOLUTION,
  ONE_SOLUTION,
  INFINITELY_MANY_SOLUTIONS
};

typedef std::pair<double, std::optional<std::unordered_map<int, double>>> solutionPair; 

struct Solution;

std::ostream& operator<< (std::ostream& out, Solution s);

struct Solution {
  SolutionType type;
  std::vector<solutionPair> m_solutions;
  std::set<int> freeVariables;
  int num_free_variables;

  auto get_compute_function() {
    return [&](std::initializer_list<double> lst) -> std::vector<double> {
      if (lst.size() != freeVariables.size()) {
        throw std::runtime_error("Passed in incorrect number of values for free variables. You need to pass in " + std::to_string(freeVariables.size()) + " values!");
      }
      // create mapping from list passed in to free variables
      std::unordered_map<int, double> freeVariableMap; 
      auto setIt = freeVariables.begin();
      for (auto it = lst.begin(); it != lst.end(); ++it) {
        freeVariableMap[*setIt] = *it;
        ++setIt;
      }
      std::vector<double> realSolution;
      for (std::size_t i = 0; i < m_solutions.size(); ++i) {
        if (freeVariables.find(i) != freeVariables.end()) {
          realSolution.emplace_back(freeVariableMap[i]);
          continue;
        }
        auto& solution = m_solutions[i];
        if (!solution.second.has_value()) {
          realSolution.emplace_back(solution.first);
        } else {
          auto map = solution.second.value();
          double acc = solution.first;
          for (const auto& [key, val]: map) {
            acc += val * freeVariableMap[key];
          }
          realSolution.emplace_back(acc);
        }
      }
      return realSolution;
    };
  }

  std::string toString() const {
    std::stringstream ss;
    ss << (*this);
    return ss.str();
  }
};

#endif
