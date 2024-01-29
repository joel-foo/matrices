#ifndef SOLUTION_H
#define SOLUTION_H 

#include <set>
#include <sstream> 
#include <string>
#include <unordered_map>

#include "Vectors.h"

enum class SolutionType {
  NO_SOLUTION,
  ONE_SOLUTION,
  INFINITELY_MANY_SOLUTIONS
};

struct SystemSolution;

std::ostream& operator<< (std::ostream& out, SystemSolution s);

struct VariableSolution {
  double val;
  std::unordered_map<int, double> variable_count_map; 
  
  bool contains_free_variables() const {
    return !variable_count_map.empty();
  }
};

struct SystemSolution {
  SolutionType type;
  std::vector<VariableSolution> m_solutions;
  std::set<int> free_variables;
  int num_free_variables;

  auto get_compute_function() {
    return [&](std::initializer_list<double> lst) -> std::vector<double> {
      if (lst.size() != free_variables.size()) {
        throw std::runtime_error("Passed in incorrect number of values for free variables. You need to pass in " + std::to_string(free_variables.size()) + " values!");
      }
      // create mapping from list passed in to free variables
      std::unordered_map<int, double> free_variable_val_map; 
      auto setIt = free_variables.begin();
      for (auto it = lst.begin(); it != lst.end(); ++it) {
        free_variable_val_map[*setIt] = *it;
        ++setIt;
      }
      std::vector<double> res;
      for (std::size_t i = 0; i < m_solutions.size(); ++i) {
        if (free_variables.contains(i)) {
          res.emplace_back(free_variable_val_map[i]);
          continue;
        }
        auto& solution = m_solutions[i];
        if (!solution.contains_free_variables()) {
          res.emplace_back(solution.val);
        } else {
          auto& map = solution.variable_count_map;
          double acc = solution.val;
          for (const auto& [key, val]: map) {
            acc += val * free_variable_val_map[key];
          }
          res.emplace_back(acc);
        }
      }
      return res;
    };
  }

  bool contains_free_variables() const {
    return !free_variables.empty();
  };

  std::string toString() const {
    std::stringstream ss;
    ss << (*this);
    return ss.str();
  }
};

#endif
