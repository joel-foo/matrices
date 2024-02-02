#include <cmath>
#include <set>
#include <unordered_map>
#include <vector>

#include "FloatingPoint.h"
#include "Matrix.h"
#include "Solution.h"
#include "Vectors.h"


namespace linalg {

int compute_max_in_col(int r, int c, const Matrix<double>& m) {
  const TwoDVector<double>& matrix = m.m_matrix;
  double max_val = std::abs(matrix[r][c]);
  int row = r;
  for (int i = r + 1; i < m.m_rows; ++i) {
    if (std::abs(matrix[i][c]) >= max_val) {
      max_val = std::abs(matrix[i][c]);
      row = i; 
    } 
  }
  return row;
}

// * Choice of pivot: " In any case, choosing the largest possible absolute value of the pivot improves the numerical stability of the algorithm, when floating point is used for representing numbers."
Matrix<double> gaussian_elimination(const Matrix<double>& m) {
  auto result = m; //copy assignment is invoked
  auto& matrix = result.m_matrix;
  int r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    // finds the row which gives the max abs value in cth column (from rth row inclusive)
    int row_i = compute_max_in_col(r, c, matrix);
    // if max pivot == 0, move to next col
    if (matrix[row_i][c] == 0) {
      ++c;
      continue;
    }
    std::swap(matrix[r], matrix[row_i]);
    // make everything below result[r][c] 0
    for (int i = r + 1; i < result.m_rows; ++i) {
      double factor = matrix[i][c] / matrix[r][c];
      for (int j = c; j < result.m_cols; ++j) {
        matrix[i][j] -= factor * matrix[r][j];
        round_if_below_threshold(matrix[i][j]);
      }
    }
    ++r;
    ++c;
  }
  return matrix;
}

Matrix<double> gauss_jordan_elimination(const Matrix<double>& m) {
  auto result = gaussian_elimination(m);
  auto& matrix = result.m_matrix;
  int c = 0;
  // make all leading pivot elements 1 
  for (auto it = matrix.begin(); it != matrix.end(); ++it) {
    auto& row = *it;
    while (row[c] == 0) {
      ++c;
    }
    // apply scaling factor
    if (row[c] != 1) {
      double factor = row[c];
      for (int j = c; j < result.m_cols; ++j) {
        row[j] /= factor;
        round_if_below_threshold(row[j]);
      }
    }
    // make all elements above the pivot in each row 0
    auto prev = it;
    while (prev != matrix.begin()) {
      --prev;
      auto& currRow = *prev;
      if (row[c] == 0)
        continue;
      double factor = -currRow[c];
      for (int j = c; j < result.m_cols; ++j) {
        currRow[j] += factor * row[j];
        round_if_below_threshold(currRow[j]);
      }
    }
  }
  return result;
}

int get_rank(const Matrix<double>& m) {
  auto result = gaussian_elimination(m);
  auto& matrix = result.m_matrix;
  int r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    if (matrix[r][c] != 0) {
      ++r;
    }
    ++c;
  }
  return r;
}


TwoDVector<double> get_row_space(const Matrix<double>& m) {
  auto result = gauss_jordan_elimination(m);
  auto& matrix = result.m_matrix;
  int rank = get_rank(m);
  TwoDVector<double> row_space(matrix.begin(), matrix.begin() + rank);
  return row_space;
}

TwoDVector<double> get_col_space(const Matrix<double>& m) {
  TwoDVector<double> col_space;
  auto result = gauss_jordan_elimination(m);
  auto& matrix = result.m_matrix;
  std::vector<int> col_idxs;
  int r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    if (matrix[r][c] != 0) {
      col_idxs.emplace_back(c);
      ++r;
    } 
    ++c;
  }
  for (auto i: col_idxs) {
    std::vector<double> col_vec;
    for (int j = 0; j < result.m_rows; ++j) {
      col_vec.emplace_back(matrix[j][i]);
    }
    col_space.emplace_back(std::move(col_vec));
  }
  return col_space;
}

Matrix<double> inverse(const Matrix<double>& matrix) {
  if (matrix.m_rows != matrix.m_cols) {
    throw std::runtime_error("Number of rows must equal number of columns!");
  }
  if (get_rank(matrix) != matrix.m_rows) {
    throw std::runtime_error("Matrix is not invertible!");
  }
  int rank = matrix.m_rows;
  auto augmented_matrix = matrix.concat(IdentityMatrix<double>(rank));
  auto result = gauss_jordan_elimination(augmented_matrix);
  TwoDVector<double> vec;
  for (int i = 0; i < rank; i++) {
    std::vector<double> v(result.m_matrix[i].begin() + rank, result.m_matrix[i].end());
    vec.emplace_back(std::move(v));
  }
  return vec;
}

// solves the equation Ax = b 
// => given m equations, n unknowns.
// => A is a m x n matrix, x is n x 1, b is a m x 1 column matrix. 
Solution::SystemSolution solve_linear_system(const Matrix<double>& A, const Matrix<double>& b) {
  // a special case: assuming A is nxn, if A is invertible => x = A^-1(b), x is unique IFF A is invertible (t)
  // proof: https://www.quora.com/How-do-I-show-that-AX-B-has-a-unique-solution-if-and-only-if-Matrix-A-is-invertible.
  using namespace Solution;
  if (A.m_rows != b.m_rows) {
    throw std::runtime_error("A and b need to have the same number of rows");
  }
  auto augmented_matrix = A.concat(b);
  auto RREF = gauss_jordan_elimination(augmented_matrix);
  int num_variables = A.m_cols;
  SolutionType solutionType = get_solution_type(RREF);
  switch (solutionType) {
    case SolutionType::NO_SOLUTION:
      return SystemSolution{.type = solutionType};
    case SolutionType::ONE_SOLUTION:
    {
      std::vector<VariableSolution> solutions;
      for (int i = 0; i < num_variables; ++i) {
        solutions.emplace_back(VariableSolution{.val=RREF.m_matrix[i][RREF.m_cols - 1]});
      }
      return SystemSolution{.type = solutionType, .m_solutions=solutions};
    }
    case SolutionType::INFINITELY_MANY_SOLUTIONS:
    {
      int r = 0, c = 0;
      auto& matrix = RREF.m_matrix;

      std::set<int> pivot_cols;
      std::set<int> non_pivot_cols;
      std::unordered_map<int, int> pivot_col_to_row_map;

      while (r < RREF.m_rows && c < RREF.m_cols - 1) {
        if (matrix[r][c] != 0) {
          pivot_cols.insert(c);
          pivot_col_to_row_map[c] = r;
          ++r;
        }
        ++c;
      }

      for (int i = 0; i < RREF.m_cols - 1; ++i) {
        if (!pivot_cols.contains(i)) {
          non_pivot_cols.insert(i);
        }
      }

      std::unordered_map<int, VariableSolution> solutionMap; 
      std::vector<VariableSolution> solutions;

      // iterate backwards - bottom up DP :)
      for (auto it = pivot_cols.rbegin(); it != pivot_cols.rend(); ++it) {
        int currCol = *it; 
        int currRow = pivot_col_to_row_map[currCol];
        double rhsVal = matrix[currRow][RREF.m_cols - 1];
        //value component
        for (int c = currCol + 1; c < RREF.m_cols - 1; ++c) {
          if (matrix[currRow][c] == 0) continue;
          // if non pivot, easy! just deduct the free variable.
          auto& currMap = solutionMap[currCol].variable_count_map;
          if (non_pivot_cols.contains(c)) {
            currMap[c] -= matrix[currRow][c];
          } else {
            // if pivot, do some math. 
            solutionMap[currCol].val -= solutionMap[c].val * matrix[currRow][c];
            for (const auto& [key, val]: currMap) {
              if (val == 0) continue;
              currMap[key] -= val * matrix[currRow][c];
            }
          }
        }
        solutionMap[currCol].val += rhsVal;
      }
      for (int i = 0; i < RREF.m_cols - 1; ++i) {
        solutions.emplace_back(std::move(solutionMap[i]));
      } 
      return SystemSolution{.type = solutionType, .m_solutions=solutions, .free_variables=non_pivot_cols, .num_free_variables=static_cast<int>(non_pivot_cols.size())};
    }
  }
}

// A is a REF or RREF form. 
Solution::SolutionType get_solution_type(const Matrix<double>& A) {
  using namespace Solution;
  // inconsistent if if the last column of a row-echelon form of the augmented matrix is a pivot column
  int r = 0, c = 0;
  auto& matrix = A.m_matrix;
  while (r < A.m_rows && c < A.m_cols) {
    if (matrix[r][c] != 0) {
      if (c == A.m_cols - 1) return SolutionType::NO_SOLUTION;
      ++r;
    }
    ++c;
  }
  int num_variables = A.m_cols - 1; 
  int non_zero_rows = r;
  // only one solution -> every col is a pivot col except the last
  // OR alternatively, a consistent linear system has only one solution if the number of variables in the linear system is equal to the number of nonzero rows in a row-echelon form of the augmented matrix.
  // infinitely many solutions -> at least non pivot col that is not the last col
  return num_variables == non_zero_rows ? SolutionType::ONE_SOLUTION: SolutionType::INFINITELY_MANY_SOLUTIONS;
}

}
