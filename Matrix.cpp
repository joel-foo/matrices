#include "FloatingPoint.h"
#include "Matrix.h"
#include "Solution.h"

std::size_t maxCol(std::size_t r, std::size_t c, const Matrix<double>& m) {
  const TwoDVector<double>& matrix = m.m_matrix;
  double max_val = std::abs(matrix[r][c]);
  std::size_t row = r;
  for (std::size_t i = r + 1; i < m.m_rows; i++) {
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
  std::size_t r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    // finds the row which gives the max abs value in cth column (from rth row inclusive)
    std::size_t row_i = maxCol(r, c, matrix);
    // if max pivot == 0, move to next col
    if (matrix[row_i][c] == 0) {
      ++c;
      continue;
    }
    std::swap(matrix[r], matrix[row_i]);
    // make everything below result[r][c] 0
    for (std::size_t i = r + 1; i < result.m_rows; ++i) {
      double factor = matrix[i][c] / matrix[r][c];
      for (std::size_t j = c; j < result.m_cols; ++j) {
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
  std::size_t c = 0;
  // make all leading pivot elements 1 
  for (auto it = matrix.begin(); it != matrix.end(); ++it) {
    auto& row = *it;
    while (row[c] == 0) {
      ++c;
    }
    // apply scaling factor
    if (row[c] != 1) {
      double factor = row[c];
      for(std::size_t j = c; j < result.m_cols; ++j) {
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
      for (std::size_t j = c; j < result.m_cols; ++j) {
        currRow[j] += factor * row[j];
        round_if_below_threshold(currRow[j]);
      }
    }
  }
  return result;
}

int getRank(const Matrix<double>& m) {
  auto result = gaussian_elimination(m);
  auto& matrix = result.m_matrix;
  std::size_t r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    if (matrix[r][c] != 0) {
      ++r;
    }
    ++c;
  }
  return static_cast<int>(r);
}


TwoDVector<double> getRowSpace(const Matrix<double>& m) {
  auto result = gauss_jordan_elimination(m);
  auto& matrix = result.m_matrix;
  int rank = getRank(m);
  TwoDVector<double> rowSpace(matrix.begin(), matrix.begin() + rank);
  return rowSpace;
}

TwoDVector<double> getColSpace(const Matrix<double>& m) {
  TwoDVector<double> colSpace;
  auto result = gauss_jordan_elimination(m);
  auto& matrix = result.m_matrix;
  std::vector<std::size_t> colIdxs;
  std::size_t r = 0, c = 0;
  while (r < result.m_rows && c < result.m_cols) {
    if (matrix[r][c] != 0) {
      colIdxs.emplace_back(c);
      ++r;
    } 
    ++c;
  }
  for (auto i: colIdxs) {
    std::vector<double> colVec;
    for(std::size_t j = 0; j < result.m_rows; ++j) {
      colVec.emplace_back(matrix[j][i]);
    }
    colSpace.emplace_back(std::move(colVec));
  }
  return colSpace;
}

Matrix<double> inverse(const Matrix<double>& matrix) {
  if(matrix.m_rows != matrix.m_cols) {
    throw std::runtime_error("Number of rows must equal number of columns!");
  }
  if(getRank(matrix) != matrix.getRows()) {
    throw std::runtime_error("Matrix is not invertible!");
  }
  std::size_t rank = matrix.m_rows;
  auto augmentedMatrix = matrix.concat(IdentityMatrix<double>(static_cast<int>(rank)));
  auto result = gauss_jordan_elimination(augmentedMatrix);
  TwoDVector<double> vec;
  for(std::size_t i = 0; i < rank; i++) {
    std::vector<double> v(result.m_matrix[i].begin() + static_cast<long>(rank), result.m_matrix[i].end());
    vec.emplace_back(std::move(v));
  }
  return vec;
}

// solves the equation Ax = b 
// => given m equations, n unknowns.
// => A is a m x n matrix, x is n x 1, b is a m x 1 column matrix. 
Solution solve_linear_system(const Matrix<double>& A, const Matrix<double>& b) {
  // a special case: assuming A is nxn, if A is invertible => x = A^-1(b), x is unique IFF A is invertible (t)
  // proof: https://www.quora.com/How-do-I-show-that-AX-B-has-a-unique-solution-if-and-only-if-Matrix-A-is-invertible.
  auto augmentedMatrix = A.concat(b);
  auto RREF = gauss_jordan_elimination(augmentedMatrix);
  std::size_t numVariables = A.m_cols;
  SolutionType solutionType = get_solution_type(RREF);
  switch (solutionType) {
    case SolutionType::NO_SOLUTION:
      return Solution{.type = solutionType};
    case SolutionType::ONE_SOLUTION:
    {
      std::vector<solutionPair> solutions;
      for (std::size_t i = 0; i < numVariables; ++i) {
        solutions.emplace_back(std::make_pair(RREF.m_matrix[i][RREF.m_cols - 1], std::nullopt));
      }
      return Solution{.type = solutionType, .m_solutions=solutions};
    }
    case SolutionType::INFINITELY_MANY_SOLUTIONS:
    {
      std::size_t r = 0, c = 0;
      auto& matrix = RREF.m_matrix;

      std::set<std::size_t> pivot_cols;
      std::set<std::size_t> non_pivot_cols;
      std::unordered_map<std::size_t, std::size_t> pivotColToRowMap;

      while (r < RREF.m_rows && c < RREF.m_cols - 1) {
        if (matrix[r][c] != 0) {
          pivot_cols.insert(c);
          pivotColToRowMap[c] = r;
          ++r;
        }
        ++c;
      }

      for(std::size_t i = 0; i < RREF.m_cols - 1; ++i) {
        if (pivot_cols.find(i) == pivot_cols.end()) {
          non_pivot_cols.insert(i);
        }
      }

      std::unordered_map<std::size_t, solutionPair> solutionMap; 
      std::vector<solutionPair> solutions;

      // iterate backwards - bottom up DP :)
      for (auto it = pivot_cols.rbegin(); it != pivot_cols.rend(); ++it) {
        auto currCol = *it; 
        auto currRow = pivotColToRowMap[currCol];
        auto rhsVal = matrix[currRow][RREF.m_cols - 1];
        //value component
        for(std::size_t c = currCol + 1; c < RREF.m_cols - 1; ++c) {
          if (matrix[currRow][c] == 0) continue;
          // if non pivot, easy! just deduct the free variable.
          auto& opt = solutionMap[currCol].second;
          if (!opt.has_value()) {
            opt.emplace(std::unordered_map<std::size_t, double>());
          }
          auto& currMap = opt.value();
          if (non_pivot_cols.find(c) != non_pivot_cols.end()) {
            currMap[c] -= matrix[currRow][c];
          } else {
            // if pivot, do some math. 
            solutionMap[currCol].first -= solutionMap[c].first * matrix[currRow][c];
            for(const auto& [key, val]: currMap) {
              if (val == 0) continue;
              currMap[key] -= val * matrix[currRow][c];
            }
          }
        }
        solutionMap[currCol].first += rhsVal;
      }
      for(std::size_t i = 0; i < RREF.m_cols - 1; ++i) {
        solutions.emplace_back(solutionMap[i]);
      } 
      return Solution{.type = solutionType, .m_solutions=solutions, .freeVariables=non_pivot_cols, .num_free_variables=static_cast<int>(non_pivot_cols.size())};
    }
  }
}

// A is a REF or RREF form. 
SolutionType get_solution_type(const Matrix<double>& A) {
  // inconsistent if if the last column of a row-echelon form of the augmented matrix is a pivot column
  std::size_t r = 0, c = 0;
  auto& matrix = A.m_matrix;
  while (r < A.m_rows && c < A.m_cols) {
    if (matrix[r][c] != 0) {
      if (c == A.m_cols - 1) return SolutionType::NO_SOLUTION;
      ++r;
    }
    ++c;
  }
  std::size_t numVariables = A.m_cols - 1; 
  std::size_t nonZeroRows = r;
  // only one solution -> every col is a pivot col except the last
  // OR alternatively, a consistent linear system has only one solution if the number of variables in the linear system is equal to the number of nonzero rows in a row-echelon form of the augmented matrix.
  // infinitely many solutions -> at least non pivot col that is not the last col
  return numVariables == nonZeroRows ? SolutionType::ONE_SOLUTION: SolutionType::INFINITELY_MANY_SOLUTIONS;
}

