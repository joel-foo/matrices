#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <initializer_list>
#include <set>
#include <vector>

#include "Solution.h"
#include "Vectors.h"

template <typename T> 
using TwoDList = std::initializer_list<std::initializer_list<T>>;

//forward declare class
template <class T> 
class Matrix; 

//predeclare friend templates
template <class T> 
std::ostream& operator<<(std::ostream& out, const Matrix<T>& m);

template <class T> 
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs);

template <class T> 
Matrix<T> operator*(const Matrix<T>& lhs, const int rhs);

template <class T> 
Matrix<T> operator*(const int lhs, const Matrix<T>& rhs);

template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs);

template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs);

template <class T> 
class Matrix {
  protected:
    std::size_t m_rows;
    std::size_t m_cols;
    TwoDVector<T> m_matrix;

  public: 
    Matrix(int rows, int cols, T init=0);
    Matrix(const TwoDVector<T>& matrix);
    Matrix(const TwoDList<T>& list);

    // why need <>? : https://isocpp.org/wiki/faq/templates#template-friends
    friend std::ostream& operator<< <>(std::ostream& out, const Matrix<T>& m);

    // https://stackoverflow.com/questions/4421706/what-are-the-basic-rules-and-idioms-for-operator-overloading
    //addition
    Matrix& operator+=(const Matrix<T>& rhs);
    friend Matrix operator+ <>(const Matrix<T>& lhs, const Matrix<T>& rhs);

    //scalar multiplication
    Matrix& operator*=(const int rhs);
    friend Matrix operator* <>(const Matrix<T>& lhs, const int rhs);
    friend Matrix operator* <>(const int lhs, const Matrix<T>& rhs);

    //matrix multiplication
    Matrix& operator*=(const Matrix<T>& rhs);
    friend Matrix operator* <>(const Matrix<T>& lhs, const Matrix<T>& rhs);

    //subtraction
    Matrix& operator-=(const Matrix<T>& rhs);
    friend Matrix operator- <>(const Matrix<T>& lhs, const Matrix<T>& rhs);

    //transpose
    Matrix transpose() const;

    //concat 
    //axis = 0 -> horizontal; 1 -> vertical
    Matrix concat(const Matrix<T>& rhs, int axis=0) const;

    //flatten 1D
    std::vector<T> flatten() const;

    int getRows() const {return static_cast<int>(m_rows);}
    int getCols() const {return static_cast<int>(m_cols);}

    // get element at (r,c)
    const T& operator()(int r, int c) const;

    // set element at (r,c)
    T& operator()(int r, int c);

    bool operator==(const Matrix& other) const{
      return m_matrix == other.m_matrix;
    }

    //friend non member functions below:
    friend std::size_t maxCol(std::size_t r, std::size_t c, const Matrix<double>& matrix);

    //produces a matrix of Row Echelon Form (REF), note: all zero rows are guaranteed to be at bottom of both the REF and RREF. 
    friend Matrix<double> gaussian_elimination(const Matrix<double>& matrix);

    //produces a matrix of Reduced Row Echelon Form (RREF). note: RREF is unique while REF is not. 
    friend Matrix<double> gauss_jordan_elimination(const Matrix<double>& matrix);

    friend int getRank(const Matrix<double>& m);

    // row space and col space are obtained from RREF for simplicity. 
    friend TwoDVector<double> getRowSpace(const Matrix<double>& m);
    friend TwoDVector<double> getColSpace(const Matrix<double>& m);

    friend Matrix<double> inverse(const Matrix<double>& m);

    friend SolutionType get_solution_type(const Matrix<double>& A);

    friend Solution solve_linear_system(const Matrix<double>& A, const Matrix<double>& b);
};


template <class T>
class IdentityMatrix: public Matrix<T> {
  public:
    IdentityMatrix(int size) 
      : Matrix<T>{size, size} 
    {
      for (std::size_t i = 0; i < this->m_rows; ++i) {
        this->m_matrix[i][i] = 1;
      }
    }
};

template<class T>
Matrix<T>::Matrix(int rows, int cols, T init) {
  if (rows <= 0) {
    throw std::runtime_error("rows must be > 0!");
  }
  if (cols <= 0) {
      throw std::runtime_error("rows must be > 0!");
  }
  m_rows = static_cast<std::size_t>(rows);
  m_cols = static_cast<std::size_t>(cols);
  m_matrix = TwoDVector<T>(m_rows, std::vector<T>(m_cols, init));
}

template<class T>
Matrix<T>::Matrix(const TwoDVector<T>& matrix) {
  if (matrix.empty()) {
    throw std::runtime_error("Matrix cannot be empty");
  }
  std::set<std::size_t> cols;
  for (const auto& row: matrix) {
    cols.insert(row.size());
    if (cols.size() > 1) {
      throw std::runtime_error("Inconsistent number of columns in matrix!");
    }
  }
  m_rows = matrix.size();
  m_cols = matrix[0].size();
  m_matrix = matrix; //invokes copy assignment
}

template<class T>
Matrix<T>::Matrix(const TwoDList<T>& matrix) {
  if (matrix.size() == 0) {
    throw std::runtime_error("Matrix cannot be empty");
  }
  std::set<std::size_t> cols;
  for(const auto& row: matrix) {
    cols.insert(row.size());
    if (cols.size() > 1) {
      throw std::runtime_error("Inconsistent number of columns in matrix!");
    }
    m_matrix.emplace_back(row);
  }
  m_rows = matrix.size();
  m_cols = *(cols.begin());
}

template<class T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& m) {
  for (const auto& row: m.m_matrix) {
    for (const auto& v: row) {
      out << v << " ";
    }
    out << "\n";
  }
  return out;
}

//member functions return *this so we can chain multiple operations
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix& rhs) {
  if (m_rows != rhs.m_rows) {
    throw std::runtime_error("Number of rows must be equal!");
  }
  if (m_cols != rhs.m_cols) {
    throw std::runtime_error("Number of columns must be equal!");
  }
  for (std::size_t i = 0; i < m_rows; ++i) {
    for (std::size_t j = 0; j < m_cols; ++j) {
      m_matrix[i][j] += rhs.m_matrix[i][j];
    }
  }
  return *this;
}

// implement + in terms of += 
template<class T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
  Matrix<T> result = lhs;
  result += rhs;
  return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const int rhs) {
  for (std::size_t i = 0; i < m_rows; ++i) {
    for (std::size_t j = 0; j < m_cols; ++j) {
      m_matrix[i][j] *= rhs;
    }
  }
  return *this;
}

template<class T>
Matrix<T> operator*(const Matrix<T>& lhs, const int rhs) {  
  Matrix<T> result = lhs;
  result *= rhs;
  return result;
}

template<class T>
Matrix<T> operator*(const int lhs, const Matrix<T>& rhs) {
  return rhs * lhs;
}

// here, we implement *= using * instead, since * creates a matrix of new dimensions. 
template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
  *this = *this * rhs;
  return *this;
}

template<class T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
  if (lhs.m_cols != rhs.m_rows) {
    throw std::runtime_error("Rows and cols dimension mismatch!");
  }  
  // normal method
  Matrix<T> result {lhs.getRows(), rhs.getCols()};
  for (std::size_t i = 0; i < result.m_rows; ++i) {
    for (std::size_t j = 0; j < result.m_cols; ++j) {
      for (std::size_t k = 0; k < lhs.m_cols; ++k) {
        result.m_matrix[i][j] += lhs.m_matrix[i][k] * rhs.m_matrix[k][j];
      } 
    }
  }
  // other alternatives...
  return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& rhs) {
  if (m_rows != rhs.m_rows) {
    throw std::runtime_error("Number of rows must be equal!");
  }
  if (m_cols != rhs.m_cols) {
    throw std::runtime_error("Number of columns must be equal!");
  }
  for (std::size_t i = 0; i < m_rows; ++i) {
    for (std::size_t j = 0; j < m_cols; ++j) {
      m_matrix[i][j] -= rhs.m_matrix[i][j];
    }
  }
  return *this;
}

// implement - in terms of -= 
template<class T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
  Matrix<T> result = lhs;
  result -= rhs;
  return result;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> result {getCols(), getRows()};
  for (std::size_t i = 0; i < m_cols; ++i) {
    for (std::size_t j = 0; j < m_rows; ++j) {
      result.m_matrix[i][j] = m_matrix[j][i];
    }
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::concat(const Matrix<T>& rhs, int axis) const {
  // horizontal
  if (!axis) {
    if (m_rows != rhs.m_rows) {
      throw std::runtime_error("Fail to join horizontally due to mismatched row number");
    }
    Matrix<T> result = *this;
    for (std::size_t i = 0; i < m_rows; ++i) {
      auto& leftRow = result.m_matrix[i], rightRow = rhs.m_matrix[i];
      leftRow.insert(leftRow.end(), rightRow.begin(), rightRow.end());
    }
    result.m_cols += rhs.m_cols;
    return result;
  } else {
     // vertical
    if (m_cols != rhs.m_cols) {
      throw std::runtime_error("Fail to join vertically due to mismatched col number");
    }
    Matrix<T> result = *this;
    auto& leftMatrix = result.m_matrix, rightMatrix = rhs.m_matrix;
    leftMatrix.insert(leftMatrix.end(), rightMatrix.begin(), rightMatrix.end());
    result.m_rows += rhs.m_rows;
    return result;
  }
}

template<class T>
std::vector<T> Matrix<T>::flatten() const {
  std::vector<T> vec;
  for(const auto& row: m_matrix) {
    vec.insert(vec.end(), row.begin(), row.end());
  }
  return vec;
}

template<class T>
const T& Matrix<T>::operator()(int r, int c) const {
  if (r < 0 or r >= getRows()) {
    throw std::runtime_error("Invalid row");
  }
  if (c < 0 or c >= getCols()) {
    throw std::runtime_error("Invalid column");
  }
  return m_matrix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
}

template<class T>
T& Matrix<T>::operator()(int r, int c) {
  if (r < 0 or r >= getRows()) {
    throw std::runtime_error("Invalid row");
  }
  if (c < 0 or c >= getCols()) {
    throw std::runtime_error("Invalid column");
  }
  return m_matrix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
}

#endif
