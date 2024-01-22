#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <iostream>
#include <initializer_list>
#include <set>
#include <vector>

// typedefs do not work with template paramteters, but using does!
template <typename T> 
using TwoDVector = std::vector<std::vector<T>>;

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
  private:
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

    std::size_t getRows() const {return m_rows;}
    std::size_t getCols() const {return m_cols;}

    // get element at (r,c)
    const T& operator()(int r, int c) const;

    // set element at (r,c)
    T& operator()(int r, int c);

    bool operator==(const Matrix& other) const{
      return m_matrix == other.m_matrix;
    }
};


template <class T>
class IdentityMatrix: public Matrix<T> {
  public:
    IdentityMatrix(int size) 
      : Matrix<T>{size, size} 
    {
      for (std::size_t i {0}; i < this->m_rows; ++i) {
        for (std::size_t j {0}; j < this->m_cols; ++j) {
          if (i == j) 
            this->m_matrix[i][j] = 1;
        }
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
  if (matrix.size() == 0) {
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
  if (this->m_rows != rhs.m_rows) {
    throw std::runtime_error("Number of rows must be equal!");
  }
  if (this->m_cols != rhs.m_cols) {
    throw std::runtime_error("Number of columns must be equal!");
  }
  for (std::size_t i {0}; i < this->m_rows; ++i) {
    for (std::size_t j {0}; j < this->m_cols; ++j) {
      this->m_matrix[i][j] += rhs.m_matrix[i][j];
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
  for (std::size_t i {0}; i < this->m_rows; ++i) {
    for (std::size_t j {0}; j < this->m_cols; ++j) {
      this->m_matrix[i][j] *= rhs;
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
  Matrix<T> result {static_cast<int>(lhs.m_rows), static_cast<int>(rhs.m_cols)};
  for (std::size_t i {0}; i < result.m_rows; ++i) {
    for (std::size_t j {0}; j < result.m_cols; ++j) {
      for (std::size_t k {0}; k < lhs.m_cols; ++k) {
        result.m_matrix[i][j] += lhs.m_matrix[i][k] * rhs.m_matrix[k][j];
      } 
    }
  }
  // other alternatives...
  return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& rhs) {
  if (this->m_rows != rhs.m_rows) {
    throw std::runtime_error("Number of rows must be equal!");
  }
  if (this->m_cols != rhs.m_cols) {
    throw std::runtime_error("Number of columns must be equal!");
  }
  for (std::size_t i {0}; i < this->m_rows; ++i) {
    for (std::size_t j {0}; j < this->m_cols; ++j) {
      this->m_matrix[i][j] -= rhs.m_matrix[i][j];
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
  Matrix<T> result {static_cast<int>(this->m_cols), static_cast<int>(this->m_rows)};
  for (std::size_t i {0}; i < this->m_cols; ++i) {
    for (std::size_t j{0}; j < this->m_rows; ++j) {
      result.m_matrix[i][j] = this->m_matrix[j][i];
    }
  }
  return result;
}

template<class T>
const T& Matrix<T>::operator()(int r, int c) const {
  if (r < 0 or r >= static_cast<int>(m_rows)) {
    throw std::runtime_error("Invalid row");
  }
  if (c < 0 or c >= static_cast<int>(m_cols)) {
    throw std::runtime_error("Invalid column");
  }
  return m_matrix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
}

template<class T>
T& Matrix<T>::operator()(int r, int c) {
  if (r < 0 or r >= static_cast<int>(m_rows)) {
    throw std::runtime_error("Invalid row");
  }
  if (c < 0 or c >= static_cast<int>(m_cols)) {
    throw std::runtime_error("Invalid column");
  }
  return m_matrix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
}

#endif
