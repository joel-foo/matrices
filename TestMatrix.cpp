#include <cassert>
#include <iostream> 
#include <vector>

#include "Matrix.h"

typedef Matrix<int> Mint;

void test_add_matrices() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{2,4,6},{8,10,12},{14,16,18}};
  assert((m1 + m1) == expected);
}

void test_subtract_matrices() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected (3,3,0);
  assert((m1 - m1) == expected);
}

void test_multiply_scalar() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{2,4,6},{8,10,12},{14,16,18}};
  assert((m1 * 2) == expected);
}

void test_multiply_matrices() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  Mint m2 {{9,8,7},{6,5,4},{3,2,1}};
  Mint expected {{30,24,18},{84,69,54},{138,114,90}};
  assert((m1 * m2) == expected);
}

void test_transpose_matrix() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{1,4,7},{2,5,8},{3,6,9}};
  assert(m1.transpose() == expected);
}

void test_get_matrix() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  assert(m1(2,2) == 9);

  bool exceptionThrown = false;
  try {
    m1(3,3);
  } catch (const std::runtime_error& error) {
    exceptionThrown = true;
  }
  assert(exceptionThrown);
}

void test_set_matrix() {
  Mint m1 {{1,2,3},{4,5,6},{7,8,9}};
  m1(2,2) = 10;
  assert(m1(2,2) == 10);

  bool exceptionThrown = false;
  try {
    m1(3,3) = 10;
  } catch (const std::runtime_error& error) {
    exceptionThrown = true;
  }
  assert(exceptionThrown);
}


int main() {
  test_add_matrices();
  test_subtract_matrices();
  test_multiply_scalar();
  test_multiply_matrices();
  test_transpose_matrix();
  test_get_matrix();
  test_set_matrix();
  std::cout << "all tests pass!" << "\n";
}
