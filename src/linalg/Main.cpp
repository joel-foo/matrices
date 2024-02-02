#include <iostream>

#include "EquationSolverTest.h"
#include "MatrixTest.h"
#include "Theorems.h"

int main() {

  // basic matrix functionality
  MatrixTest::run();

  // gaussian, gauss jordan, systems of linear equations
  EquationSolverTest::run();

  linalg::Matrix<double> m (2,3);
  std::cout << linalg::get_rank(m) << "\n";

  // fun theorem stuff
  Theorems::test();

  std::cout << "all tests pass!" << "\n";

  return 0;
}
