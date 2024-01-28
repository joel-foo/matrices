#include "EquationSolverTest.h"
#include "MatrixTest.h"
#include "Theorems.h"

int main() {

  // basic matrix functionality
  MatrixTest::run();

  // gaussian, gauss jordan, systems of linear equations
  EquationSolverTest::run();

  // fun theorem stuff
  Theorems::test();

  std::cout << "all tests pass!" << "\n";

  return 0;
}
