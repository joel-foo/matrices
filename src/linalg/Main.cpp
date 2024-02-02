#include <iostream>

#include "EquationSolverTest.h"
#include "Theorems.h"

int main() {

  // gaussian, gauss jordan, systems of linear equations
  EquationSolverTest::run();

  // fun theorem stuff
  Theorems::test();

  std::cout << "all tests pass!" << "\n";

  return 0;
}
