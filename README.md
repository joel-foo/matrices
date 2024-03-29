### Matrices in C++

- Fundamental operations:

  - Matrix addition
  - Matrix subtraction
  - Matrix and scalar multiplication
  - Transpose
  - Concatenating two matrices either horizontally or vertically

- Other operations:

  - Gaussian elimination
  - Gauss jordan elimination
  - Rank of matrix
  - Row space of matrix
  - Column space of matrix
  - **Solving linear systems (any kind!)**

- This project uses **CMake** and **Google Tests**.

### Build using make

```bash
# in the root project directory
cmake -S . -B build

cd build

make

# to run tests
make test / ctest
```

### Example

_After building the project using make, the following example program can found in build/src directory. Navigate to that directory and run ./Demo to execute the program._

Consider a system of equations in 5 variables: x1 - x5, which is represented as a coefficient matrix A and column vector of constants b. Note that the types of these matrices must be parametrised with `Double`, in the cases where solving equations are involved.

An example program of solving a system with an infinite number of solutions:

```cpp
#include <iostream>

#include "Matrix.h"

// make your own aliases
typedef linalg::Matrix<double> Md;

int main() {

  using std::cout;

  // A is the coefficient matrix
  Md A {{0,2,2,1,-2},{0,0,1,1,1},{0,0,0,0,2}};

  // b is the column matrix of constants
  Md b {{2},{3},{4}};

  linalg::Solution::SystemSolution solution = linalg::solve_linear_system(A, b);

  // prints "x1:a1, x2:2+0.5a2, x3:1-a2, x4:a2, x5:2 where a1, a2 are arbitrary parameters"
  cout << solution << "\n";

  // prints 2, which is the number of free variables (a1, a2)
  cout << solution.num_free_variables << "\n";

  // fun is a function that takes in a list of initialiser values (each corresponding to the free variables in their natural order (a1, a2, ..., an) and returns the values of all the variables for that particular combination of values
  auto fun = solution.get_compute_function();

  // setting a1 = 1, a2 = 2
  // prints: (1, 3, -1, 2, 2), which corresponds to x1: 1, x2: 3, x3: -1, x4: 2, x5: 2
  cout << fun({1, 2}) << "\n";

  // prints: (2, 3.5, -2, 3, 2)
  cout << fun({2, 3}) << "\n";

  return 0;
}
```
