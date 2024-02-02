#include <cassert>

#include "Matrix.h"

typedef linalg::Matrix<int> Mint;
typedef linalg::Matrix<double> Md;

namespace MatrixTest {
    void test_add_matrices() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      Mint expected {{2,4,6},{8,10,12},{14,16,18}};
      assert(A + A == expected);
    }

    void test_subtract_matrices() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      Mint expected (3,3,0);
      assert(A - A == expected);
    }

    void test_multiply_scalar() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      Mint expected {{2,4,6},{8,10,12},{14,16,18}};
      assert(A * 2 == expected);
    }

    void test_multiply_matrices() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      Mint B {{9,8,7},{6,5,4},{3,2,1}};
      Mint expected {{30,24,18},{84,69,54},{138,114,90}};
      assert(A * B == expected);
    }

    void test_transpose_matrix() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      Mint expected {{1,4,7},{2,5,8},{3,6,9}};
      assert(A.transpose() == expected);
    }

    void test_get_matrix() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      assert(A(2,2) == 9);

      bool exceptionThrown = false;
      try {
        A(3,3);
      } catch (const std::runtime_error& error) {
        exceptionThrown  = true;
      }
      assert(exceptionThrown);
    }

    void test_set_matrix() {
      Mint A {{1,2,3},{4,5,6},{7,8,9}};
      A(2,2) = 10;
      assert(A(2,2) == 10);

      bool exceptionThrown = false;
      try {
        A(3,3) = 10;
      } catch (const std::runtime_error& error) {
        exceptionThrown = true;
      }
      assert(exceptionThrown);
    }

    void test_gaussian_elimination() {
      Md A {{2,1,-1,8},{-3,-1,2,-11},{-2,1,2,-3}};
      Md expected {{-3,-1,2,-11}, {0,5/3.0,2/3.0,13/3.0}, {0,0,0.2,-0.2}};
      assert(gaussian_elimination(A) == expected);
    }

    void test_gauss_jordan_elimination() {
      Md A {{2,1,-1,8},{-3,-1,2,-11},{-2,1,2,-3}};
      Md expected {{1,0,0,2},{0,1,0,3},{0,0,1,-1}};
      assert(gauss_jordan_elimination(A) == expected);
    }

    void test_get_rank() {
      Md A {{1,2,1},{-2,-3,1},{3,5,0}};
      assert(get_rank(A) == 2);
    }

    void test_get_row_space(){
      Md A {{1,2,1},{-2,-3,1},{3,5,0}};
      TwoDVector<double> expected {{1,0,-5}, {0,1,3}};
      assert(get_row_space(A) == expected);
    }

    void test_get_col_space() {
      Md A {{1,2,1},{-2,-3,1},{3,5,0}};
      TwoDVector<double> expected {{1,0,0}, {0,1,0}};
      assert(get_col_space(A) == expected);
    }

    void test_inverse() {
      // nxn
      Md A {{2,-1,0}, {-1,2,-1}, {0,-1,2}};
      Md expected {{0.75,0.5,0.25}, {0.5,1,0.5}, {0.25,0.5,0.75}};
      assert(inverse(A) == expected);

      // not n x n
      bool exceptionThrown = false;
      A = {{2,-1,0},{-1,2,-1}};
      try {
        inverse(A);
      } catch (const std::runtime_error& error) {
        exceptionThrown = true;
      }
      assert(exceptionThrown);
    }

    void run() {
      test_add_matrices();
      test_subtract_matrices();
      test_multiply_scalar();
      test_multiply_matrices();
      test_transpose_matrix();
      test_get_matrix();
      test_set_matrix();
      test_gaussian_elimination();
      test_gauss_jordan_elimination();
      test_get_rank();
      test_get_row_space();
      test_get_col_space();
      test_inverse();
  }
}
