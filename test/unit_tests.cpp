#include <gtest/gtest.h>
#include <vector>

#include "Matrix.h"

typedef linalg::Matrix<int> Mint;
typedef linalg::Matrix<double> Md;

// Test fundamental matrix operations

TEST(MatrixTest, AddMatrices) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{2,4,6},{8,10,12},{14,16,18}};
  EXPECT_EQ(A + A, expected);
}

TEST(MatrixTest, SubtractMatricesInt) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected (3, 3, 0);
  EXPECT_EQ(A - A, expected);
}

TEST(MatrixTest, SubtractMatricesDouble) {
  Md A {{1,2,3},{4,5,6},{7,8,9}};
  Md B {{1.1,1.2,1.3},{2.1,2.2,2.3},{3.1,3.2,3.3}};
  Md expected {{-0.1,0.8,1.7},{1.9,2.8,3.7},{3.9,4.8,5.7}};
  EXPECT_EQ(A - B, expected);
}

TEST(MatrixTest, MultiplyMatrixByScalar) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{3,6,9},{12,15,18},{21,24,27}};
  EXPECT_EQ(A * 3, expected);
}

TEST(MatrixTest, MultiplyMatrices) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  Mint B {{9,8,7},{6,5,4},{3,2,1}};
  Mint expected {{30,24,18},{84,69,54},{138,114,90}};
  EXPECT_EQ(A * B, expected);
}

TEST(MatrixTest, TransposeMatrix) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  Mint expected {{1,4,7},{2,5,8},{3,6,9}};
  EXPECT_EQ(A.transpose(), expected);
}


TEST(MatrixTest, GetElementFromMatrix) {
  Md A {{1,2,3},{4,5,6},{7,8,9}};
  EXPECT_EQ(A(2, 2), 9);

  EXPECT_THROW(A(3, 3), std::runtime_error);
}

TEST(MatrixTest, SetElementInMatrix) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  A(2, 2) = 10;
  EXPECT_EQ(A(2, 2), 10);

  EXPECT_THROW(A(3, 3) = 10, std::runtime_error);
}

TEST(MatrixTest, FlattenMatrix) {
  Mint A {{1,2,3},{4,5,6},{7,8,9}};
  std::vector<int> expected({1,2,3,4,5,6,7,8,9});
  EXPECT_EQ(A.flatten(), expected);
}

// Test more exciting operations

TEST(MatrixTest, GaussianElimination) {
  Md A {{2,1,-1,8},{-3,-1,2,-11},{-2,1,2,-3}};
  Md expected {{-3,-1,2,-11}, {0,5/3.0,2/3.0,13/3.0}, {0,0,0.2,-0.2}};
  EXPECT_EQ(linalg::gaussian_elimination(A), expected);
}

TEST(MatrixTest, GaussJordanElimination) {
  Md A {{2,1,-1,8},{-3,-1,2,-11},{-2,1,2,-3}};
  Md expected {{1,0,0,2},{0,1,0,3},{0,0,1,-1}};
  EXPECT_EQ(linalg::gauss_jordan_elimination(A), expected);
}

TEST(MatrixTest, Rank) {
  Md A {{1,2,1},{-2,-3,1},{3,5,0}};
  EXPECT_EQ(linalg::get_rank(A), 2);
}

TEST(MatrixTest, RowSpace){
  Md A {{1,2,1},{-2,-3,1},{3,5,0}};
  TwoDVector<double> expected {{1,0,-5}, {0,1,3}};
  EXPECT_TRUE(linalg::get_row_space(A) == expected);
}

TEST(MatrixTest, ColSpace){
  Md A {{1,2,1},{-2,-3,1},{3,5,0}};
  TwoDVector<double> expected {{1,0,0}, {0,1,0}};
  EXPECT_TRUE(linalg::get_col_space(A) == expected);
}

TEST(MatrixTest, Inverse){
  Md A {{2,-1,0}, {-1,2,-1}, {0,-1,2}};
  Md expected {{0.75,0.5,0.25}, {0.5,1,0.5}, {0.25,0.5,0.75}};
  EXPECT_EQ(linalg::inverse(A), expected);

  A = {{2,-1,0},{-1,2,-1}};
  EXPECT_THROW(linalg::inverse(A), std::runtime_error);
}
