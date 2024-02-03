#include <gtest/gtest.h>
#include <vector>

#include "Matrix.h"

typedef linalg::Matrix<double> Md;

using namespace linalg::Solution;

// checking types of solutions

TEST(EquationSolverTest, TypeNoSolution) {
 // REF of augmented matrix of inconsistent system
  Md A {{3,2,3,4},{0,0,1,1},{0,0,0,2}};
  EXPECT_EQ(get_solution_type(A), SolutionType::NO_SOLUTION);
}

TEST(EquationSolverTest, TypeOneSolution) {
  // REF of augmented matrix of consistent system (with one soln)
  Md A {{1,2,3,4},{0,2,0,1},{0,0,-1,2}};
  Md B {{1,1,2,3,4},{0,2,0,1,-1},{0,0,4,-1,2},{0,0,0,-1,2},{0,0,0,0,0}};
  EXPECT_EQ(get_solution_type(A), SolutionType::ONE_SOLUTION);
  EXPECT_EQ(get_solution_type(B), SolutionType::ONE_SOLUTION);
}

TEST(EquationSolverTest, TypeInfiniteSolutions) {
  // REF of augmented matrix of consistent system (with infinitely many solns)
  Md A {{5,1,2,3,4},{0,0,-1,0,1},{0,0,0,1,2}};
  Md B {{0,1,2,3,4},{0,0,-1,0,1},{0,0,0,1,2}};
  EXPECT_EQ(get_solution_type(A), SolutionType::INFINITELY_MANY_SOLUTIONS);
  EXPECT_EQ(get_solution_type(B), SolutionType::INFINITELY_MANY_SOLUTIONS);
}

// solving the linear systems

TEST(EquationSolverTest, NoSolution) {
    Md a {{1,2,1},{2,-2,2},{4,8,4}};
    Md b {{1}, {2}, {-4}};
    auto S = solve_linear_system(a, b);
    EXPECT_EQ(S.num_free_variables, 0);
    EXPECT_TRUE(S.m_solutions.empty());
}


TEST(EquationSolverTest, InfiniteSolutions_System1) {
    Md a {{1,1,2}, {3,3,6}};
    Md b {{1}, {3}};    
    auto S = solve_linear_system(a, b);
    auto fun = S.get_compute_function();

    EXPECT_EQ(S.num_free_variables, 2);
    EXPECT_EQ(S.toString(), "x1:1-a1-2a2, x2:a1, x3:a2 where a1, a2 are arbitrary parameters");
    EXPECT_TRUE(fun({1,2}) == std::vector<double>({-4,1,2}));
}

TEST(EquationSolverTest, InfiniteSolutions_System2) {
    Md a {{0,0,2,4,2},{1,2,4,5,3},{-2,-4,-5,-4,3}};
    Md b {{8},{-9},{6}};
    auto S = solve_linear_system(a, b);
    auto fun = S.get_compute_function();

    EXPECT_EQ(S.num_free_variables, 2);
    EXPECT_EQ(S.toString(), "x1:-29-2a1+3a2, x2:a1, x3:8-2a2, x4:a2, x5:-4 where a1, a2 are arbitrary parameters");
    EXPECT_TRUE(fun({1,2}) == std::vector<double>({-25,1,4,2,-4}));
}

TEST(EquationSolverTest, OneSolution_System1) {
    Md a {{1,0,0,0},{1,1,1,1},{1,3,9,27},{1,4,16,64}};
    Md b {{10},{7},{-11},{-14}};
    auto S = solve_linear_system(a, b);
    // auto fun = S.get_compute_function();

    EXPECT_EQ(S.num_free_variables, 0);
    EXPECT_EQ(S.toString(), "x1: 10, x2: 2, x3: -6, x4: 1");
}

// analytical solution using inverse of an invertible n x n matrix; should be same result as above test
TEST(EquationSolverTest, OneSolution_System1_Inverse) {
    Md a {{1,0,0,0},{1,1,1,1},{1,3,9,27},{1,4,16,64}};
    Md b {{10},{7},{-11},{-14}};
    auto S = solve_linear_system(a, b);

    // verify that for nxn matrix, we can get x = A^-1*b directly.
    std::vector<double> x = (inverse(a) * b).flatten();
    std::vector<double> y;
    // resize = allocation + instance creation; reserve = allocation only
    y.resize(x.size());
    std::transform(S.m_solutions.begin(), S.m_solutions.end(), y.begin(), [](auto p){return p.val;});
    EXPECT_TRUE(x == y);
}


TEST(EquationSolverTest, OneSolution_System2) {
    Md a {{1,1,2}, {1,-1,-1}, {1,1,-1}};
    Md b {{1}, {0}, {2}};
    auto S = solve_linear_system(a, b);
    // auto fun = S.get_compute_function();

    EXPECT_EQ(S.num_free_variables, 0);
    EXPECT_EQ(S.toString(), "x1: 0.667, x2: 1, x3: -0.333");
}
