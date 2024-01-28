#include <algorithm>

#include "EquationSolverTest.h"
#include "Matrix.h"
#include "Solution.h"

typedef Matrix<double> Md;

namespace EquationSolverTest {
    void test_get_solution_type() {
    // REF of augmented matrix of inconsistent system
    Md A = {{3,2,3,4},{0,0,1,1},{0,0,0,2}};
    assert(get_solution_type(A) == SolutionType::NO_SOLUTION);

    // REF of augmented matrix of consistent system (with one soln)
    Md B = {{1,2,3,4},{0,2,0,1},{0,0,-1,2}};
    Md C = {{1,1,2,3,4},{0,2,0,1,-1},{0,0,4,-1,2},{0,0,0,-1,2},{0,0,0,0,0}};
    assert(get_solution_type(B) == SolutionType::ONE_SOLUTION);
    assert(get_solution_type(C) == SolutionType::ONE_SOLUTION);

    // REF of augmented matrix of consistent system (with infinitely many solns)
    Md D = {{5,1,2,3,4},{0,0,-1,0,1},{0,0,0,1,2}};
    Md E = {{0,1,2,3,4},{0,0,-1,0,1},{0,0,0,1,2}};
    assert(get_solution_type(D) == SolutionType::INFINITELY_MANY_SOLUTIONS);
    assert(get_solution_type(E) == SolutionType::INFINITELY_MANY_SOLUTIONS);
  }

  void test_solve_linear_system() {

    // no solution
    Md a {{1,2,1},{2,-2,2},{4,8,4}};
    Md b {{1}, {2}, {-4}};
    auto S = solve_linear_system(a, b);
    assert(S.freeVariables.empty());
    assert(S.m_solutions.empty());

    // infinite solutions
    a = {{1,1,2}, {3,3,6}};
    b = {{1}, {3}};    
    S = solve_linear_system(a, b);
    auto fun = S.compute_solution();

    assert(S.freeVariables.size() == 2);
    assert(S.toString() == "x1:1-a1-2a2, x2:a1, x3:a2 where a1, a2 are arbitrary parameters");
    assert(fun({1,2}) == std::vector<double>({-4,1,2}));

    a = {{0,0,2,4,2},{1,2,4,5,3},{-2,-4,-5,-4,3}};
    b = {{8},{-9},{6}};
    S = solve_linear_system(a, b);
    assert(S.freeVariables.size() == 2);
    assert(S.toString() == "x1:-29-2a1+3a2, x2:a1, x3:8-2a2, x4:a2, x5:-4 where a1, a2 are arbitrary parameters");
    assert(S.compute_solution()({1,2}) == std::vector<double>({-25,1,4,2,-4}));

    // only one solution
    a = {{1,0,0,0},{1,1,1,1},{1,3,9,27},{1,4,16,64}};
    b = {{10},{7},{-11},{-14}};
    S = solve_linear_system(a, b);
    assert(S.freeVariables.empty());
    assert(S.toString() == "x1: 10, x2: 2, x3: -6, x4: 1");

    // verify that for nxn matrix, we can get x = A^-1*b directly.
    auto x = (inverse(a) * b).flatten();
    std::vector<double> y;
    // resize = allocation + instance creation; reserve = allocation only
    y.resize(x.size());
    std::transform(S.m_solutions.begin(), S.m_solutions.end(), y.begin(), [](auto p){return p.first;});
    assert(x == y);

    a = {{1,1,2}, {1,-1,-1}, {1,1,-1}};
    b = {{1}, {0}, {2}};
    S = solve_linear_system(a, b);
    assert(S.freeVariables.empty());
    assert(S.toString() == "x1: 0.667, x2: 1, x3: -0.333");
  }

  void run() {
    test_get_solution_type();
    test_solve_linear_system();
  }
}
