#include <cassert>

#include "Matrix.h"

typedef linalg::Matrix<double> Md;

namespace Theorems {
  namespace RankTheorems {
    // A: 3 x 4 (m x n) 
    Md A {{1,1,0,-2},{2,0,2,2},{4,1,3,1}};
    // B: 4 x 3 (n x k)
    Md B {{1,2,0},{2,3,8},{3,3,6},{7,2,5}};
    Md AB = A*B;
    int rank_A = get_rank(A), rank_B = get_rank(B), rank_AB = get_rank(A*B);
    void test() {
      assert(rank_A == 2);

      // rank(A) <= min(m, n)
      assert(rank_A <= std::min(A.getRows(), A.getCols()));

      // only zero matrix has rnak 0 
      Md zeroMatrix (2,3);
      assert(get_rank(zeroMatrix) == 0);

      // rank(AB) <= min(rank(A), rank(B))
      assert(rank_AB <= std::min(rank_A, rank_B));

      // slyvester's rank inequality: rank(A) + rank(B) - n <= rank(AB)
      assert(rank_A + rank_B - A.getCols() <= rank_AB);

      // subaddivity: rank(A + C) <= rank(A) + rank(C), both A and C of same dimension. 
      Md C = {{3,4,2,5},{9,10,2,3},{3,3,11,0}};
      assert(get_rank(A + C) <= rank_A + get_rank(C));

      // rank(A^T) = rank(A) = rank(AA^T) = rank(A^TA)
      assert(rank_A == get_rank(A.transpose()));
      assert(rank_A == get_rank(A * A.transpose()));
      assert(rank_A == get_rank(A.transpose() * A));
    }
  }
  void verify() {
    RankTheorems::test();
  }
}
