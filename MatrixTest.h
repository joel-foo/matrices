#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

namespace MatrixTest {

    void test_add_matrices();
    void test_subtract_matrices();
    void test_multiply_scalar();
    void test_multiply_matrices();
    void test_transpose_matrix();
    void test_get_matrix();
    void test_set_matrix();

    void test_gaussian_elimination();
    void test_gauss_jordan_elimination();
    void test_get_rank();
    void test_get_row_space();
    void test_get_col_space();
    void test_inverse();

    void run();
}

#endif
