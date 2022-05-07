//
// Created by evgen on 12.03.2022.
//

#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/matrix/ThreeDiagonal/ThreeDiagonalMatrix.hpp"

TEST(THREEDIAGONAL, CSR_CAST) {
    auto a = Slae::Matrix::ThreeDiagonalMatrix::ThreeDiagonal(5, 1, 2, 3);
    auto b = a.toCSR();
    for (int i = 0; i < a.rows(); i++) {
        if (i != 0 && i != a.rows() -1) {
            ASSERT_NEAR(b(i, i - 1), 1, 1e-5);
            ASSERT_NEAR(b(i, i), 2, 1e-5);
            ASSERT_NEAR(b(i, i + 1), 3, 1e-5);
        } else if (i == 0) {
            ASSERT_NEAR(b(i, i), 2, 1e-5);
            ASSERT_NEAR(b(i, i + 1), 3, 1e-5);
        } else if (i == a.rows() - 1) {
            ASSERT_NEAR(b(i, i - 1), 1, 1e-5);
            ASSERT_NEAR(b(i, i), 2, 1e-5);
        }
    }
}
