//
// Created by evgen on 18.04.2022.
//

#include "gtest/gtest.h"
#include "my_project/vectors/Vectors.hpp"
#include "my_project/matrix/dense/Densematrix.hpp"

Slae::Matrix::DenseMatrix<double> inverse(const Slae::Matrix::DenseMatrix<double> &a) {
    Slae::Matrix::DenseMatrix<double> inverse(3, 3, 0);
    if (a.get_row_size() == 3 && a.get_col_size() == 3) {
        double det = a(0, 0) * a(1, 1) * a(2, 2) + a(1, 0) * a(2, 1) * a(0, 2) + a(0, 1) * a(1, 2) * a(2, 0) -
                     a(0, 2) * a(1, 1) * a(2, 0) - a(0, 0) * a(1, 2) * a(2, 1) - a(2, 2) * a(0, 1) * a(1, 0);
        if (std::abs(det - 0) > 1e-10) {
            inverse(0, 0) = (1 / det) * (a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2));
            inverse(1, 1) = (1 / det) * (a(0, 0) * a(2, 2) - a(2, 0) * a(0, 2));
            inverse(2, 2) = (1 / det) * (a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1));
            inverse(0, 1) = (-1 / det) * (a(1, 0) * a(2, 2) - a(2, 0) * a(1, 2));
            inverse(0, 2) = (1 / det) * (a(1, 0) * a(2, 1) - a(2, 0) * a(1, 1));
            inverse(1, 2) = (-1 / det) * (a(0, 0) * a(2, 1) - a(2, 0) * a(0, 1));
            inverse(1, 0) = inverse(0, 1);
            inverse(2, 0) = inverse(0, 2);
            inverse(2, 1) = inverse(1, 2);
        }
    }
    return inverse;
}

TEST(MATRIX, DENSEMATRIX) {
    auto A = Slae::Matrix::DenseMatrix<double>::Diag(3,1);
    std::cout << A << '\n';
    std::cout << inverse(A) << '\n';
    std::cout << A * inverse(A) << '\n';
    auto B = Slae::Matrix::DenseMatrix<double>(4, 1, 1);
    std::vector<double> b(4, 1);
//    std::cout << (A * b) << '\n';
    std::vector<double> d = sequence(1., 13., 1.);
    Slae::Matrix::DenseMatrix<double> r(3, 3, d);
    r.write_col({1, 2, 3}, 1);
    std::cout << r;
    ASSERT_TRUE(true);
}
