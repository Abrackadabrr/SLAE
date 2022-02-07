//
// Created by evgen on 06.02.2022.
//

#include "gtest/gtest.h"
#include "my_project/SlaeBaseException.hpp"
#include "my_project/matrix/ThreeDiagonalMatrix.hpp"
#include "my_project/solvers/ThreeDiagonalSolver.hpp"
#include "my_project/vectors/Vectors.h"

using std::cout;
using std::endl;


TEST(matrix, matrix_hi) {
    try {
        Slae::Matrix::ThreeDiagonalMatrix a = Slae::Matrix::ThreeDiagonalMatrix::ThreeIdentity(100000, 1, 5, 1);
//        Slae::Matrix::ThreeDiagonalMatrix a = Slae::Matrix::ThreeDiagonalMatrix::Identity(100000);
        std::vector<double> b = sequence<double>(1, 10001, 1);
        std::vector<double> res = Slae::Solvers::solveThreeDiagonal(a,b);
        cout << res;
    } catch (const Slae::SlaeBaseExceptionCpp &err) {
        auto a = err.what();
        cout << a << endl;
    }
}
