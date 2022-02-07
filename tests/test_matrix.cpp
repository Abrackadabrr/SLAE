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

/*
 * Link to solution of this system
 * https://matrixcalc.org/slu.html#solve-using-Gaussian-elimination(%7B%7B3,1,0,0,1%7D,%7B1,3,1,0,1%7D,%7B0,1,3,1,1%7D,%7B0,0,1,3,2%7D%7D)
 */
TEST(threediagonalmatrix, matrix_passed)
{
    try {
        int n = 4;
        Slae::Matrix::ThreeDiagonalMatrix a = Slae::Matrix::ThreeDiagonalMatrix::ThreeDiagonal(n, 1, 3, 1);
        std::vector<double> b{1, 1, 1, 2};

        std::vector<double> sol{14, 13, 2, 36};
        sol = sol*(1./55);
        std::vector<double> res = Slae::Solvers::solveThreeDiagonal(a,b);
        for (int i = 0; i < n-1; ++i) {
            ASSERT_NEAR(res[i], sol[i], 10e-10);
        }
        cout << res << endl;
    }
    catch (const Slae::SlaeBaseExceptionCpp &err) {
        cout << err.what() << endl;
    }
}

TEST(threediagonalmatrix, matrix_failure)
{
    try {
        int n = 4;
        Slae::Matrix::ThreeDiagonalMatrix a = Slae::Matrix::ThreeDiagonalMatrix::ThreeDiagonal(n, 1, 3, 1);
        std::vector<double> b{1, 1, 1, 2};

        std::vector<double> sol{12, 13, 2, 36};  // solution from lin above
        sol = sol*(1./55);
        std::vector<double> res = Slae::Solvers::solveThreeDiagonal(a,b);
        for (int i = 0; i < n-1; ++i) {
            ASSERT_NEAR(res[i], sol[i], 10e-10);
        }
        cout << res << endl;
    }
    catch (const Slae::SlaeBaseExceptionCpp &err) {
        cout << err.what() << endl;
    }
}
