//
// Created by evgen on 03.05.2022.
//

#include "gtest/gtest.h"
#include "my_project/solvers/Gauss.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/matrix/dense/Densematrix.hpp"
#include "my_project/matrix/ThreeDiagonal/ThreeDiagonalMatrix.hpp"
#include "my_project/vectors/Norm.hpp"
#include "my_project/solvers/ConjugateGradient.h"
#include "my_project/solvers/GMRES.h"
#include "my_project/emaxples/DirihletPiosson.h"

double laplas(double x, double y) {
    return -std::sin(x) - std::cos(y);
}

double border(double x, double y) {
    return std::sin(x) + std::cos(y);
}

TEST(GMRES, GaussDownShift) {
    double tolerance = 1e-10;
    auto A = Slae::Matrix::CSR<double>::Diag(1, 10);
    auto denseA = Slae::Matrix::DenseMatrix<double>::Diag(10, 3);
    std::vector<double> b(10, 1);
    ASSERT_NEAR(norm(GaussDownshift(A, b) - b, NormType::SecondNorm), 0, tolerance);
    ASSERT_NEAR(norm(GaussDownshift(denseA, b) - (b * (1. / 3)), NormType::SecondNorm), 0, tolerance);

    std::vector<double> data = {1, 2, 0, 3};
    b = {1, 3};
    ASSERT_NEAR(norm(GaussDownshift(Slae::Matrix::DenseMatrix<double>(2, 2, data), b) - std::vector<double>{-1, 1},
                     NormType::SecondNorm), 0, tolerance);
}

TEST(GMRES, Gmres) {
    double tolerance = 1e-12;
    int n = 71;
    auto slae = Slae::Examples::PoissonDirichletProblem(n, 1./(n-1), border, laplas);
    auto A = slae.first;
    auto b = slae.second;
    auto resGMRES = Slae::Solvers::GMRESM(n, A, b, std::vector<double>((n-2)*(n-2), 0), tolerance);
    std::cout << "gmres ends" << std::endl;
    auto resGS = Slae::Solvers::ConjugateGradient(A, b, std::vector<double>((n-2)*(n-2), 0), tolerance);
    ASSERT_NEAR(norm(resGMRES - resGS), 0, 1e-9);
}
