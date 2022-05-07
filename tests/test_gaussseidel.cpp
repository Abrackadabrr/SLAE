//
// Created by evgen on 23.03.2022.
//

#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/solvers/GaussSeidel.hpp"
#include "my_project/matrix/ThreeDiagonal/ThreeDiagonalMatrix.hpp"


class GAUSSSEIDEL_TESTS : public ::testing::Test {
protected:
    Slae::Matrix::CSR<double> a1{};
    Slae::Matrix::CSR<double> a2{};
    std::vector<double> b1{};
    std::vector<double> b2{};
    std::vector<double> answer1{};
    std::vector<double> answer2{};
    double tolerance = 1e-5;

    void SetUp() override {
        a1 = Slae::Matrix::CSR<double>::Diag(1, 100);
        a2 = Slae::Matrix::ThreeDiagonalMatrix::ThreeDiagonal(100, 1, 2, 1).toCSR();
        b1 = std::vector<double>(100, 1.101);
        b2 = std::vector<double>(100, 0);
        answer1 = std::vector<double>(100, 1.101);
        answer2 = std::vector<double>(100,0);
    }
};

TEST_F(GAUSSSEIDEL_TESTS, GAUSSSEIDEL){
    std::vector<double> res = Slae::Solvers::GaussSeidel(a1, b1, std::vector<double>(100, 1), tolerance/10);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);
    res = Slae::Solvers::GaussSeidel(a2, b2, std::vector<double>(100, 100), tolerance*1e-4);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};