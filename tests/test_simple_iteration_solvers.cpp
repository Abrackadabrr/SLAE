//
// Created by evgen on 20.02.2022.
//
#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/solvers/Jacobi.hpp"
#include "my_project/solvers/GaussSeidel.hpp"
#include "my_project/solvers/SimpleIteration.hpp"
#include "my_project/matrix/ThreeDiagonal/ThreeDiagonalMatrix.hpp"
#include "my_project/solvers/SteepestGradiendDescent.h"
#include "my_project/solvers/ChebyshevSimpleIteration.h"

class ITERATION_METHODS_TESTS : public ::testing::Test {
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

TEST_F(ITERATION_METHODS_TESTS, JACOBI){
    std::vector<double> res = Slae::Solvers::Jacobi(a1, b1, std::vector<double>(100, 1), tolerance/10);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);
    res = Slae::Solvers::GaussSeidel(a2, b2, std::vector<double>(100, 100), tolerance*1e-4);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};

TEST_F(ITERATION_METHODS_TESTS, GAUSSSEIDEL){
    std::vector<double> res = Slae::Solvers::GaussSeidel(a1, b1, std::vector<double>(100, 1), tolerance/10);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);
    res = Slae::Solvers::GaussSeidel(a2, b2, std::vector<double>(100, 100), tolerance*1e-4);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};

TEST_F(ITERATION_METHODS_TESTS, SIMPLEITERMETHOD){
    std::vector<double> res = Slae::Solvers::SimpleIteration(a1, b1, std::vector<double>(100, 1), tolerance/10, 0.3);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);

    res = Slae::Solvers::SimpleIteration(a2, b2, std::vector<double>(100, 100), tolerance*1e-4, 0.3);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};

TEST_F(ITERATION_METHODS_TESTS, FAST_GAUSSSEIDEL){
    std::vector<double> res = Slae::Solvers::FastGaussSeidel(a1, b1, std::vector<double>(100, 1), tolerance/10, 0.3);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);

    res = Slae::Solvers::FastGaussSeidel(Slae::Matrix::ThreeDiagonalMatrix::ThreeDiagonal(5, 1, 2, 1).toCSR(), vector<double>(5 ,0), std::vector<double>(5, 4), tolerance*1e-4, 0.75);
//    res = Slae::Solvers::FastGaussSeidel(a2, b2, std::vector<double>(100, 100), tolerance*1e-4, 0.75);
//    EXPECT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
    std::cout << res << '\n';
};


TEST_F(ITERATION_METHODS_TESTS, SGD){
    std::vector<double> res = Slae::Solvers::SteepestGradientDescent(a1, b1, std::vector<double>(100, 1), tolerance/10);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);

    res = Slae::Solvers::SteepestGradientDescent(a2, b2, std::vector<double>(100, 100), tolerance*1e-4);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};


TEST_F(ITERATION_METHODS_TESTS, FASTSIMPLEITERMETHOD){
    std::vector<double> res = Slae::Solvers::FastSimpleIteration<double, 10>(a1, b1, std::vector<double>(100, 1.), {0.1, 1.}, tolerance/10);
    ASSERT_NEAR(norm(answer1 - res, NormType::InfNorm), 0, tolerance);

    res = Slae::Solvers::FastSimpleIteration<double, 10>(a2, b2, std::vector<double>(100, 100), {0.1, 10}, tolerance*1e-4);
    ASSERT_NEAR(norm(answer2 - res, NormType::InfNorm), 0, tolerance);
};