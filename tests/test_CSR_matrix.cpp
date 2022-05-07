//
// Created by evgen on 19.02.2022.
//

#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRExeption.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/vectors/Vectors.hpp"
#include "my_project/vectors/Norm.hpp"
#include "my_project/matrix/utility/Element.hpp"
#include <vector>

class CSRMATRIX_TESTS: public ::testing::Test {
protected:
    Slae::Matrix::CSR<double> m1{};
    Slae::Matrix::CSR<double> m2{};
    std::vector<double> b{};
    double tolerance = 0.00001;

    void SetUp() override {
        m1 = Slae::Matrix::CSR<double>::Diag(4, 4);
        std::vector<double> values{1, 2, 1, 3, 1, 2, 1};
        std::vector<unsigned> row{0, 2, 5, 7, 7};
        std::vector<unsigned> col{1, 3, 1, 2, 3, 1, 2};
        m2 = Slae::Matrix::CSR<double>(values, col, row, 4, 4);
        b = std::vector<double>{1, 1, 1, 1};
    }
};

TEST_F(CSRMATRIX_TESTS, EXCEPTION) {
    try {
        throw Slae::Matrix::CSRException("csrexception");
    } catch (const Slae::Matrix::CSRException& err) {
        ASSERT_STREQ("csrexception", err.what());
    }
}

TEST_F(CSRMATRIX_TESTS, ACCESS_OPERATOR) {
    ASSERT_EQ(m1(0, 0), 4) << "Operator () doesnt work properly";
    ASSERT_EQ(m2(0, 0), 0) << "Operator () doesnt work properly";
    ASSERT_EQ(m1(2, 1), 0) << "Operator () doesnt work properly";
    ASSERT_EQ(m2(1, 3), 1) << "Operator () doesnt work properly";
}

TEST_F(CSRMATRIX_TESTS, MULTIPLICATION_BY_VECTOR){
    ASSERT_NEAR(norm(m1*b - std::vector<double>{4, 4, 4, 4}, NormType::InfNorm), 0, tolerance) << "Multiplication doesnt work properly";
    ASSERT_NEAR(norm(m2*b - std::vector<double>{3, 5, 3, 0}, NormType::InfNorm), 0, tolerance) << "Multiplication doesnt work properly";
    ASSERT_NEAR(b*(m2*b) - 11, 0, tolerance) << "Multiplication doesnt work properly";
}

TEST_F(CSRMATRIX_TESTS, CONSTRUCTOR_OF_SET) {
    std::vector<double> values = {4, 4, 4, 4};
    std::vector<unsigned> i = {0, 1, 2, 3};
    std::vector<unsigned> j = {0, 1, 2, 3};
    std::set<Slae::Matrix::Element<double>> elements{};
    for(int k = 0; k < values.size(); k++) {elements.insert(Slae::Matrix::Element<double>{i[k], j[k], values[k]});}
    Slae::Matrix::CSR<double> matrix(4, 4, elements);

    ASSERT_EQ(m1.n_values(), matrix.n_values());
    ASSERT_EQ(m1.height(), matrix.height());
    ASSERT_EQ(m1.width(), matrix.width());

    for (int k = 0; k < m1.height(); k++)
        for (int (l) = 0;  (l) < m1.width(); ++(l)) {
            ASSERT_NEAR(m1(k, l), matrix(k, l), std::numeric_limits<double>::epsilon());
        }

    values = std::vector<double>{1, 2, 1, 3, 1, 2, 1};
    i = std::vector<unsigned>{0, 0, 1, 1, 1, 2, 2};
    j = std::vector<unsigned>{1, 3, 1, 2, 3, 1, 2};
    elements.clear();
    for(int k = 0; k < values.size(); k++) {elements.insert(Slae::Matrix::Element<double>{i[k], j[k], values[k]});}
    matrix = Slae::Matrix::CSR(4, 4, elements);

    ASSERT_EQ(m2.n_values(), matrix.n_values());
    ASSERT_EQ(m2.height(), matrix.height());
    ASSERT_EQ(m2.width(), matrix.width());

    for (int k = 0; k < m2.height(); k++)
        for (int (l) = 0;  (l) < m2.width(); ++(l)) {
            ASSERT_NEAR(m2(k, l), matrix(k, l), std::numeric_limits<double>::epsilon());
        }
}

TEST_F(CSRMATRIX_TESTS, TRANSPOSE) {
    auto m1T = m1.transpose();
    for (int i = 0; i < m1._width; i++)
        for (int j = 0; j < m1._height; j++)
            ASSERT_NEAR(m1(i, j), m1T(j, i), 1e-10);
    auto m2T = m2.transpose();
    for (int i = 0; i < m2._width; i++)
        for (int j = 0; j < m2._height; j++)
            ASSERT_NEAR(m2(i, j), m2T(j, i), 1e-10);
    std::vector<double> values{1, 2, 1, 3, 1, 2, 1, 1};
    std::vector<unsigned> row{0, 2, 5, 7, 7, 8};
    std::vector<unsigned> col{1, 3, 1, 2, 3, 1, 2, 2};
    auto m3 = Slae::Matrix::CSR<double>(values, col, row, 5, 4);
    auto m3T = m3.transpose();
    for (int i = 0; i < m3._width; i++)
        for (int j = 0; j < m3._height; j++)
            ASSERT_NEAR(m3(i, j), m3T(j, i), 1e-10);
}
