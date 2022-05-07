//
// Created by evgen on 22.03.2022.
//

#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "set"
#include "my_project/matrix/utility/Element.hpp"
#include "my_project/vectors/Vectors.hpp"
#include "my_project/solvers/ChebyshevSimpleIteration.h"

Slae::Matrix::CSR<double> give_matrix(unsigned size){
    std::set<Slae::Matrix::Element<double>> tree;
    std::vector<double> values = sequence(-1., 1. + 2./size, 2./size);
    std::cout << values <<'\n';
    unsigned ind = 0;
    for (auto value: values) {
        if (std::abs(value) > 1e-3) {tree.insert({ind, ind, value}); ind++;}
    }
    return Slae::Matrix::CSR<double>(size, size, tree);
}

TEST(A, D) {
    unsigned size = 100;
    auto A = give_matrix(size);
    std::cout << A << '\n';
    std::vector<double> b(size, 0);
    std::vector<double> initial(size, 10);
    auto res = Slae::Solvers::FastSimpleIterationWithSteps<double, 4>(A, b, initial, {-1, 1}, 1e-5);
    std::cout << res.first << '\n';
}
