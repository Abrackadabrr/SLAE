//
// Created by evgen on 23.03.2022.
//


#include "gtest/gtest.h"
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/matrix/utility/EigenValues.h"
#include "my_project/emaxples/DirihletPiosson.h"

double laplas(double x, double y) {
    return 0;
}

double border(double x, double y) {
    return 10;
}

TEST(A, S) {
    double tolerance = 1e-11;  // calculation tolerance
    unsigned size = 20;  // amount of nodes per side
    double h = 1./size;

    auto slae = Slae::Examples::PoissonDirichletProblem(size, h, border, laplas);
    auto A = (-1) * slae.first * size * size;

    size = size - 2;  // now size*size is amount of rows in slae.first
    double l_max = (8/(h*h))*pow((std::sin((size - 1)*M_PI/(2*size))), 2);

    auto l_max_alg = Slae::Matrix::getMaxEigenValue(A, tolerance);
    std::cout << l_max  << '\n' << l_max_alg;
}
