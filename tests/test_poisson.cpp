//
// Created by evgen on 21.03.2022.
//

#include "gtest/gtest.h"
#include "my_project/solvers/GaussSeidel.hpp"
#include "my_project/solvers/ChebyshevSimpleIteration.h"
#include "my_project/emaxples/DirihletPiosson.h"
#include "ostream"

double laplas(double x, double y) {
    return 0;
}

double border(double x, double y){
    if (std::abs(y) < 1e-5|| std::abs(y - 1)  < 1e-5 || std::abs(x)  < 1e-5 || std::abs(x - 1)  < 1e-5) return 1;
    return 0;
}

TEST(EXAMPLES, POISSON) {
    double l = 1;
    unsigned size = 30;
    double h = l/(size-1);
    auto slae = Slae::Examples::PoissonDirichletProblem(size, h, border, laplas);
    std::vector<double> initial(slae.second.size(), 0.);
    std::ofstream of("/media/evgen/Big_disc/MIPT/2nd level/Chapter 4/SLAE/cmake-build-debug/matrix.txt");
    Slae::Matrix::printToFile(of, slae.first);
    of.close();
}
