//
// Created by evgen on 23.03.2022.
//

#ifndef MY_PROJECT_CONJUGATEGRADIENT_H
#define MY_PROJECT_CONJUGATEGRADIENT_H

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "vector"
#include "my_project/vectors/Norm.hpp"

namespace Slae::Solvers {
    std::pair<std::vector<double>, std::vector<std::vector<double>>>
    ConjugateGradientWithSteps(const Slae::Matrix::CSR<double> &A, const std::vector<double> &b,
                      const std::vector<double> &initial, double tolerance) {
        std::vector<double> x = initial;
        std::vector<double> r = A * initial - b;
        std::vector<double> d = r;

        std::vector<std::vector<double>> steps{};
        steps.reserve(1000);

        double alpha = 0;
        double beta = 0;
        double r_squared = 0;

        while (norm(r, NormType::InfNorm) > tolerance) {
            r_squared = (r * r);
            alpha = r_squared/ (d * (A * d));

            x = x - alpha * d; // GHANGE
            r = A * x - b;

            beta = (r * r)/r_squared;
            d = r + beta*d;

            steps.push_back(x);
        }

        return {x, steps};
    }

    std::vector<double>
    ConjugateGradient(const Slae::Matrix::CSR<double> &A, const std::vector<double> &b,
                               const std::vector<double> &initial, double tolerance) {
        std::vector<double> x = initial;
        std::vector<double> r = A * initial - b;
        std::vector<double> d = r;

        double alpha = 0;
        double beta = 0;
        double r_squared = 0;

        while (norm(r, NormType::InfNorm) > tolerance) {
            r_squared = (r * r);
            alpha = r_squared/ (d * (A * d));

            x = x - alpha * d;
            r = A * x - b;

            beta = (r * r)/r_squared;
            d = r + beta*d;
        }

        return x;
    }
}


#endif //MY_PROJECT_CONJUGATEGRADIENT_H
