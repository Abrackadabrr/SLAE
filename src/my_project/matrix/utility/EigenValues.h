//
// Created by evgen on 23.03.2022.
//

#ifndef MY_PROJECT_EIGENVALUES_H
#define MY_PROJECT_EIGENVALUES_H

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/vectors/Norm.hpp"
#include "my_project/vectors/Vectors.hpp"

namespace Slae::Matrix {

    template<typename T>
    T getMaxEigenValue(Slae::Matrix::CSR<T> A, T tolerance) {
        std::vector<T> u(A.height(), 10);
//        double previous_norm = norm(u);
//        u = A * u;
//        double current_norm = norm(u);
//        double current_approx = current_norm / previous_norm * previous_norm;
        double previous_approx = 0;
        auto u_next = (A * u) * (1/norm(u, NormType::SecondNorm));
        double current_approx = norm(u_next, NormType::SecondNorm) / norm(u, NormType::SecondNorm);

        while (std::abs(current_approx - previous_approx) > tolerance) {
            u = u_next;
            previous_approx = current_approx;

            u_next = A * u *(1/norm(u, NormType::SecondNorm));
            current_approx = norm(u_next, NormType::SecondNorm) / norm(u, NormType::SecondNorm);
        }
        return current_approx;
    }
}

#endif //MY_PROJECT_EIGENVALUES_H
