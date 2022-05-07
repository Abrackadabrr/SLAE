//
// Created by evgen on 19.03.2022.
//

#ifndef MY_PROJECT_STEEPESTGRADIENDDESCENT_H
#define MY_PROJECT_STEEPESTGRADIENDDESCENT_H

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/vectors/Vectors.hpp"
#include "my_project/vectors/Norm.hpp"
#include "iostream"

namespace Slae::Solvers {
    template<typename T>
    std::vector <T>
    SteepestGradientDescent(const Slae::Matrix::CSR<T> &A, const std::vector <T> &b, const std::vector <T> initial,
                    const T &tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in Jacobi isn't squared");

        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        std::vector <T> x = initial;
        std::vector<double> r = A*x - b;
        std::vector<double> Ar{};
        T znam = 0;
        T tao = 0;
        while (norm(r) > tolerance) {
            Ar = A * r;
            znam = (r*Ar);
            tao = (r*r) / znam;
            x = x - tao * r; // CHANGE
            r = r - tao * Ar;
        }
        return x;
    }


    template<typename T>
    std::pair<std::vector<T>, std::vector<std::vector<T>>>
    SteepestGradientDescentWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector <T> &b, const std::vector <T> initial,
                            const T &tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in Jacobi isn't squared");

        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }
        std::vector<std::vector<double>> steps{}; steps.reserve(1000);
        std::vector <T> x = initial;
        std::vector<double> r = A*x - b;
        std::vector<double> Ar{};
        T znam = 0;
        T tao = 0;
        while (norm(r) > tolerance) {
            Ar = A * r;
//            std::cout << norm(r) << '\n';
            znam = (r*Ar);
            tao = (r*r) / znam;
//            std::cout << tao << '\n';
            x = x - tao * r;
            r = r - tao * Ar;
            steps.push_back(x);
        }
        return {x, steps};
    }
}
#endif //MY_PROJECT_STEEPESTGRADIENDDESCENT_H
