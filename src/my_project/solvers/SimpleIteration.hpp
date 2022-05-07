//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_SIMPLEITERATION_HPP
#define SLAE_SIMPLEITERATION_HPP
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/SlaeBaseException.hpp"
#include "my_project/vectors/Norm.hpp"

namespace Slae::Solvers {
    template<typename T>
    std::vector<T>
    SimpleIteration(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> initial,
                    const T &tolerance, const T &tao) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in Jacobi isn't squared");

        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        std::vector<T> x = initial;
        auto r = A * x - b;
        while (norm(r) > tolerance) {
            x = x - tao*r;
            r = A * x - b;
        }
        return x;
    }

    template<typename T>
    std::pair<std::vector<T>, std::vector<std::vector<T>>>
    SimpleIterationWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> initial,
                    const T &tolerance, const T &tao) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in SimpleIteration isn't squared");

        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        std::vector<std::vector<T>> steps{}; steps.reserve(100);

        std::vector<T> x = initial;
        std::vector<T> r = A * x - b;
        while (norm(r) > tolerance) {
            x = x - tao*r;
            r = A * x - b;
            steps.push_back(x);
        }
        return {x, steps};
    }
}
#endif//SLAE_SIMPLEITERATION_HPP
