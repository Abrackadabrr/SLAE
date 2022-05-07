//
// Created by evgen on 19.02.2022.
//

#ifndef SLAE_JACOBI_HPP
#define SLAE_JACOBI_HPP

#include "my_project/SlaeBaseException.hpp"
#include <sstream>
#include "my_project/vectors/Norm.hpp"

using std::vector;

namespace Slae::Solvers {
    template<typename T>
    vector<T>
    Jacobi(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial, T tolerance) {
        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in Jacobi isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;
        vector<T> x_next(b.size());

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            for (unsigned i = 0; i < b.size(); ++i) {
                auto sum = static_cast<T>(0);
                for (int j = A._row_index[i]; j < A._row_index[i + 1]; ++j)
                    if (i != A._col[j]) sum += A._value[j] * x[A._col[j]];
                x_next[i] = ((b[i] - sum) / A(i, i));
            }
            x = x_next;
            r = (A * x) - b;
        }
        return x;
    }

    template<typename T>
    std::pair<vector<T>, vector<vector<T>>>
    JacobiWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial,
                          T tolerance) {
        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in Jacobi isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;
        vector<T> x_next(b.size());
        vector<vector<T>> steps{};
        steps.reserve(1000);

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            for (unsigned i = 0; i < b.size(); ++i) {
                auto sum = static_cast<T>(0);
                for (int j = A._row_index[i]; j < A._row_index[i + 1]; ++j)
                    if (i != A._col[j]) sum += A._value[j] * x[A._col[j]];
                x_next[i] = ((b[i] - sum) / A(i, i));
            }
            x = x_next;
            r = (A * x) - b;

            // add an error
            steps.push_back(x);
        }
        return {x, steps};
    }
}


#endif//SLAE_JACOBI_HPP
