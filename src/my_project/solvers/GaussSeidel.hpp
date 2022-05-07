//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_GAUSSSEIDEL_HPP
#define SLAE_GAUSSSEIDEL_HPP

#include "vector"

using std::vector;

#include "my_project/matrix/sparse/CSRmatrix.h"

namespace Slae::Solvers {
    template<typename T>
    vector<T> GaussSeidelIteration(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, std::vector<T> x) {
        for (unsigned i = 0; i < b.size(); ++i) {
            auto sum = static_cast<T>(0);
            for (int j = A._row_index[i]; j < A._row_index[i + 1]; ++j)
                if (i != A._col[j]) sum += A._value[j] * x[A._col[j]];
            x[i] = ((b[i] - sum) / A(i, i));
        }
        return x;
    }

    template<typename T>
    vector<T> GaussSeidelIterationUpsideDown(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, std::vector<T> x) {
        for (int i = b.size() - 1; i >= 0; i--) {
            auto sum = static_cast<T>(0);
            for (int j = A._row_index[i]; j < A._row_index[i + 1]; ++j)
                if (i != A._col[j]) sum += A._value[j] * x[A._col[j]];
            x[i] = ((b[i] - sum) / A(i, i));
        }
        return x;
    }

    template<typename T>
    vector<T>
    GaussSeidel(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial, T tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in GaussSeidel isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) { // DOUBLE NO COMPARE WITH 0
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            x = Slae::Solvers::GaussSeidelIteration(A, b, x);
            r = (A * x) - b;
        }
        return x;
    }

    template<typename T>
    vector<T>
    SymmetricGaussSeidel(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial,
                         T tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in GaussSeidel isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            x = Slae::Solvers::GaussSeidelIteration(A, b, x);
            x = Slae::Solvers::GaussSeidelIterationUpsideDown(A, b, x);
            r = (A * x) - b;
        }
        return x;
    }

    template<typename T>
    vector<T>
    FastGaussSeidel(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial, T tolerance,
                    T rho) {
        vector<T> previous = initial;
        vector<T> current = Slae::Solvers::GaussSeidelIteration(A, b, initial);
        vector<T> next = current;

        vector<T> r = A * next - b;
        auto mu0 = static_cast<T>(1.);
        auto mu1 = static_cast<T>(1. / rho);
        auto mu = (static_cast<T>(2.) / rho) * mu1 - mu0;

        while (norm(r) > tolerance) {
            next = ((static_cast<T>(2.) * mu1) / (rho * mu)) * Slae::Solvers::GaussSeidelIteration(A, b, current) -
                   (mu0 / mu) * previous;
            r = (A * next) - b;

            // refresh all states
            previous = current;
            current = next;
            mu0 = mu1;
            mu1 = mu;
            mu = (static_cast<T>(2.) / rho) * mu1 - mu0;
        }
        return next;
    }


    template<typename T>
    std::pair<vector<T>, vector<vector<T>>>
    GaussSeidelWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial,
                         T tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in GaussSeidel isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;
        vector<vector<T>> steps{};
        steps.reserve(100);

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            x = Slae::Solvers::GaussSeidelIterationUpsideDown(A, b, x);
            r = (A * x) - b;

            // add an error
            steps.push_back(x);
        }
        return {x, steps};
    }

    template<typename T>
    std::pair<vector<T>, vector<vector<T>>>
    SymmetricGaussSeidelWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial,
                                  T tolerance) {

        if (A._height != A._width) throw Slae::SlaeBaseExceptionCpp("Matrix in GaussSeidel isn't squared");
        for (unsigned i = 0; i < A._height; i++)
            if (A(i, i) == 0) {
                std::stringstream buff;
                buff << "Zero on position (" << i << ", " << i << ")";
                throw Slae::SlaeBaseExceptionCpp(buff.str());
            }

        vector<T> r(b.size());
        vector<T> x = initial;
        vector<vector<T>> steps{};
        steps.reserve(1000);
        steps.push_back(initial);

        r = (A * x) - b;
        while (norm(r) > tolerance) {
            x = Slae::Solvers::GaussSeidelIteration(A, b, x);
            steps.push_back(x);

            x = Slae::Solvers::GaussSeidelIterationUpsideDown(A, b, x);
            r = (A * x) - b;

            steps.push_back(x);
        }
        return {x, steps};
    }

    template<typename T>
    std::pair<vector<T>, vector<vector<T>>>
    FastGaussSeidelWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b,
                             const std::vector<T> &initial, T tolerance, T rho) {
        vector<T> previous = initial;
        vector<T> current = Slae::Solvers::GaussSeidelIterationUpsideDown(A, b,
                                                                          Slae::Solvers::GaussSeidelIteration(A, b,
                                                                                                              initial));
        vector<T> next = current;
        vector<vector<T>> steps{};
        steps.reserve(100);
        steps.push_back(initial);
        steps.push_back(initial);

        vector<T> r = A * next - b;
        auto mu0 = static_cast<T>(1.);
        auto mu1 = static_cast<T>(1. / rho);
        auto mu = (static_cast<T>(2.) / rho) * mu1 - mu0;

        while (norm(r) > tolerance) {
            next = ((static_cast<T>(2.) * mu1) / (rho * mu)) *
                   Slae::Solvers::GaussSeidelIterationUpsideDown(A, b,
                                                                 Slae::Solvers::GaussSeidelIteration(A, b, current)) -
                   (mu0 / mu) * previous;
            r = (A * next) - b;

            // refresh all states
            previous = current;
            current = next;
            mu0 = mu1;
            mu1 = mu;
            mu = (static_cast<T>(2.) / rho) * mu1 - mu0;

            // add current step
            steps.push_back(current);
            steps.push_back(current);
        }
        return {next, steps};
    }

}
#endif//SLAE_GAUSSSEIDEL_HPP
