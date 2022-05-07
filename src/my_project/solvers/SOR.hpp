//
// Created by vladimir on 30.03.2022.
//

#ifndef MY_PROJECT_SOR_HPP
#define MY_PROJECT_SOR_HPP

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/vectors/Norm.hpp"
#include "my_project/vectors/Vectors.hpp"

namespace Slae::Solvers {
    template<typename T>
    std::vector<T>
    SuccessiveOverRelaxation(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState,
                             const T &tolerance) {
        double w = 0.8;
        std::vector<T> r = A * initialState - b;
        std::vector<T> currentState = initialState;
        T sum;
        auto norma = norm(r);
        while (norm(r) > tolerance) {
            for (int i = 0; i < A.height(); ++i) {
                sum = static_cast<T>(0);
                for (int k = 0; k < A.height(); ++k) {
                    if (k != i) sum += A(i, k) * currentState[k];
                }
                currentState[i] = (1 - w) * currentState[i] + w * (b[i] - sum) / A(i, i);
            }
            r = A * currentState - b;
        }
        return currentState;
    }


    template<typename T>
    std::pair<std::vector<T>, std::vector<std::vector<T>>>
    SuccessiveOverRelaxationWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b,
                                      const std::vector<T> &initialState,
                                      const T &tolerance, double w_) {
        double w = w_;
        std::vector<std::vector<T>> steps{};
        std::vector<T> r = A * initialState - b;
        std::vector<T> currentState = initialState;
        T sum;
        auto norma = norm(r);
        while (norm(r) > tolerance) {
            for (int i = 0; i < A.height(); ++i) {
                sum = static_cast<T>(0);
                for (int k = 0; k < A.height(); ++k) {
                    if (k != i) sum += A(i, k) * currentState[k]; // CHANGE
                }
                currentState[i] = (1 - w) * currentState[i] + w * (b[i] - sum) / A(i, i);
            }
            r = A * currentState - b;
            steps.push_back(currentState);
        }
        return {currentState, steps};
    }
} //Slae::Solvers

#endif //MY_PROJECT_SOR_HPP
