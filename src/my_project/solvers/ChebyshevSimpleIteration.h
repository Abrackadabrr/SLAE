//
// Created by evgen on 20.03.2022.
//

#ifndef MY_PROJECT_CHEBYSHEVSIMPLEITERATION_H
#define MY_PROJECT_CHEBYSHEVSIMPLEITERATION_H

#include <vector>
#include <cmath>
#include "my_project/vectors/Norm.hpp"

template<typename T, int powOf2>
class ChebyshevPolynomial {
private:
public:
    std::vector<T> CalculateRoots(std::pair<T, T> section) {
        int polyOrder = std::pow(2, powOf2);
        std::vector<T> roots(polyOrder);
        for (int i = 1; i <= polyOrder; i++) {
            roots[i - 1] = (section.first + section.second) / 2. + (section.second - section.first) / 2 *
                                                                   static_cast<T>(std::cos(
                                                                           (2 * i - 1) * M_PI_2 / polyOrder));
        }
        std::vector<int> idx = {0, 1};
        std::vector<int> next_idx;
        int currOrder = 2;
        for (int i = 1; i < powOf2; ++i) {
            currOrder = pow(2, i + 1); // CHANGE !!!
            next_idx.resize(currOrder);
            for (int j = 0; j < currOrder - 1; j += 2) {
                next_idx[j] = idx[j / 2];
                next_idx[j + 1] = currOrder - 1 - next_idx[j];
            }
            idx = next_idx;
        }

        std::vector<T> result(roots.size());
        for (int i = 0; i < result.size(); ++i) {
            result[i] = roots[idx[i]];
        }
        return result;
    }
};

namespace Slae::Solvers {
    template<typename T, int powOf2>
    std::vector<T>
    FastSimpleIteration(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initialState,
                        const std::pair<T, T> &borders, const T &tolerance) {

        ChebyshevPolynomial<T, powOf2> pol;
        std::vector<T> tao_roots = pol.CalculateRoots(borders);
        std::vector<T> currentState = initialState;
        std::vector<T> r = A * currentState - b;
        while (norm(r) > tolerance) {
            for (int i = 0; i < tao_roots.size(); i++) {
                currentState = currentState - 1 / tao_roots[i] * r;  // CHANGE
                r = A * currentState - b;
            }
        }
        return currentState;
    }


    template<typename T, int powOf2>
    std::pair<std::vector<T>, std::vector<std::vector<T>>>
    FastSimpleIterationWithSteps(const Slae::Matrix::CSR<T> &A, const std::vector<T> &b,
                                 const std::vector<T> &initialState,
                                 const std::pair<T, T> &borders, const T &tolerance) {

        ChebyshevPolynomial<T, powOf2> pol;
        std::vector<std::vector<T>> steps{};
        steps.reserve(1000);
        steps.push_back(initialState);
        std::vector<T> tao_roots = pol.CalculateRoots(borders);
        std::vector<T> currentState = initialState;
        std::vector<T> r = A * currentState - b;
        while(norm(r) > tolerance) {
            for (int i = 0; i < tao_roots.size(); i++) {
                currentState = currentState - (1 / tao_roots[i]) * r;
                steps.push_back(currentState);
                r = A * currentState - b;
            }
        }

        return {currentState, steps};
    }
};  //Slae::Solvers

#endif //MY_PROJECT_CHEBYSHEVSIMPLEITERATION_H
