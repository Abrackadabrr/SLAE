//
// Created by evgen on 03.05.2022.
//

#ifndef MY_PROJECT_GAUSS_H
#define MY_PROJECT_GAUSS_H

#include "my_project/matrix/sparse/CSRmatrix.h"

template<typename T, typename MatrixT>
std::vector<T> GaussDownshift(MatrixT A, std::vector<T> b) {
    int n = b.size();
    std::vector<T> res(n, 0);
    res[n - 1] = b[n - 1] / A(n - 1, n - 1);
    for (int i = n - 2; i >= 0; i--) {
        T sum = static_cast<T>(b[i]);
        for (int k = n - 1; k > i; k--) {
            sum -= A(i, k) * res[k];
        }
        res[i] = (1 / A(i, i)) * sum;
    }
    return res;
}

#endif //MY_PROJECT_GAUSS_H
