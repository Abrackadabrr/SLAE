//
// Created by evgen on 05.02.2022.
//

#ifndef MY_PROJECT_THREEDIAGONALSOLVER_HPP
#define MY_PROJECT_THREEDIAGONALSOLVER_HPP

#include "my_project/matrix/ThreeDiagonalMatrix.hpp"
#include "my_project/defines.h"

namespace Slae::Solvers {
    std::vector<double>
    solveThreeDiagonal(const Slae::Matrix::ThreeDiagonalMatrix &matrix, const std::vector<double> &col);
}

#endif //MY_PROJECT_THREEDIAGONALSOLVER_HPP
