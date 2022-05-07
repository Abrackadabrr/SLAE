//
// Created by evgen on 05.02.2022.
//

#ifndef MY_PROJECT_THREEDIAGONALSOLVER_HPP
#define MY_PROJECT_THREEDIAGONALSOLVER_HPP

#include "my_project/matrix/ThreeDiagonal/ThreeDiagonalMatrix.hpp"
#include "my_project/Defines.h"

namespace Slae::Solvers {
    /** @brief The method solves the system of equations using the sweep method
    * Solves a system of linear algebraic equations using the sweep method. You can learn about the sweep method
    * from ... (source link)
    *
    * @param matrix threediagonal matrix
     *
    * @param right_side column
    *
    * @warning Matrix instance must have the property of diagonal dominance.
    * If it have not, instead of real result of SLAE, you wil be given by nan or smth like that
    *
    * @throw SlaeBaseExceptionCpp is thrown if the number of matrix rows and column height do not match
    * @return The result of SLAE
    */
    std::vector<double>
    solveThreeDiagonal(const Slae::Matrix::ThreeDiagonalMatrix &matrix, const std::vector<double> &right_side);
}

#endif //MY_PROJECT_THREEDIAGONALSOLVER_HPP
