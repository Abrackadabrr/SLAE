//
// Created by evgen on 05.02.2022.
//

#ifndef MY_PROJECT_THREEDIAGONALMATRIX_HPP
#define MY_PROJECT_THREEDIAGONALMATRIX_HPP

#include <vector>
#include <array>
#include <sstream>
#include "my_project/SlaeBaseException.hpp"
#include "my_project/Defines.h"

namespace Slae::Matrix{
    class ThreeDiagonalMatrix {
    private:
        std::vector<std::array<double, 3>> data_;
    public:
        static ThreeDiagonalMatrix Zero(int n);
        static ThreeDiagonalMatrix Identity(int n);
        static ThreeDiagonalMatrix ThreeDiagonal(int n, double val1, double val2, double val3);
        
        double & operator()(int i, int j);
        [[nodiscard]] const double & operator()(int i, int j) const;
        explicit ThreeDiagonalMatrix(int n);
        [[nodiscard]] int rows() const noexcept;
    };
}  //  namespace Slae::Matrix

#endif //MY_PROJECT_THREEDIAGONALMATRIX_HPP
