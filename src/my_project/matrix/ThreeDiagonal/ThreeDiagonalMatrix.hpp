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
#include "my_project/matrix/sparse/CSRmatrix.h"
#include "set"
#include "my_project/matrix/utility/Element.hpp"

namespace Slae::Matrix{
    /**
     * @brief Class of special type of matrix: matrix with three diagonal
     *
     * Elements in matrix almost equal zero, except of elements on the main diagonal
     * and diagonals that are above and below it. So there is a special type of holding these values.
     */
    class ThreeDiagonalMatrix {
    private:
        std::vector<std::array<double, 3>> data_;
    public:

        /** @brief ThreeDiagonalMatrix static 'constructor'
         * Creates three-diagonal matrix with size 'size' and fill it by zeros
         *
         * @param size Matrix size
         *
         * @return zero-filled size-sized ThreeDiagonalMatrix instance
         */
        static ThreeDiagonalMatrix Zero(int n);

        /** @brief ThreeDiagonalMatrix static 'constructor'
         * Creates three-diagonal matrix with size 'size' and fill it by ones on main diagonal
         *
         * @param size Matrix size
         *
         * @return ThreeDiagonalMatrix instance
         */
        static ThreeDiagonalMatrix Identity(int n);

        /** @brief ThreeDiagonalMatrix static 'constructor'
         * Creates three-diagonal matrix with size 'size' and fill it by given values on diagonals
         *
         * @param size Matrix size
         *
         * @param var1 Figure on "under main" diagonal
         *
         * @param var2 Figure on "main" diagonal
         *
         * @param var1 Figure on "up main" diagonal
         *
         * @return ThreeDiagonalMatrix instance
         */
        static ThreeDiagonalMatrix ThreeDiagonal(int n, double val1, double val2, double val3);

       /** @brief ThreeDiagonalMatrix class constructor
        *
        * Creates three-diagonal matrix with size 'size'
        * @warning You shouldn't use this instance without set all of elements with special methods
        * @see  operator()(int i, int j) Zero(int n) ThreeDiagonal(int n, double val1, double val2, double val3)
        *
        * @param size Matrix instance size
        */
        explicit ThreeDiagonalMatrix(int n);

        /** @brief Access operator
         *
         * @return &(i, j) element of three-diagonal matrix; i means number of row, j belongs to {0, 1, 2}
         *
         * @param i  Number of row, belongs to {0,..., matrix size - 1}
         *
         * @param j  One of diagonal elements in a row i; j belongs to {0, 1, 2} === {below, main, above}
         *
         * @warning If you try to access (0, 0) or (matrix size - 1, 2) elements, exception will be thrown
         * because elements with these indexes don't exist.
         *
         * @throw SlaeBaseExceptionCpp if indexes are incorrect
         */
        double & operator()(int i, int j);

        /** @brief Access operator
         *
         * @return const &(i, j) element of three-diagonal matrix; i means number of row, j belongs to {0, 1, 2}
         *
         * @param i  Number of row, belongs to {0,..., matrix size - 1}
         *
         * @param j  One of diagonal elements in a row i; j belongs to {0, 1, 2} === {below, main, above}
         *
         * @warning If you try to access (0, 0) or (matrix size - 1, 2) elements, exception will be thrown
         * because elements with these indexes don't exist.
         *
         * @throw SlaeBaseExceptionCpp if indexes are incorrect
         */
        [[nodiscard]] const double & operator()(int i, int j) const;

        /** @brief
         *
         * @return matrix row-size
        */
        [[nodiscard]] int rows() const noexcept;

        /**
         * @brief Cast to CSR format
         *
         * @return CSR three diagonal matrix
         */
         [[nodiscard]] CSR<double> toCSR(){
             std::set<Slae::Matrix::Element<double>> elements;
             for (unsigned i = 0; i < data_.size(); i++) {
                 if (i != 0 && i != data_.size() - 1) {
                     elements.insert({i, i - 1, data_[i][0]});
                     elements.insert({i, i, data_[i][1]});
                     elements.insert({i, i + 1, data_[i][2]});
                 } else if (i == 0) {
                     elements.insert({i, i, data_[i][1]});
                     elements.insert({i, i + 1, data_[i][2]});
                 } else if (i == data_.size()-1) {
                     elements.insert({i, i - 1, data_[i][0]});
                     elements.insert({i, i, data_[i][1]});
                 }
             }
            return CSR<double>(data_.size(), data_.size(), elements);
         }
    };
}  //  namespace Slae::Matrix

#endif //MY_PROJECT_THREEDIAGONALMATRIX_HPP
