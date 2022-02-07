//
// Created by evgen on 05.02.2022.
//

#include "ThreeDiagonalSolver.hpp"


std::vector<double>
Slae::Solvers::solveThreeDiagonal(const Slae::Matrix::ThreeDiagonalMatrix &matrix, const std::vector<double> &right_side) {
    int n = matrix.rows();

    if (n != right_side.size()) {
        std::stringstream buff;
        buff << "Matrix's and vector's size arent same." << ERR_INFO;
        throw SlaeBaseExceptionCpp(buff.str());
    }

    std::vector<std::array<double, 2>> extra_params{};
    extra_params.reserve(n);

    std::array<double, 2> pre_params{- matrix(0, 2) / matrix(0, 1),
                                     right_side[0] / matrix(0, 1)};

    extra_params.push_back(pre_params);  // this operation doesnt take a lot of time because extra_params.capacity == n

    for (int i = 1; i < n - 1; ++i) {
        pre_params = {- matrix(i, 2) / (matrix(i, 0) * pre_params[0] + matrix(i, 1)),
                      (right_side[i] - matrix(i, 0) * pre_params[1]) / (matrix(i, 0) * pre_params[0] + matrix(i, 1))};
        extra_params.push_back(pre_params);
    }

    std::vector<double> result(n);
    result[n - 1] = (right_side[n - 1] - matrix(n - 1, 0) * pre_params[1]) / (matrix(n - 1, 0) * pre_params[0] + matrix(n - 1, 1));
    for (int i = n - 2; i >= 0; i--) {
        const auto& params = extra_params[i];
        result[i] = params[0] * result[i + 1] + params[1];
    }
    return result;
}
