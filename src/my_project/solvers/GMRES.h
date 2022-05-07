//
// Created by evgen on 30.04.2022.
//

#ifndef MY_PROJECT_GMRES_H
#define MY_PROJECT_GMRES_H

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "my_project/matrix/dense/Densematrix.hpp"

template<typename T>
bool isOrtoSystem(std::vector<std::vector<T>> vectors) {
    for (int i = 0; i < vectors.size(); i++)
        for (int j = 0; j < vectors.size(); j++)
            if (i != j && std::abs(vectors[i] * vectors[j]) > 1e-10) return false;
    return true;
}

template<typename T>
void NextBasisVector(Slae::Matrix::DenseMatrix<T> &V, Slae::Matrix::DenseMatrix<T> &H, const Slae::Matrix::CSR<T> &A,
                     int i) {
    auto t = A * V.get_col(i - 1);
    std::vector<double> sum(V.get_row_size(), 0);
    for (int k = 0; k < i; k++) {
        auto v_k = V.get_col(k);
        H(k, i - 1) = v_k * t;
        sum = sum + H(k, i - 1) * v_k;
    }
    t = t - sum;
    H(i, i - 1) = Norm<NormType::SecondNorm>(t);
    if (std::abs(H(i, i - 1)) > 1e-10)
        V.write_col(t * (1 / H(i, i - 1)), i);
}

template<typename T>
void Rotate(Slae::Matrix::DenseMatrix<T> &H, std::vector<T> &right_side, int ind,
            std::vector<std::array<T, 2>> &rotations_before, T tolerance) {
    int i = ind;
    int j = i + 1;

    for (int k = 0; k < i; k++) {
        T xi = H(k, i);
        T xj = H(k + 1, i);
        T cos = rotations_before[k][0];
        T sin = rotations_before[k][1];
        H(k, i) = xi * cos - xj * sin;
        H(k + 1, i) = xi * sin + xj * cos;
    }

    T xi = H(i, i);
    T xj = H(j, i);

    T cos, sin;
    if (std::abs(xj) > 1e-10) {
        T relationship = -xj / xi;
        T sign = relationship / std::abs(relationship);
        cos = std::sqrt(1 / (1 + relationship * relationship));
        sin = sign * std::sqrt(1 / (1 + ((1 / relationship) * (1 / relationship))));
    } else {
        rotations_before.push_back({static_cast<T>(1), static_cast<T>(0)});
        return;
    }
    H(i, i) = xi * cos - xj * sin;
    H(j, i) = xi * sin + xj * cos; /* H(j, i) = 0; */
    rotations_before.push_back({cos, sin});
    xi = right_side[i];
    xj = right_side[j];
    right_side[i] = xi * cos - xj * sin;
    right_side[j] = xi * sin + xj * cos;
}

namespace Slae::Solvers {
/***
 * Generalised minimal residual method iteration
 * @tparam T - value type
 * @param m - number of iterations
 * @param A - matrix of slae
 * @param b - right size of slae
 * @param initial - initial approximation
 * @param tolerance - tolerance
 * @return solution to slae
 */
    template<typename T>
    std::pair<std::vector<T>, T>
    GMRESIteration(int m, const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial, T tolerance) {
        int n = b.size();
        Slae::Matrix::DenseMatrix<T> H(m + 1, m, 0);
        Slae::Matrix::DenseMatrix<T> V(n, m + 1, 0);
        auto r0 = A * initial - b;
        std::vector<std::array<double, 2>> rotations{};
        rotations.reserve(m);
        std::vector<T> right_side(m + 1, 0);
        right_side[0] = Norm<NormType::SecondNorm>(r0);
        V.write_col(r0 * (1 / right_side[0]), 0);
        int i = 1;
        T normr = std::abs(right_side[0]);
        while(i < m+1 && normr > tolerance) {
            NextBasisVector(V, H, A, i);
            Rotate(H, right_side, i - 1, rotations, tolerance);
            normr = std::abs(right_side[i]);
            i++;
        }
        auto e = std::vector<T>(right_side.begin(), right_side.begin() + i - 1);
        std::vector<double> y = GaussDownshift(H, e);
        std::vector<double> res(n, 0);
        for (int k = 0; k < y.size(); k++) {
            res = res + y[k] * V.get_col(k);
        }
        return {initial - res, normr};
    }

    template<typename T>
    std::vector<T> GMRESM(int m, const Slae::Matrix::CSR<T> &A, const std::vector<T> &b, const std::vector<T> &initial, T tolerance){
        std::pair<std::vector<T>, T> current = GMRESIteration(m , A, b, initial, tolerance);
        while (current.second > tolerance){
            current = GMRESIteration(m, A, b, current.first, tolerance);
        }
        return current.first;
    }

    template<typename T, int m>
    std::vector<T> GMRESwithDebug(Slae::Matrix::CSR<T> A, std::vector<T> b, std::vector<T> initial, T tolerance) {
        int n = b.size();
        Slae::Matrix::DenseMatrix<T> H(m + 1, m, 0);  // hessenberg matrix
        Slae::Matrix::DenseMatrix<T> V(n, m + 1, 0);  // matrix contains basis vectors in krylov subspace

        auto r0 = A * initial - b;  // initial residual operation takes O(amount of A.values.size)

        std::vector<std::array<double, 2>> rotations{};
        rotations.reserve(m);

        // creation right side of Least Square problem
        // initial is xi0
        std::vector<T> right_side(m + 1, 0);
        right_side[0] = norm(r0, NormType::SecondNorm);

        V.write_col(r0 * (1 / right_side[0]), 0);  // first vector in krylov subspace(1), op takes O(n)

        // start loop for m iterations of gmres
        // i means number of xi. so xi_i - xi_0 is in krylov_subspace_i
        // H_i matches krylov_subspace_i but it is required to evaluate basis in krylov_subspace_i+1 to find H_i
        for (int i = 1; i < m + 1; i++) {
            // build H matrix
            NextBasisVector(V, H, A, i);
            // now i have actual V(n x i+1) - basis in krylov subspace(i+1)
            // and H(i+1, i) matrix - 'transition matrix' between Vi (previous step) and Vi+1 that i currently have.
            // H is matrix i want to solve LS problem. Keep in mind that its actual size is (i+1 x i)
            // To evaluate norm(ri) i should rotate H and e1 vector
            // Due to i is amount of columns in current H i should rotate i-1 column.
            std::cout << "Iteration " << i << '\n' << H << '\n';
            Rotate(H, right_side, i - 1, rotations, tolerance);
            std::cout << "Residual is " << std::abs(right_side[i]) << '\n';
//            if (std::abs(right_side[i]) < tolerance) break;
            std::cout << "Basis i+1\n" << V << '\n' /*<< "RS:" << right_side << '\n'*/;

            // checking is system right
            std::vector<std::vector<T>> vectors(i);
            for (int k = 0; k < i; k++) {
                vectors[k] = V.get_col(k);
            }
//            std::cout << "Is system orthogonal " << std::boolalpha << isOrtoSystem(vectors) << '\n';
        }
        H.deleteLastRow();
        auto e = std::vector<T>(right_side.begin(), right_side.begin() + m);
        std::vector<double> y = GaussDownshift(H, e);
        std::vector<double> res(n, 0);
        for (int i = 0; i < y.size(); i++) {
            res = res + y[i] * V.get_col(i);
        }
        return initial - res;
    }
};  // Slae::Solvers


#endif //MY_PROJECT_GMRES_H
