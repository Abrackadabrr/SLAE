//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_JACOBI_HPP
#define SLAE_JACOBI_HPP
#include "../sparse/CSR.hpp"

template<typename T>
std::vector<T> Yacobi(const CSR<T> &A, const std::vector<T> &b);

#endif//SLAE_JACOBI_HPP
