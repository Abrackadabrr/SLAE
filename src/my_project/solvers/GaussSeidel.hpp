//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_GAUSSSEIDEL_HPP
#define SLAE_GAUSSSEIDEL_HPP
#include "../sparse/CSR.hpp"

template<typename T>
std::vector<T> GaussSeidel(const CSR<T> &A, const std::vector<T> &b);

#endif//SLAE_GAUSSSEIDEL_HPP
