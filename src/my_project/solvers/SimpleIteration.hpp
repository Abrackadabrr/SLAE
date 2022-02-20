//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_SIMPLEITERATION_HPP
#define SLAE_SIMPLEITERATION_HPP
#include "../sparse/CSR.hpp"
template<typename T>
std::vector<T> SimpleIteration(const CSR<T> &A, const std::vector<T> &b, const T &tao);

#endif//SLAE_SIMPLEITERATION_HPP
