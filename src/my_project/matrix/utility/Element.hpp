//
// Created by evgen on 12.02.2022.
//

#ifndef SLAE_ELEMENT_HPP
#define SLAE_ELEMENT_HPP
#include <cstdio>

namespace Slae::Matrix {
    template<typename T>
    struct Element {
        std::size_t i;
        std::size_t j;
        T value;

        bool operator<(Element<T> const &rgh) const;
    };

    template<typename T>
    bool Element<T>::operator<(const Element<T> &rgh) const {
        return i < rgh.i or (i == rgh.i and this->j < rgh.j);
    }
}  //  Slae::Matrix

#endif//SLAE_ELEMENT_HPP
