//
// Created by petrov on 12.02.2022.
//

#ifndef SLAE_DENSEMATRIX_HPP
#define SLAE_DENSEMATRIX_HPP
#include <vector>
#include <set>
#include "../utility/Element.hpp"
template<typename T>
class DenseMatrix{
public:
    using elm_t = T;
    using idx_t = std::size_t;

private:

    std::vector<T> matrix;
    idx_t H, W;

public:
    DenseMatrix(const idx_t &h, const idx_t& w);

    DenseMatrix(const idx_t &h, const idx_t& w, const std::set<Triplet<T>>& in);

    elm_t& operator()(const idx_t& i, const idx_t& j);
    const elm_t& operator()(const idx_t& i, const idx_t& j) const;
    [[nodiscard]] const idx_t& sizeH() const;
    [[nodiscard]] const idx_t& sizeW() const;

    void swap(const idx_t& first, const idx_t& second);

    void deleteLastRow();

};
#endif//SLAE_DENSEMATRIX_HPP
