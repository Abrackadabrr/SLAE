//
// Created by evgen on 06.02.2022.
//

#ifndef MY_PROJECT_VECTORS_H
#define MY_PROJECT_VECTORS_H

#include <algorithm>
#include <fstream>
#include <numeric>

/*
 * Returns value_t that is result of scalar producing of vectors
 */
    template<typename value_t>
    value_t operator*(const std::vector<value_t>& first, const std::vector<value_t>& second) {
        return std::inner_product(first.begin(), first.end(), second.begin(), 0,
                                  std::plus<value_t>(), std::multiplies<value_t>());
    }

/*
 * Returns sum of vectors like in linear algebra
 */
    template<typename value_t>
    std::vector<value_t> operator+(const std::vector<value_t>& first, const std::vector<value_t>& second) {
        std::vector<value_t> result(first.size());
        std::transform(first.begin(), first.end(), second.begin(), result.begin(), std::plus<>());
        return result;
    }

/*
 * Production of vector and number
 */
    template<typename value_t>
    std::vector<value_t> operator*(const std::vector<value_t>& first, value_t number) {
        std::vector<value_t> result;
        result.reserve(first.size());
        for (auto& i : first)
            result.push_back(i*number);
        return std::move(result);
    }

/*
 * Production of number and vector
 */
    template<typename value_t>
    value_t operator*(value_t number, const std::vector<value_t>& first) {
        return first*number;
    }

/*
 * Returns sequence of instance between start ans stop with step between each neighboring elements
 */
    template<typename value_type>
    std::vector<value_type> sequence(value_type start, value_type stop, value_type step)
    {
        std::vector<value_type> result;
        result.reserve(static_cast<int>((stop - start)/step));
        for (int i = 0, stop_int = static_cast<int>((stop - start)/step); i < stop_int; ++i) {
            result.emplace_back(static_cast<value_type>(start + i * step));
        }
        return std::move(result);
    }

/*
 * Output operator
 * value_t is required to have the similar output operator
 */
    template<typename value_t>
    std::ostream &operator << (std::ostream &os, const std::vector<value_t>& v)
    {
        os << "[";
        for(int i = 0; i < v.size() - 1; i++)
        {os << v[i] <<" ";};
        os << v[v.size() - 1];
        os << ']';
        return os;
    }

/*
 * Output function for csv type
 * value_t is required to have an output operator
 */
    template<typename value_t>
    void printCsvFormat(std::ostream &os, const std::vector<value_t>& v)
    {
        for(int i = 0; i < v.size() - 1; i++)
        {os << v[i] <<";";};
        os << v[v.size() - 1];
    }


#endif //MY_PROJECT_VECTORS_H
