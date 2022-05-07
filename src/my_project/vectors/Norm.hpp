//
// Created by petrov on 19.02.2022.
//

#ifndef SLAE_NORM_HPP
#define SLAE_NORM_HPP
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

/**
 * Class describing type norms
 */
enum class NormType{
    InfNorm = 0,
    FirstNorm = 1,
    SecondNorm = 2
};

/**
 * @first, second - vectors
 * @return value_t that is result of scalar producing of vectors
 */
template<typename value_t>
value_t operator*(const std::vector<value_t>& first, const std::vector<value_t>& second) {
    return std::inner_product(first.begin(), first.end(), second.begin(), static_cast<value_t>(0),\
    std::plus<>(), std::multiplies<>());
}

/**
 *
 * @tparam value_t  type of elements in vector
 * @param v vector
 * @param type_n type of norm to calculate
 * @return value_t figure that describe type_n norm
 */
template<typename value_t>
value_t norm(const std::vector<value_t> v, NormType type_n=NormType::InfNorm){
    auto norm = static_cast<value_t>(0);
    if (type_n == NormType::InfNorm) {
        auto abs_compare = [](value_t a, value_t b) {return std::abs(a) < std::abs(b);};
        norm = std::abs(*std::max_element(v.begin(), v.end(), abs_compare));
    }
    if (type_n == NormType::FirstNorm) {
        auto addition = [](value_t a, value_t b) {return std::abs(a) + std::abs(b);};
        norm = std::accumulate(v.begin(), v.end(), static_cast<value_t>(0), addition);
    }
    else if (type_n == NormType::SecondNorm) {
        norm = std::sqrt(v*v);
    }
    return norm;
}

template<NormType type, typename value_t>
value_t Norm(const std::vector<value_t>& v) {return static_cast<value_t>(0);};

template<>
double Norm<NormType::InfNorm>(const std::vector<double>& v){
    auto abs_compare = [](double a, double b) {return std::abs(a) < std::abs(b);};
    return std::abs(*std::max_element(v.begin(), v.end(), abs_compare));
}

template<>
double Norm<NormType::FirstNorm>(const std::vector<double>& v){
    auto addition = [](double a, double b) {return a + std::abs(b);};
    return std::accumulate(v.begin(), v.end(), static_cast<double>(0), addition);
}

template<>
double Norm<NormType::SecondNorm>(const std::vector<double>& v) {
    return std::sqrt(v*v);
}

#endif//SLAE_NORM_HPP
