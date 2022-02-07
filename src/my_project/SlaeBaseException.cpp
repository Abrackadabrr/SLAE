//
// Created by kuznetsov on 04.02.2022.
//

#include "SlaeBaseException.hpp"

namespace Slae {

SlaeBaseExceptionCpp::SlaeBaseExceptionCpp(const char *message) noexcept
    : message_(message) {}

SlaeBaseExceptionCpp::SlaeBaseExceptionCpp(const std::string &message) noexcept
    : SlaeBaseExceptionCpp(message.c_str()) {}

const char *SlaeBaseExceptionCpp::what() const noexcept { return message_.c_str(); }
} // namespace Slae
