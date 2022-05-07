//
// Created by evgen on 19.02.2022.
//

#include "CSRExeption.h"

namespace Slae::Matrix {
    CSRException::CSRException(const char *message) noexcept
            : SlaeBaseExceptionCpp(message){}

    CSRException::CSRException(const std::string &message) noexcept
            : SlaeBaseExceptionCpp(message) {}

} // namespace Slae
