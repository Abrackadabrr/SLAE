//
// Created by evgen on 19.02.2022.
//
#ifndef CSR_EXCEPTION
#define CSR_EXCEPTION

#include <exception>
#include <sstream>
#include <string>
#include "my_project/SlaeBaseException.hpp"


namespace Slae::Matrix {
    class CSRException : public Slae::SlaeBaseExceptionCpp {
    public:
        /** Constructor (C++ STL strings).
           *  @param message The error message.
         */
        explicit CSRException(const std::string &message) noexcept;

        /** Constructor (C strings).
           *  @param message C-style string error message.
           *                 The string contents are copied upon construction.
           *                 Hence, responsibility for deleting the char* lies
           *                 with the caller.
         */
        explicit CSRException(const char *message) noexcept;

        /** Returns a pointer to the (constant) error description.
           *  @return A pointer to a const char*. The underlying memory
           *          is in posession of the Exception object. Callers must
           *          not attempt to free the memory.
         */
    };

} // namespace Slae::Matrix

#endif // CSR_EXCEPTION
