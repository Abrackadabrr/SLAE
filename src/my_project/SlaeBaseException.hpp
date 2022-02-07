//
// Created by kuznetsov on 04.02.2022.
//

#ifndef SEM_1_MYEXCEPTION_HPP
#define SEM_1_MYEXCEPTION_HPP

#include <exception>
#include <sstream>
#include <string>


namespace Slae {

class SlaeBaseExceptionCpp : public std::exception {
  /* Базовый класс исключения в пакетах С++ */
protected:
  /* Поля класса */
  std::string message_;

public:
  /* Конструкторы класса */
  /** Constructor (C++ STL strings).
     *  @param message The error message.
   */
  explicit SlaeBaseExceptionCpp(const std::string& message) noexcept;
  /** Constructor (C strings).
     *  @param message C-style string error message.
     *                 The string contents are copied upon construction.
     *                 Hence, responsibility for deleting the char* lies
     *                 with the caller.
   */
  explicit SlaeBaseExceptionCpp(const char* message) noexcept;

  /*  Методы класса */
  /** Returns a pointer to the (constant) error description.
     *  @return A pointer to a const char*. The underlying memory
     *          is in posession of the Exception object. Callers must
     *          not attempt to free the memory.
   */
  [[nodiscard]] const char* what() const noexcept override;
};


} // namespace Slae

#endif // SEM_1_MYEXCEPTION_HPP
