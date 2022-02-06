//
// Created by kuznetsov on 04.02.2022.
//

#include "gtest/gtest.h"
#include <my_project/SlaeBaseException.hpp>

void throwException() { throw Slae::SlaeBaseExceptionCpp("Hi"); }

TEST(EXCEPTION, EXCEPTION_HI) {
  bool isCought = false;
  try {
    throwException();
  } catch (const Slae::SlaeBaseExceptionCpp &err) {
    isCought = true;
  }
  ASSERT_TRUE(isCought);
}
