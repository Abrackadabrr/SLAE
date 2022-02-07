cmake_minimum_required(VERSION 3.17)
project(my_project)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
message("TESTING: ${WITH_TESTS}")
