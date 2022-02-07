cmake_minimum_required(VERSION 3.17)
project(my_project)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE MY_PROJECT_SOURCES *.h *.cpp *.hpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

add_library(my_project ${MY_PROJECT_SOURCES})
target_link_libraries(my_project PUBLIC )

target_include_directories(my_project PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

