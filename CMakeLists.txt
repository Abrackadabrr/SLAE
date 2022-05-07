cmake_minimum_required(VERSION 3.17)
project(my_project)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)
message("TESTING: ${WITH_TESTS}")
if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

message("PLOTTING: ${WITH_PLOTS}")
if (WITH_PLOTS)
    enable_testing()
    add_subdirectory(plots)
endif ()

message("VTK_EXAMPLES: ${WITH_PLOTS}")
if (WITH_PLOTS)
    enable_testing()
    add_subdirectory(vtk_examples)
endif ()
