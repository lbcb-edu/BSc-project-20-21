cmake_minimum_required(VERSION 3.16.3)

project(white_mapper VERSION 1.0)

configure_file(src/projectControl.h.in src/projectControl.h)

add_library(white_align src/white_alignment.hpp src/white_alignment.cpp)
add_library(white_minim src/white_minimizers.hpp src/white_minimizers.cpp)

add_executable(white_mapper src/white_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
target_include_directories("${PROJECT_NAME}" PUBLIC "${PROJECT_BINARY_DIR}")

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(white_mapper bioparser white_align white_minim)

option(align_tests "Align tests" ON)
if (align_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(white_mapper_align_test test/white_alignment_test.cpp)
  target_link_libraries(white_mapper_align_test PUBLIC white_align GTest::Main)
endif()