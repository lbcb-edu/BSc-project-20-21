cmake_minimum_required(VERSION 3.10)

project(brown_mapper)

add_library(brown_align brown_alignment.hpp brown_alignment.cpp)
add_library(brown_minim brown_minimizer.hpp brown_minimizer.cpp)

add_executable(brown_mapper brown_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
target_include_directories("${PROJECT_NAME}" PUBLIC "${PROJECT_BINARY_DIR}")

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(brown_mapper bioparser brown_align brown_minim)

option(align_tests "Align tests" ON)
if (align_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(brown_mapper_test tests/brown_align_tests.cpp)
  target_link_libraries(brown_mapper_test PUBLIC brown_align GTest::Main)
endif()
