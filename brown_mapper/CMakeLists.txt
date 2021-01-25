cmake_minimum_required(VERSION 3.10)

project(brown_mapper)

add_executable(brown_mapper brown_mapper.cpp)

add_subdirectory(include)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

target_link_libraries(brown_mapper bioparser thread_pool brown_align brown_minim)

option(align_tests "Align tests" ON)
if (align_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(brown_mapper_align_test tests/brown_align_tests.cpp)
  target_link_libraries(brown_mapper_align_test PUBLIC brown_align GTest::Main)
endif()

option(minimizers_tests "minimizer tests" ON)
if (minimizers_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(brown_mapper_minim_test tests/brown_minimizers_tests.cpp)
  target_link_libraries(brown_mapper_minim_test PUBLIC brown_minim GTest::Main)
endif()
