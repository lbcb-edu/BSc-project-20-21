cmake_minimum_required(VERSION 3.9)

project(blonde_mapper LANGUAGES CXX)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT TARGET bioparser)
  add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif ()

# Project subdirs 
add_subdirectory(src)
# add_subdirectory(test)

option(blonde_alignment_tests "blonde alignment tests" ON)
if (blonde_alignment_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_test test/blonde_alignment_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test PUBLIC blonde_alignment GTest::Main)
endif()