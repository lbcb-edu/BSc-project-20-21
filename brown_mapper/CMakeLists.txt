cmake_minimum_required(VERSION 3.10)

project(brown_mapper)

add_library(brown_align brown_alignment.hpp brown_alignment.cpp)

add_executable(brown_mapper brown_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

target_include_directories("${PROJECT_NAME}" PUBLIC "${PROJECT_BINARY_DIR}")

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(brown_mapper bioparser brown_align)
