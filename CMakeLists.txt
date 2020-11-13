cmake_minimum_required (VERSION 3.12)
project (blue_mapper VERSION 0.4.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable("${PROJECT_NAME}" src/blue_mapper.cpp)

configure_file(include/blue_mapper.hpp.in include/blue_mapper.hpp)

target_include_directories(
  "${PROJECT_NAME}" PUBLIC
  "${PROJECT_BINARY_DIR}/include"
  "${PROJECT_SOURCE_DIR}/blue_alignment/include"
  )

add_subdirectory(blue_alignment)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(blue_mapper PUBLIC blue_alignment bioparser)
