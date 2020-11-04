cmake_minimum_required(VERSION "3.15.0")

project(orange_mapper VERSION 1.0)

configure_file(
    "src/config.h.in"
    "src/config.h"
)

add_executable(${PROJECT_NAME} src/orange_mapper.cpp)

find_package(ZLIB REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE
  ZLIB::ZLIB)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(orange_mapper bioparser)