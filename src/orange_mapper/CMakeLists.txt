cmake_minimum_required(VERSION "3.10.0")

project(orange_mapper VERSION 0.1.0)

add_executable(orange_mapper orange_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ZLIB)               
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()


