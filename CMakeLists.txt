cmake_minimum_required(VERSION "3.19.0")

project(orange_mapper VERSION 0.1.0)

add_executable(${PROJECT_NAME} src/orange_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(orange_mapper bioparser)