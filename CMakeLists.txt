cmake_minimum_required (VERSION 3.16)
project (blue_mapper VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME} src/blue_mapper.cpp)
