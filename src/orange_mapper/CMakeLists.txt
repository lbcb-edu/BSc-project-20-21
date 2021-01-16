cmake_minimum_required(VERSION "3.10.0")

project(orange_mapper VERSION 0.1.0)
add_executable(${PROJECT_NAME} orange_mapper.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ZLIB)               
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

add_subdirectory(orange_alignment)
add_subdirectory(orange_minimizer)

target_link_libraries(${PROJECT_NAME} bioparser orange_alignment)
target_link_libraries(${PROJECT_NAME} bioparser orange_minimizer)

# Setup testing
#enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIR})

#add_executable(_test
 #   orange_mapper_test.cpp
#)
#target_link_libraries(test orange_mapper GTest::Main)

