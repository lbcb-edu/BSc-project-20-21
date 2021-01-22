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



option(orange_minimizers_build_tests "Build unit tests for orange minimizers" ON)
if (orange_minimizers_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_minimizers_test
    test/orange_minimizers_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_minimizers_test PUBLIC orange_minimizers GTest::Main)
  
endif()

option(orange_alignment_build_tests "Build unit tests for orange alignment" ON)
if (orange_alignment_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_alignment_test
    test/orange_alignment_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_alignment_test PUBLIC orange_alignment GTest::Main)
  
endif()


