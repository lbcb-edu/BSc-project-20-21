cmake_minimum_required (VERSION 3.12)
project (blue_mapper VERSION 0.9.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_executable("${PROJECT_NAME}" src/blue_mapper.cpp)
configure_file(include/blue_mapper.hpp.in include/blue_mapper.hpp)

target_include_directories(
  "${PROJECT_NAME}" PUBLIC
  "${PROJECT_BINARY_DIR}/include"
  )

add_subdirectory(blue_minimizers)
add_subdirectory(blue_alignment)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(blue_mapper PUBLIC blue_minimizers blue_alignment bioparser)


option(blue_alignment_build_tests "Build unit tests for blue alignment" ON)
if (blue_alignment_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_alignment_test
    test/blue_alignment_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_alignment_test PUBLIC blue_alignment GTest::Main)
  
endif()

option(blue_minimizers_build_tests "Build unit tests for blue minimizers" ON)
if (blue_minimizers_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_minimizers_test
    test/blue_minimizers_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_minimizers_test PUBLIC blue_minimizers GTest::Main)
  
endif()
