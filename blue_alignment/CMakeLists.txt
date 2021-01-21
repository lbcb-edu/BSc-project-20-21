add_library(blue_alignment "blue_alignment.cpp")
target_include_directories(blue_alignment PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

option(blue_alignment_build_tests "Build unit tests for blue alignment" ON)
if (blue_alignment_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_alignment_test
    test/blue_alignment_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_alignment_test PUBLIC blue_alignment GTest::Main)
  
endif()