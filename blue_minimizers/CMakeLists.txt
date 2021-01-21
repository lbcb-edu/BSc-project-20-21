add_library(blue_minimizers "blue_minimizers.cpp")
target_include_directories(blue_minimizers PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

option(blue_minimizers_build_tests "Build unit tests for blue minimizers" ON)
if (blue_minimizers_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_minimizers_test
    test/blue_minimizers_test.cpp
  )
  target_link_libraries(${PROJECT_NAME}_minimizers_test PUBLIC blue_minimizers GTest::Main)
  
endif()