cmake_minimum_required(VERSION 3.2.2)

# set the project name and version
project(pink_mapper VERSION 0.3.0)

configure_file(Pink_mapperConfig.h.in Pink_mapperConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(ALIGN_INCLUDES ${PROJECT_SOURCE_DIR}/../include/align)
set(MINIMIZERS_INCLUDES ${PROJECT_SOURCE_DIR}/../include/minimizers)
set(PARSER_INCLUDES ${PROJECT_SOURCE_DIR}/../include/bioparser)
set(HELP_INCLUDES ${PROJECT_SOURCE_DIR}/../include/help)
file(GLOB SRC_FILES *.cpp)
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/../test/*.cpp)

# add the executable
option(build_mapper "Build mapper" ON)
if(build_mapper)

include_directories(${PARSER_INCLUDES})
include_directories(${ALIGN_INCLUDES})
include_directories(${MINIMIZERS_INCLUDES})
include_directories(${HELP_INCLUDES})

add_executable(${PROJECT_NAME} pink_mapper.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

find_package(ZLIB)               
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
endif()


option(mapper_build_tests "Build mapper unit tests" OFF)

if (mapper_build_tests)
  find_package(GTest REQUIRED)
  enable_testing()

  
  include_directories(${ALIGN_INCLUDES})
  include_directories(${MINIMIZERS_INCLUDES})
  
  
  add_executable(pink_mapper_test ${TEST_SRC_FILES})
  target_link_libraries(pink_mapper_test
    GTest::Main)

  add_test(UnitTests pink_mapper_test)
endif ()





