cmake_minimum_required(VERSION 3.2.2)

# set the project name and version
cmake_policy(SET CMP0048 NEW)
project(pink_mapper VERSION 0.2.0)

configure_file(Pink_mapperConfig.h.in Pink_mapperConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(ALIGN_INCLUDES ${PROJECT_SOURCE_DIR}/../include/align)
set(PARSER_INCLUDES ${PROJECT_SOURCE_DIR}/../include/bioparser)
file(GLOB SRC_FILES *.cpp)
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/../test/*.cpp)

# add the executable
option(build_mapper "Build build mapper" ON)
if(build_mapper)

include_directories(${PARSER_INCLUDES})
include_directories(${ALIGN_INCLUDES})

add_executable(${PROJECT_NAME} pink_mapper.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

find_package(ZLIB)               
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
endif()


option(alignment_build_tests "Build alignment unit tests" OFF)

if (alignment_build_tests)
  enable_testing()
  find_package(GTest REQUIRED)
  

  
  include_directories(${ALIGN_INCLUDES})
  
  
  add_executable(pink_alignment_test ${TEST_SRC_FILES})
  target_link_libraries(pink_alignment_test ${GTEST_LIBRARIES})
endif ()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/help.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)