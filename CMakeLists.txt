cmake_minimum_required(VERSION 3.18.2)

project(white_mapper VERSION 1.0)

configure_file(src/projectControl.h.in src/projectControl.h)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(${PROJECT_NAME} src/white_mapper.cpp)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(white_mapper bioparser)