cmake_minimum_required(VERSION 3.9)

project(blonde_mapper)

add_executable(blonde_mapper blonde_mapper.cpp)
if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(blonde_mapper bioparser)