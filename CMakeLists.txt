cmake_minimum_required(VERSION 3.9)

project(blonde_mapper)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT TARGET bioparser)
  add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif ()

# Project subdirs 
add_subdirectory(src)