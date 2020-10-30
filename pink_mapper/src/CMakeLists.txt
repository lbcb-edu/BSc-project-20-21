cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(pink_mapper VERSION 0.1.0)

configure_file(Pink_mapperConfig.h.in Pink_mapperConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(pink_mapper pink_mapper.cpp)

target_include_directories(pink_mapper PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

find_package(ZLIB)
                          
target_link_libraries(pink_mapper ZLIB::ZLIB)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/help.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)


