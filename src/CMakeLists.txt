add_executable(blonde_mapper blonde_mapper.cpp)

add_library(blonde_alignment blonde_alignment.h blonde_alignment.cpp)

target_link_libraries(
	blonde_mapper 
		blonde_alignment
		bioparser
)