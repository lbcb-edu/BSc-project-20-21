add_executable(blonde_mapper blonde_mapper.cpp)

add_library(blonde_alignment blonde_alignment.h blonde_alignment.cpp)
add_library(blonde_minimizers blonde_minimizers.h blonde_minimizers.cpp)

target_link_libraries(
	blonde_mapper 
		blonde_alignment
		blonde_minimizers
		bioparser
)
