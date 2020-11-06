add_executable(blonde_mapper blonde_mapper.cpp)

add_library(blonde_alignment INTERFACE)

target_link_libraries(
	blonde_mapper 
		blonde_alignment
		bioparser
)