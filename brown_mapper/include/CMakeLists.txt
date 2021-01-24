add_library(brown_align brown_alignment.hpp brown_alignment.cpp)
add_library(brown_minim brown_minimizer.hpp brown_minimizer.cpp)

if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()

if (NOT TARGET thread_pool)
  add_subdirectory(thread_pool EXCLUDE_FROM_ALL)
endif ()

