## Dependencies
find_package(MPI)


set(PARALLEL_SOURCES
  cfmpi.cpp
)

## ######################
## Build rules
## ######################
add_library(parallel SHARED ${PARALLEL_SOURCES})


target_include_directories(parallel
                           PRIVATE 
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           ${MPI_CXX_INCLUDE_PATH})

