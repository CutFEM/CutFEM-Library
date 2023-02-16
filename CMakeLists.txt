cmake_minimum_required(VERSION 3.20)
project(libCutFEM VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## ##############################
## Options
## ##############################

option(USE_MPI "Perform a parallel computation" ON)
option(USE_OMP "Perform a parallel omp computation" OFF)
option(USE_UMFPACK "Found the UMFpack library" ON)
option(USE_LAPACK "Found the lapacke library" OFF)
option(USE_MUMPS "Found the MUMPS library" ON)

option(CUTFEM_BUILD_DEBUG "Build to allow debugging" OFF)
option(CUTFEM_BUILD_TESTS "Build and run library tests" OFF)
option(CUTFEM_BUILD_PYTHON_WRAPPER "Build a library for Python" OFF)
option(CUTFEM_BUILD_LIBRARY "Build the libraries" ON)
option(CUTFEM_BUILD_MAIN "Build the cpp example" ON)
option(CUTFEM_CREATE_DOCS "Build library documentation (requires Doxygen and Graphviz/Dot to be installed)" ON)

if(USE_MPI)
  set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/open-mpi/4.1.4_2/bin/mpicxx")
endif()






## ##############################
## Configuration
## ##############################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -std=c++2a -fPIC -O3")

message("compiler ${CMAKE_CXX_COMPILER}")

if(${CUTFEM_BUILD_DEBUG})
set(CMAKE_BUILD_TYPE Debug)
endif()


# if(${CUTFEM_CREATE_DOCS})
#     add_subdirectory(Documentation/)
# endif()


add_subdirectory(cpp/)

