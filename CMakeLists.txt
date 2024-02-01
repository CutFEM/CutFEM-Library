cmake_minimum_required(VERSION 3.16)
cmake_minimum_required(VERSION 3.16)
project(libCutFEM VERSION 0.1 LANGUAGES CXX)

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
option(CUTFEM_BUILD_LIBRARY "Build the libraries" OFF)
option(CUTFEM_BUILD_MAIN "Build main file" OFF)
option(CUTFEM_BUILD_FORCE_PROJECT "Build force project" OFF)
option(CUTFEM_BUILD_EXAMPLE "Build the cpp example" ON)
option(CUTFEM_CREATE_DOCS "Build library documentation (requires Doxygen and Graphviz/Dot to be installed)" OFF)
option(CUTFEM_RAPID_YAML "Use rapid yaml" OFF)

configure_file(cutFEMConfig.in cutFEMConfig.h)


## ##############################
## Configuration
## ##############################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")



if(${CUTFEM_BUILD_DEBUG})
set(CMAKE_BUILD_TYPE Debug)
endif()


if(${CUTFEM_CREATE_DOCS})
    add_subdirectory(Documentation/)
endif()


add_subdirectory(cpp/)

