cmake_minimum_required(VERSION 3.10)

project(CutFEM VERSION 1.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

## ##############################
## Options
## ##############################

option(USE_MPI "Perform a parallel computation" ON)
option(USE_OMP "Perform a parallel omp computation" OFF)
option(USE_UMFPACK "Found the UMFpack library" ON)
option(USE_LAPACK "Found the lapacke library" OFF)
option(USE_MUMPS "Found the MUMPS library" ON)

option(CUTFEM_BUILD_TESTS "Build and run library tests" ON)
option(CUTFEM_BUILD_PYTHON_WRAPPER "Build a library for Python" ON)
option(CUTFEM_BUILD_MAIN "Build the cpp example" ON)

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

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -fPIC  -std=c++17 -O3 ")


# To allow debugging
set(CMAKE_BUILD_TYPE Debug)

if(USE_MPI)
  find_package(MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  add_subdirectory(parallel)
  list(APPEND EXTRA_LIBS parallel)
endif()
add_subdirectory(common)
list(APPEND EXTRA_LIBS common)
add_subdirectory(FESpace)
list(APPEND EXTRA_LIBS FESpace)
add_subdirectory(solver)
list(APPEND EXTRA_LIBS solver)
add_subdirectory(problem)
list(APPEND EXTRA_LIBS problem)


# if(${CUTFEM_BUILD_TESTS})
# add_subdirectory(test/)
# endif()


if(${CUTFEM_BUILD_PYTHON_WRAPPER})
add_library(_darcy SHARED wrap_lib/lib_darcy.cpp)
target_link_libraries(_darcy PUBLIC ${EXTRA_LIBS})
endif()

#message("${EXTRA_INCLUDES}")
#message("${UMFPACK_INCLUDES}")

if(${CUTFEM_BUILD_MAIN})


set(EX_NAME surfactant)
add_executable(${EX_NAME}  mainFiles/surfactant_cg.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

set(EX_NAME mesh)
add_executable(${EX_NAME}  mainFiles/test_mesh.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

set(EX_NAME darcy)
add_executable(${EX_NAME}  mainFiles/darcy.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

set(EX_NAME bulkCD)
add_executable(${EX_NAME}  mainFiles/bulk_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME cgbulkCD)
add_executable(${EX_NAME}  mainFiles/cg_bulk_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME cgcoupledCD)
add_executable(${EX_NAME}  mainFiles/cg_coupled_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )


target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#                        
set(EX_NAME limiter)
add_executable(${EX_NAME}  mainFiles/limiter_example.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})


# set(EX_NAME showMesh)
# add_executable(${EX_NAME}  mainFiles/showMesh.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

# set(EX_NAME shearBubble)
# add_executable(${EX_NAME}  mainFiles/shearBubble.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

# set(EX_NAME stokesMixed)
# add_executable(${EX_NAME}  mainFiles/stokesRT.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})



endif()