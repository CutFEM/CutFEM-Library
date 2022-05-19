cmake_minimum_required(VERSION 3.10)


set(CMAKE_C_COMPILER "/usr/local/Cellar/gcc/11.1.0/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/11.1.0/bin/g++-11")

project(CutFEM VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


option(USE_MPI "Perform a parallel computation" ON)
option(USE_OMP "Perform a parallel omp computation" ON)
option(USE_UMFPACK "Found the UMFpack library" ON)
option(USE_LAPACK "Found the lapacke library" OFF)
option(USE_MUMPS "Found the MUMPS library" ON)


if(USE_MPI)
  set(CMAKE_CXX_COMPILER "/usr/local/Cellar/open-mpi/4.1.1_2/bin/mpicxx")
endif()
## ##############################
## Configuration
## ##############################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -fPIC -std=c++17  -DCHECK_KN -O3 -w -fopenmp -D_Float128=__float128")

message( " cxx compiler ${CMAKE_CXX_COMPILER}")


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
list(APPEND EXTRA_LIBS -lstdc++fs)


#
set(EX_NAME stokesMixed)
add_executable(${EX_NAME}  mainFiles/stokesRT.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})
#
set(EX_NAME peiFu)
add_executable(${EX_NAME}  mainFiles/peiFu.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})

set(EX_NAME darcy)
add_executable(${EX_NAME}  mainFiles/darcy.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                       "${PROJECT_BINARY_DIR}"
                       ${EXTRA_INCLUDES})

#
set(EX_NAME convDiff)
add_executable(${EX_NAME}  mainFiles/statConvectionDiffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})


set(EX_NAME testMarker)
add_executable(${EX_NAME}  mainFiles/testMarker.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})



#
set(EX_NAME bulkSurface)
add_executable(${EX_NAME}  mainFiles/bulkSurface.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})
#
set(EX_NAME bulkSurface3)
add_executable(${EX_NAME}  mainFiles/bulkSurface3.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})

set(EX_NAME mainTest)
add_executable(${EX_NAME}  mainFiles/mainTest.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            ${EXTRA_INCLUDES})

set(EX_NAME surfactant2)
add_executable(${EX_NAME}  mainFiles/surfactant.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})

#
set(EX_NAME poisson)
add_executable(${EX_NAME}  mainFiles/poisson.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})
