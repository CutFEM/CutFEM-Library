cmake_minimum_required(VERSION 3.10)


set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/15.0.3/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/15.0.3/bin/clang++")


project(CutFEM VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


option(USE_MPI "Perform a parallel computation" OFF)
option(USE_OMP "Perform a parallel omp computation" OFF)
option(USE_UMFPACK "Found the UMFpack library" ON)
option(USE_LAPACK "Found the lapacke library" OFF)
option(USE_MUMPS "Found the MUMPS library" OFF)


if(USE_MPI)
  set(CMAKE_CXX_COMPILER "/usr/local/Cellar/open-mpi/4.1.4_2/bin/mpicxx")
endif()
## ##############################
## Configuration
## ##############################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -fPIC -std=c++17  -DCHECK_KN -O3 -w  -D_Float128=__float128")

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
#list(APPEND EXTRA_LIBS -lstdc++fs)
list(APPEND EXTRA_LIBS /usr/local/opt/libomp/lib/libomp.dylib)
list(APPEND EXTRA_INCLUDES /usr/local/opt/libomp/include)


get_directory_property(COMMON_SOURCES DIRECTORY common DEFINITION COMMON_SOURCES)
get_directory_property(FESPACE_SOURCES DIRECTORY FESpace DEFINITION FESPACE_SOURCES)
get_directory_property(PROBLEM_SOURCES DIRECTORY problem DEFINITION PROBLEM_SOURCES)
get_directory_property(SOLVER_SOURCES DIRECTORY solver DEFINITION SOLVER_SOURCES)


add_library(_darcy SHARED wrap_lib/lib_darcy.cpp)
# add_executable(_darcy mainFiles/darcy.cpp)

target_include_directories(_darcy
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")
target_sources(_darcy PRIVATE ${COMMON_SOURCES} ${FESPACE_SOURCES} ${SOLVER_SOURCES} ${PROBLEM_SOURCES})

#

set(EX_NAME limiter)
add_executable(${EX_NAME}  mainFiles/limiter_example.cpp )
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
#
set(EX_NAME showMesh)
add_executable(${EX_NAME}  mainFiles/showMesh.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})
#
#
set(EX_NAME shearBubble)
add_executable(${EX_NAME}  mainFiles/shearBubble.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})

set(EX_NAME stokesMixed)
add_executable(${EX_NAME}  mainFiles/stokesRT.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                        "${PROJECT_BINARY_DIR}"
                        ${EXTRA_INCLUDES})
