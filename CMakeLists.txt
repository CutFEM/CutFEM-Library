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


set(EX_NAME showMesh)
add_executable(${EX_NAME}  mainFiles/showMesh.cpp )
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
set(EX_NAME dropChanel)
add_executable(${EX_NAME}  mainFiles/dropChanel.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME reinitialization2)
add_executable(${EX_NAME}  mainFiles/reinitialization2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
set(EX_NAME reinitialization3)
add_executable(${EX_NAME}  mainFiles/reinitialization3.cpp )
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
set(EX_NAME RBeuler2)
add_executable(${EX_NAME}  mainFiles/RBeuler2D.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
set(EX_NAME RBeuler3)
add_executable(${EX_NAME}  mainFiles/RBeuler3D.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
set(EX_NAME RBSeuler3)
add_executable(${EX_NAME}  mainFiles/RBSeuler3D.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME RBsimpsonP1)
add_executable(${EX_NAME}  mainFiles/RBsimpsonP1.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME RBsurfactantSimpson)
add_executable(${EX_NAME}  mainFiles/mainRBsurfactantSimpsonP1.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
set(EX_NAME RBSConservative2)
add_executable(${EX_NAME}  mainFiles/RBSConservative2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
#
set(EX_NAME pairBubble2)
add_executable(${EX_NAME}  mainFiles/pairBubble.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
set(EX_NAME shearBubble)
add_executable(${EX_NAME}  mainFiles/shearBubbleSurfactantConservative.cpp )
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
set(EX_NAME bulkSurfactantNavierStokes)
add_executable(${EX_NAME}  mainFiles/bulkSurfactantNavierStokes.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME testBDM1)
add_executable(${EX_NAME}  mainFiles/mainBDM1.cpp )
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

set(EX_NAME surfactant2D)
add_executable(${EX_NAME}  mainFiles/surfactant.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME curvature3)
add_executable(${EX_NAME}  mainFiles/curvature3.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME surfactant3D)
add_executable(${EX_NAME}  mainFiles/surfactant3D.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                           ${EXTRA_INCLUDES})

set(EX_NAME poisson)
add_executable(${EX_NAME}  mainFiles/poisson.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                         ${EXTRA_INCLUDES})
#
set(EX_NAME poisson3)
add_executable(${EX_NAME}  mainFiles/poisson3.cpp )
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

#set(EX_NAME surfactant)
#add_executable(${EX_NAME}  mainFiles/surfactant.cpp )
#target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#target_include_directories(${EX_NAME} PUBLIC
#                        "${PROJECT_BINARY_DIR}"
#                        ${EXTRA_INCLUDES})


#set(EX_NAME CutFEM_Paper2)
#add_executable(${EX_NAME}  mainFiles/mainTimeCutFEM.cpp )
#target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#target_include_directories(${EX_NAME} PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          ${EXTRA_INCLUDES})

set(EX_NAME iceSheet)
add_executable(${EX_NAME}  mainFiles/iceSheetMain.cpp)
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})



set(EX_NAME stokes2)
add_executable(${EX_NAME}  mainFiles/stokes.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})
#
set(EX_NAME stokes3)
add_executable(${EX_NAME}  mainFiles/stokes3.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

set(EX_NAME curvature2)
add_executable(${EX_NAME}  mainFiles/curvature2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES})

#set(EX_NAME FEM_Stokes)
#add_executable(${EX_NAME}  mainFiles/mainStokes.cpp )
#target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#target_include_directories(${EX_NAME} PUBLIC
#                         "${PROJECT_BINARY_DIR}"
#                          ${EXTRA_INCLUDES})

set(EX_NAME cutStokes2)
add_executable(${EX_NAME}  mainFiles/cutStokes2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            ${EXTRA_INCLUDES})

#
set(EX_NAME cutStokes3)
add_executable(${EX_NAME}  mainFiles/cutStokes3.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${EX_NAME} PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            ${EXTRA_INCLUDES})

#set(EX_NAME FEM_Erik)
#add_executable(${EX_NAME}  mainFiles/erikMain.cpp )
#target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#target_include_directories(${EX_NAME} PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          ${EXTRA_INCLUDES})

#set(EX_NAME FEM_Time)
#add_executable(${EX_NAME}  mainFiles/mainTimeFEM.cpp )
#target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})
#target_include_directories(${EX_NAME} PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          ${EXTRA_INCLUDES})
