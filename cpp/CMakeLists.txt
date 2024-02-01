
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -std=c++20 -fPIC -O3")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -std=c++2a -fPIC -O3")

if(USE_MPI)
  find_package(MPI REQUIRED)
  
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
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


if(${CUTFEM_BUILD_PYTHON_WRAPPER})
add_library(_pyCutFEM SHARED 
                            wrap_lib/lib_cutfem.cpp
                            wrap_lib/lib_darcy.cpp
                            wrap_lib/lib_stokes.cpp)
target_link_libraries(_pyCutFEM PUBLIC ${EXTRA_LIBS})
endif()


if(${CUTFEM_BUILD_FORCE_PROJECT})
add_subdirectory(force_project)
# list(APPEND EXTRA_LIBS force_project)
endif()



if(${CUTFEM_BUILD_EXAMPLE})
add_subdirectory(example/)
endif()


if(${CUTFEM_BUILD_TESTS})
add_subdirectory(test/)
endif()

if(${CUTFEM_BUILD_MAIN})



# set(EX_NAME risingbubble)
# add_executable(${EX_NAME}  mainFiles/rising_bubble.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME stationaryvortex)
add_executable(${EX_NAME}  mainFiles/stationary_vortex.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME taylorgreen)
add_executable(${EX_NAME}  mainFiles/taylor_green_vortex.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME henry)
add_executable(${EX_NAME}  mainFiles/algoim_henry.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME bulkCD2)
add_executable(${EX_NAME}  mainFiles/bulk_convection_diffusion2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )




set(EX_NAME kiln_bed)
add_executable(${EX_NAME}  mainFiles/kiln_bed_motion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME maxwell3d)
add_executable(${EX_NAME}  mainFiles/maxwell3D.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME statbulk)
add_executable(${EX_NAME}  mainFiles/stationary_bulk.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME paper)
add_executable(${EX_NAME}  mainFiles/paper.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME poisson)
add_executable(${EX_NAME}  mainFiles/poissonFEM.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME darcy)
add_executable(${EX_NAME}  mainFiles/darcy.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

# set(EX_NAME stokes)
# add_executable(${EX_NAME}  mainFiles/stokes_fictitious_example1.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME algoimbulk2)
add_executable(${EX_NAME}  mainFiles/algoim_bulk2.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME bulkCD)
add_executable(${EX_NAME}  mainFiles/bulk_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )
set(EX_NAME dgbulkCD)
add_executable(${EX_NAME}  mainFiles/dg_bulk_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME cgbulkCD)
add_executable(${EX_NAME}  mainFiles/cg_bulk_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )

set(EX_NAME cgcoupledCD)
add_executable(${EX_NAME}  mainFiles/cg_coupled_convection_diffusion.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS} )
# #                        
set(EX_NAME limiter)
add_executable(${EX_NAME}  mainFiles/limiter_example.cpp )
target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})


# set(EX_NAME surfactant)
# add_executable(${EX_NAME}  mainFiles/surfactant.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

# set(EX_NAME surfactantCG2)
# add_executable(${EX_NAME}  mainFiles/surfactantCG2.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})


# set(EX_NAME stokes)
# add_executable(${EX_NAME}  mainFiles/stokes_interface_problem.cpp )
# target_link_libraries(${EX_NAME} PUBLIC ${EXTRA_LIBS})

set(EX_NAME test_BDM2)
add_executable(${EX_NAME}  mainFiles/test_BDM2.cpp )
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
