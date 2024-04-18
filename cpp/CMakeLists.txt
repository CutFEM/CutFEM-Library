
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
