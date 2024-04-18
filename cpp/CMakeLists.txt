
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -std=c++2a -fPIC -O3")

if(USE_MPI)
  find_package(MPI REQUIRED)
  
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  include_directories(${MPI_CXX_INCLUDE_PATH})

  add_subdirectory(parallel)
  list(APPEND CUTFEM_LIBS parallel)
endif()

add_subdirectory(common)
list(APPEND CUTFEM_LIBS common)
add_subdirectory(FESpace)
list(APPEND CUTFEM_LIBS FESpace)
add_subdirectory(solver)
list(APPEND CUTFEM_LIBS solver)
add_subdirectory(problem)
list(APPEND CUTFEM_LIBS problem)




if(${CUTFEM_BUILD_EXAMPLE})
add_subdirectory(example/)
endif()

