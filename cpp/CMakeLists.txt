
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_SILENT "-Qunused-arguments")


if(USE_OMP)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_SILENT} ${LDFLAGS}  -L/opt/homebrew/opt/libomp/lib/ -fopenmp -std=c++2a -fPIC -O3")
include_directories(/opt/homebrew/opt/libomp/include)
else()
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}   -std=c++2a -fPIC -O3")
endif()


if(USE_MPI)
  find_package(MPI )
  include_directories(${MPI_CXX_INCLUDE_DIR})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

add_subdirectory(parallel)
list(APPEND CUTFEM_LIBS parallel)
add_subdirectory(common)
list(APPEND CUTFEM_LIBS common)
add_subdirectory(FESpace)
list(APPEND CUTFEM_LIBS FESpace)
add_subdirectory(solver)
list(APPEND CUTFEM_LIBS solver)
add_subdirectory(problem)
list(APPEND CUTFEM_LIBS problem)


add_library(cutfem INTERFACE)
target_link_libraries(cutfem INTERFACE
  parallel
  common
  FESpace
  solver
  problem
)

if(${CUTFEM_BUILD_EXAMPLE})
add_subdirectory(example/)
endif()

