
find_package(UMFPACK)
find_package(LAPACK)
find_package(MUMPS)


if(UMFPACK_LIBRARIES)
    include_directories(${UMFPACK_INCLUDES})
    list(APPEND LINK_LIBRARY ${UMFPACK_LIBRARIES})
    option(USE_UMFPACK "Found the UMFpack library" ON)
    message("Found UMFPACK")
else()
    message("Not Found UMFPACK")
endif()

if(LAPACK_LIBRARIES)
  include_directories(${LAPACK_INCLUDES})
  list(APPEND LINK_LIBRARY ${LAPACK_LIBRARIES})
#  option(USE_LAPACK "Found the lapacke library" ON)
  message("Found LAPACK")
else()
      message("Not Found LAPACK")
endif()

if(MUMPS_LIBRARIES AND USE_MPI)
  include_directories(${MUMPS_INCLUDES})
  list(APPEND LINK_LIBRARY ${MUMPS_LIBRARIES})
  option(USE_MUMPS "Found the mumps library" ON)
  message("Found MUMPS")
else()
  if(USE_MPI)
    message("Not Found MUMPS")
  elseif(MUMPS_LIBRARIES)
    include_directories(${MUMPS_INCLUDES})
    list(APPEND LINK_LIBRARY ${MUMPS_LIBRARIES})
    option(USE_MUMPS "Found the mumps library" ON)
    message("Found seq MUMPS")
  else()
    message("Not Found MUMPS")
  endif()
endif()

configure_file(../cutFEMConfig.in ../cutFEMConfig.h)

set(SOLVER_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp
)

if(USE_MUMPS AND USE_MPI)
list(APPEND  SOLVER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solverMumpsMPI.cpp)
list(APPEND LINK_LIBRARY parallel)
elseif(USE_MUMPS)
list(APPEND  SOLVER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solverMumps.cpp)
endif()

add_library(solver SHARED ${SOLVER_SOURCES})

target_include_directories(solver
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")


target_link_libraries(solver PRIVATE ${LINK_LIBRARY})
