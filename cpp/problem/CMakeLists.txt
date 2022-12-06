

set(PROBLEM_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/CutFEM_parameter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/testFunction.cpp
${CMAKE_CURRENT_SOURCE_DIR}/itemVF.cpp
${CMAKE_CURRENT_SOURCE_DIR}/problem.cpp
${CMAKE_CURRENT_SOURCE_DIR}/baseProblem.cpp
  #GenericMapping.cpp
  # ${CMAKE_CURRENT_SOURCE_DIR}/extension.cpp
  #  levelSet.cpp
  #  reinitialization.cpp

)

add_library(problem SHARED ${PROBLEM_SOURCES})

target_include_directories(problem 
INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if(USE_MPI)
target_link_libraries(problem PRIVATE parallel common FESpace solver )
else()
target_link_libraries(problem PRIVATE common FESpace solver)
endif()
#