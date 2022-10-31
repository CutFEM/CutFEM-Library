

set(PROBLEM_SOURCES
  CutFEM_parameter.cpp
  testFunction.cpp
  itemVF.cpp
  problem.cpp
  baseProblem.cpp
  #GenericMapping.cpp
  extension.cpp
  #  levelSet.cpp
  #  reinitialization.cpp

)

add_library(problem SHARED ${PROBLEM_SOURCES})

target_include_directories(problem
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(problem PRIVATE parallel common FESpace solver /usr/local/opt/libomp/lib/libomp.dylib)
