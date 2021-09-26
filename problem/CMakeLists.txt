

set(PROBLEM_SOURCES
  levelSet.cpp
  reinitialization.cpp
  #cutLevelSet.cpp
  GenericMapping.cpp
  CutFEM_parameter.cpp
  itemVF.cpp
  testFunction.cpp
)

add_library(problem SHARED ${PROBLEM_SOURCES})


target_include_directories(problem
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(problem PRIVATE parallel common FESpace solver)
