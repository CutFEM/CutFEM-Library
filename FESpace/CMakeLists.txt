
set(FESPACE_SOURCES
  QuadratureFormular.cpp
  P0.cpp
  P1.cpp
  P1dc.cpp
  P2.cpp
  P2dc.cpp
  P3dc.cpp
  P3.cpp
  RT0.cpp
  RT1.cpp
  BDM1.cpp
  finiteElement.cpp
  GTypeOfFE.cpp
  FESpace.cpp
  CutFESpace.cpp
  expression.cpp
  macroElement.cpp
  limiter.cpp
#  integrationFunFEM.cpp
#  funFEM.cpp
#  funCutFEM.cpp
#  integration.cpp
)



add_library(FESpace SHARED ${FESPACE_SOURCES})

target_include_directories(FESpace
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")
            ##               INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(FESpace PRIVATE common parallel)
