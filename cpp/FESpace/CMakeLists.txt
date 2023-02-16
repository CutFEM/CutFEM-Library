
set(FESPACE_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/QuadratureFormular.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P0.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P1.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P1dc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P2.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P2dc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P3dc.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P3.cpp
${CMAKE_CURRENT_SOURCE_DIR}/P4.cpp
${CMAKE_CURRENT_SOURCE_DIR}/RT0.cpp
  #RT1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RT2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BDM1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BDM2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/P2BR.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/transformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/finiteElement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/GTypeOfFE.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FESpace.cpp
  #CutFESpace.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/expression.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/macroElement.cpp
  #limiter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/paraview.cpp
)


add_library(FESpace SHARED ${FESPACE_SOURCES})
target_include_directories(FESpace
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")
if(USE_MPI)                   
target_link_libraries(FESpace PRIVATE common parallel)
else()
target_link_libraries(FESpace PRIVATE common)
endif()

