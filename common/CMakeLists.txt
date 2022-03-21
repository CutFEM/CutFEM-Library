
set(COMMON_SOURCES
  libmesh5.c
  GTime.cpp
  dataStruct1D.cpp
  dataStruct2D.cpp
  dataStruct3D.cpp
  GenericMesh.cpp
  Mesh1dn.cpp
  Mesh2dn.cpp
  Mesh3dn.cpp
  interface_levelSet.cpp
  cut_mesh.cpp
  cut_method.cpp
  marker.cpp
  base_interface.cpp
  Interface2dn.cpp
  Interface3dn.cpp
  SparseMatMap.cpp
  geometry.cpp
  spline.cpp

)

add_library(common SHARED ${COMMON_SOURCES})

target_include_directories(common
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")
                      ##     INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(common PRIVATE parallel)
