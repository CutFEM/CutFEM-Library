
set(COMMON_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/global.cpp
${CMAKE_CURRENT_SOURCE_DIR}/GTime.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dataStruct1D.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dataStruct2D.cpp
${CMAKE_CURRENT_SOURCE_DIR}/dataStruct3D.cpp
${CMAKE_CURRENT_SOURCE_DIR}/GenericMesh.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Mesh1dn.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Mesh2dn.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Mesh3dn.cpp
${CMAKE_CURRENT_SOURCE_DIR}/base_interface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/interface_levelSet.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cut_method.cpp
${CMAKE_CURRENT_SOURCE_DIR}/marker.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parametrization.hpp
${CMAKE_CURRENT_SOURCE_DIR}/spline.cpp
${CMAKE_CURRENT_SOURCE_DIR}/SparseMatMap.cpp
${CMAKE_CURRENT_SOURCE_DIR}/geometry.cpp
)
# message("common source files = ${COMMON_SOURCES}")
add_library(common SHARED ${COMMON_SOURCES})

target_include_directories(common
PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")
                      ##     INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
if(USE_MPI)                   
target_link_libraries(common PRIVATE parallel)
else()
target_link_libraries(common)
endif()

# target_sources(_darcy PRIVATE more_cool_stuff.cpp)


