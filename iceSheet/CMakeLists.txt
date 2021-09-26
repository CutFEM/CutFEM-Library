
set(ICESHEET_SOURCES
  markerIceSheet.cpp
  iceSheet.cpp
)

add_library(iceSheet SHARED ${ICESHEET_SOURCES})


target_include_directories(iceSheet
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(iceSheet PRIVATE parallel common FESpace solver)
