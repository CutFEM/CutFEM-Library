

add_executable(force      force.cpp)
target_link_libraries(force PUBLIC ${EXTRA_LIBS})
target_include_directories(force PUBLIC ${CMAKE_BINARY_DIR})

add_executable(force_test      force_test.cpp )
target_link_libraries(force_test PUBLIC ${EXTRA_LIBS})
target_include_directories(force_test PUBLIC ${CMAKE_BINARY_DIR})

