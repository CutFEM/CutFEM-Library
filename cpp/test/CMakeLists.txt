set (tests
        ./test_main.cpp
        )        

add_subdirectory(Catch2)
add_executable(runTests ${tests})
target_link_libraries(runTests PRIVATE Catch2::Catch2WithMain ${EXTRA_LIBS})
target_include_directories(runTests PUBLIC
"${CMAKE_CURRENT_SOURCE_DIR}"
"${PROJECT_BINARY_DIR}")

