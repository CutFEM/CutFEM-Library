
find_package(PYTHON)

message(${CMAKE_SOURCE_DIR})

function(add_example FULL_NAME)
    get_filename_component(EXE_NAME ${FULL_NAME} NAME_WLE )
    add_executable(${EXE_NAME} ${FULL_NAME})
    target_link_libraries(${EXE_NAME} PRIVATE   ${CUTFEM_LIBS}  ${PYTHON_LIBRARIES})
    target_include_directories(${EXE_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                                                    ${CMAKE_SOURCE_DIR}
                                                    ${PYTHON_INCLUDES})
endfunction()

add_example(hyberbolic_conservation_law/DG_CutFEM_ex_1.cpp)
add_example(hyberbolic_conservation_law/DG_CutFEM_ex_2.cpp)

add_example(hyberbolic_conservation_law/limiter_smooth_solution.cpp)
add_example(hyberbolic_conservation_law/limiter_accuracy_conservation_test.cpp)
add_example(hyberbolic_conservation_law/limiter_non_smooth_initial_data.cpp)


add_example(darcy_fictitious/DF_example1.cpp)


add_example(darcy_interface/DI_example1.cpp)
add_example(darcy_interface/DI_example1_3D.cpp)
add_example(darcy_interface/DI_example2.cpp)


add_example(convection_diffusion/bulk.cpp)
add_example(convection_diffusion/coupled.cpp)



if(PYTHON_FOUND)
add_example(stokes/fictitious_problem/stokes_fictitious_example1.cpp)
add_example(stokes/fictitious_problem/stokes_fictitious_vorticity_example1.cpp)
add_example(stokes/interface_problem/stokesRT_static_drop.cpp)
add_example(stokes/interface_problem/stokesRT_dynamic_drop.cpp)
endif()

                                                                                



