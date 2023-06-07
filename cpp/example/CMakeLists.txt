add_executable(limiter_paper1_ex_2 hyberbolic_conservation_law/DG_CutFEM_ex_2.cpp)
target_link_libraries( limiter_paper1_ex_2
                        PRIVATE 
                        ${EXTRA_LIBS} )


add_executable(limiter_ex_1 hyberbolic_conservation_law/limiter_smooth_solution.cpp)
target_link_libraries( limiter_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )

add_executable(limiter_ex_2 hyberbolic_conservation_law/limiter_accuracy_conservation_test.cpp)
target_link_libraries( limiter_ex_2
                        PRIVATE 
                        ${EXTRA_LIBS} )
                    
add_executable(limiter_ex_3 hyberbolic_conservation_law/limiter_non_smooth_initial_data.cpp)
target_link_libraries( limiter_ex_3
                        PRIVATE 
                        ${EXTRA_LIBS} )

add_executable(darcy_ex_1 darcy_fictitious/darcy_example1.cpp)
target_link_libraries( darcy_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )

add_executable(stokesRT_ex_1 stokes/stokes_fictitious_example1.cpp)
target_link_libraries( stokesRT_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )

add_executable(stokesV_ex_1 stokes/stokes_fictitious_vorticity_example1.cpp)
target_link_libraries( stokesV_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )                        

add_executable(NS_raisingDrop2D navier_stokes/navier_stokes_raising_drop_2D.cpp)
target_link_libraries(NS_raisingDrop2D
                        PRIVATE 
                        ${EXTRA_LIBS} )

# add_executable(NSS_ex_2 navier_stokes/navier_stokes_surfactant_example2.cpp)
# target_link_libraries(NSS_ex_2
#                         PRIVATE 
#                         ${EXTRA_LIBS} )


