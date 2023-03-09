

add_executable(darcy_ex_1 darcy/darcy_example1.cpp)
target_link_libraries( darcy_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )

# add_executable(stokes_ex_1 stokes/stokes_fictitious_example1.cpp)
# target_link_libraries( stokes_ex_1
#                         PRIVATE 
#                         ${EXTRA_LIBS} )

add_executable(NS_raisingDrop2D navier_stokes/navier_stokes_raising_drop_2D.cpp)
target_link_libraries(NS_raisingDrop2D
                        PRIVATE 
                        ${EXTRA_LIBS} )

# add_executable(NSS_ex_2 navier_stokes/navier_stokes_surfactant_example2.cpp)
# target_link_libraries(NSS_ex_2
#                         PRIVATE 
#                         ${EXTRA_LIBS} )


