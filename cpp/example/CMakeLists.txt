
#set(INCLUDE_PYTHON /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/include/python3.10/
#/usr/local/lib/python3.10/site-packages/numpy/core/include)
#set(LIB_PYTHON /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/lib/libpython3.10.dylib) 

find_package(PYTHON)

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

#Darcy fictitious example                        
add_executable(darcy_ex_1 darcy_fictitious/darcy_example1.cpp)
target_link_libraries( darcy_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )


add_executable(darcy_interface_1 darcy_interface/darcy_example1.cpp)
target_link_libraries( darcy_interface_1
                        PRIVATE 
                        ${EXTRA_LIBS} )    
                        
add_executable(darcy_interface3D_1 darcy_interface/darcy_example1_3D.cpp)
target_link_libraries( darcy_interface3D_1
                        PRIVATE 
                        ${EXTRA_LIBS} )                          

add_executable(darcy_interface_2 darcy_interface/darcy_example2.cpp)
target_link_libraries( darcy_interface_2
                        PRIVATE 
                        ${EXTRA_LIBS} )   

# Stokes fictitious examples                        
add_executable(stokesRT_ex_1 stokes/fictitious_problem/stokes_fictitious_example1.cpp)
target_link_libraries( stokesRT_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )

add_executable(stokesV_ex_1 stokes/fictitious_problem/stokes_fictitious_vorticity_example1.cpp)
target_link_libraries( stokesV_ex_1
                        PRIVATE 
                        ${EXTRA_LIBS} )
                        
# Stokes interface examples                 
add_executable(stokesRT_static_drop stokes/interface_problem/stokesRT_static_drop.cpp)
target_link_libraries( stokesRT_static_drop
                        PRIVATE 
                        ${EXTRA_LIBS}
                        ${PYTHON_LIBRARIES} )                        
target_include_directories(stokesRT_static_drop PRIVATE 
                                ${CMAKE_CURRENT_BINARY_DIR}
                                ${PYTHON_INCLUDES})     
                                
add_executable(stokesRT_dynamic_drop stokes/interface_problem/stokesRT_dynamic_drop.cpp)
target_link_libraries( stokesRT_dynamic_drop
                        PRIVATE 
                        ${EXTRA_LIBS}
                        ${PYTHON_LIBRARIES} )                        
target_include_directories(stokesRT_dynamic_drop PRIVATE 
                                ${CMAKE_CURRENT_BINARY_DIR}
                                ${PYTHON_INCLUDES})       
                                                                                          

# Navier Stokes examples
# add_executable(NS_raisingDrop2D navier_stokes/navier_stokes_raising_drop_2D.cpp)
# target_link_libraries(NS_raisingDrop2D
#                         PRIVATE 
#                         ${EXTRA_LIBS} 
#                         )
                        

# add_executable(NSS_ex_2 navier_stokes/navier_stokes_surfactant_example2.cpp)
# target_link_libraries(NSS_ex_2
#                         PRIVATE 
#                         ${EXTRA_LIBS} )


