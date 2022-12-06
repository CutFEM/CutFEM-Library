from darcy_wrapper import *
from darcy_data_example1_2D import *


fun_level_set = USER_FUN_LS(func_level_set)
fun_div = USER_FUNC(func_div)
fun_neumann = USER_FUNC(func_neumann) 
fun_phat = USER_FUNC(func_phat)
fun_velocity = USER_FUNC(func_velocity)
fun_pressure = USER_FUNC(func_pressure)

set_verbose(0)

nx = 21

for x in range(3):

    darcy = Darcy2()

    darcy.build_mesh(nx, nx, sq_SW, sq_SW, sq_LGTH, sq_LGTH)
    darcy.init_space(fun_level_set, 'BDM1')

    darcy.add_bulk_integral(fun_div)
    darcy.add_interface_integral(fun_phat)
    darcy.add_natural_BC(fun_neumann)

    darcy.set_stabilization_penalty(0.1, 0.1)
    darcy.add_macro_stabilization(0.25)

    darcy.solve_umfpack()

    # darcy.write_vtk_file('../output/example.vtk')

    error_divu_L2 = darcy.L2error_div(fun_div)
    error_u_L2 = darcy.L2error_vel(fun_velocity)
    error_p_L2 = darcy.L2error_pressure(fun_pressure)

    print(error_divu_L2)
    print(error_p_L2)
    print(error_u_L2)

    nx = 2*nx - 1
