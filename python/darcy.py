# from math import *
import numpy as np
import ctypes
from darcy_wrapper import *
from matplotlib import pyplot as plt
from darcy_data import *
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


fun_div         = USER_FUNC(func_div)
fun_neumann     = USER_FUNC(func_neumann)
fun_phat        = USER_FUNC(func_phat)
fun_velocity    = USER_FUNC(func_velocity)
fun_pressure    = USER_FUNC(func_pressure)

set_verbose(0)

darcy = Darcy2()

darcy.build_mesh(41,41,0.,0.,1.,1.)
darcy.init_space('RT0')

darcy.add_bulk_integral(fun_div)
darcy.add_interface_integral(fun_phat)
darcy.add_natural_BC(fun_neumann)

darcy.set_stabilization_penalty(0.1, 0.1)
darcy.add_macro_stabilization(0.25)

darcy.solve_umfpack()
#
# darcy.write_vtk_file('../output/scotti_nx21.vtk')
#
error_divu_L2 = darcy.L2error_div(fun_div)
error_u_L2    = darcy.L2error_vel(fun_velocity)
error_p_L2    = darcy.L2error_pressure(fun_pressure)


print(error_divu_L2)
print(error_p_L2)
print(error_u_L2)
