import sys
import os

# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
# where the current directory is present.
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)


from stokes_wrapper import *
from stokes_data_example1_2D import *
import matplotlib.pyplot as plt
import numpy as np


fun_level_set = USER_FUN_LS(func_level_set)
fun_rhs = USER_FUNC(func_rhs)
fun_div = USER_FUNC(func_div)
fun_velocity = USER_FUNC(func_velocity)
fun_pressure = USER_FUNC(func_pressure)



set_verbose(2)

nx = 11
h = np.empty(0)
err_u = np.empty(0)
err_p = np.empty(0)
err_div = np.empty(0)

stab_manual = -1
stab_classic = 0
stab_mixed = 1

element = 'BDM1'

for x in range(4):
    stokes = FictitiousStokesCutFEM()

    stokes.build_mesh(nx, nx, sq_SW, sq_SW, sq_LGTH, sq_LGTH)
    stokes.init_space(fun_level_set, element)

    stokes.add_bulk_integral(fun_rhs)
    stokes.add_interface_integral(fun_velocity)
    
    stokes.add_lagrange_multiplier_classic(0.)

    stokes.set_stabilization_Cu(1.)
    stokes.set_stabilization_Cp(0.)
    stokes.set_stabilization_Cpu(1.)
    
    stokes.add_macro_stabilization(0.5, stab_mixed)

    stokes.solve()
    
    stokes.post_process_pressure(fun_pressure)

    stokes.write_vtk_file('python/output/example_stokesRT_'+str(x)+'.vtk')
    
    error_divu_L2 = stokes.L2error_div(fun_div)
    error_u_L2 = stokes.L2error_vel(fun_velocity)
    error_p_L2 = stokes.L2error_pressure(fun_pressure)

    h = np.append(h, [1./nx])
    err_u = np.append(err_u, error_u_L2)
    err_p = np.append(err_p, error_p_L2)
    err_div = np.append(err_div, error_divu_L2)

    nx = 2*nx - 1
    
    

print(h)
print(err_p)
print(err_u)
print(err_div)

# plot
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(h, err_p, 'r*--')
ax2.plot(h, err_u, 'b^--')
ax1.plot(h, h, 'k-.', h, 3*h**2, 'k--')
ax2.plot(h, h, 'k-.', h, 3*h**2, 'k--')

for ax in fig.get_axes():
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)

fig.suptitle('Stokes Vorticity - Example 1 -'+element)
plt.show()
