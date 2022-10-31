import numpy
import ctypes as ct
# from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

lib = ct.cdll.LoadLibrary('../cpp/build/lib/lib_darcy.dylib')

# Function types for divergence, velocity, pressure
#double fun_div(const R2 P, int compInd, int dom)
USER_FUNC = ct.CFUNCTYPE(ct.c_double, ct.POINTER(ct.c_double), ct.c_int, ct.c_int)

class Darcy2(object):

    def __init__(self):
        lib.Darcy2_new.argtypes               = None
        lib.Darcy2_new.restype                = ct.c_void_p
        self.obj = lib.Darcy2_new()

    def __del__(self):
        lib.Darcy2_delete.argtypes = [ct.c_void_p]
        lib.Darcy2_delete.restype = None
        lib.Darcy2_delete(self.obj)

    def build_mesh(self,nx,ny,orx,ory,lx,ly):
        lib.Darcy2_build_mesh.restypes         = None
        lib.Darcy2_build_mesh.argtypes         = [ct.c_void_p, ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_double, ct.c_double ]
        lib.Darcy2_build_mesh(self.obj, nx,ny,orx,ory,lx,ly)

    def init_space(self, s) :
        lib.Darcy2_init_space.restypes     = None
        lib.Darcy2_init_space.argtypes     = [ct.c_void_p, ct.c_char_p]
        lib.Darcy2_init_space(self.obj, s.encode('utf-8'))

    def add_bulk_integral(self, f) :
        lib.Darcy2_add_bulk_integral.restype    = None
        lib.Darcy2_add_bulk_integral.argtypes   = [ct.c_void_p, USER_FUNC]
        lib.Darcy2_add_bulk_integral(self.obj, f)


    def add_interface_integral(self, f) :
        lib.Darcy2_add_interface_integral.restype    = None
        lib.Darcy2_add_interface_integral.argtypes   = [ct.c_void_p, USER_FUNC]
        lib.Darcy2_add_interface_integral(self.obj, f)

    def add_natural_BC(self, f) :
        lib.Darcy2_add_natural_BC.restype    = None
        lib.Darcy2_add_natural_BC.argtypes   = [ct.c_void_p, USER_FUNC]
        lib.Darcy2_add_natural_BC(self.obj, f)

    def set_stabilization_penalty(self, cu, cp) :
        lib.Darcy2_set_stabilization_penalty.restype    = None
        lib.Darcy2_set_stabilization_penalty.argtypes   = [ct.c_void_p, ct.c_double, ct.c_double]
        lib.Darcy2_set_stabilization_penalty(self.obj, cu, cp)

    def add_full_stabilization(self) :
        lib.Darcy2_add_full_stabilization.restype    = None
        lib.Darcy2_add_full_stabilization.argtypes   = [ct.c_void_p]
        lib.Darcy2_add_full_stabilization(self.obj)

    def add_macro_stabilization(self, dlt) :
        lib.Darcy2_add_macro_stabilization.restype   = None
        lib.Darcy2_add_macro_stabilization.argtypes  = [ct.c_void_p, ct.c_double]
        lib.Darcy2_add_macro_stabilization(self.obj, dlt)

    # def solve_umfpack(self) :
    #     lib.Darcy2_solve_umfpack.restype    = None
    #     lib.Darcy2_solve_umfpack.argtypes   = [ct.c_void_p]
    #     lib.Darcy2_solve_umfpack(self.obj)

    def write_vtk_file(self, s):
        lib.Darcy2_write_vtk_file.restype    = None
        lib.Darcy2_write_vtk_file.argtypes   = [ct.c_void_p, ct.c_char_p]
        lib.Darcy2_write_vtk_file(self.obj, s.encode('utf-8'))

    def L2error_div(self, f):
        lib.Darcy2_L2error_div.restype    = ct.c_double
        lib.Darcy2_L2error_div.argtypes   = [ct.c_void_p, USER_FUNC]
        return lib.Darcy2_L2error_div(self.obj, f)

    def L2error_vel(self, f):
        lib.Darcy2_L2error_u.restype    = ct.c_double
        lib.Darcy2_L2error_u.argtypes   = [ct.c_void_p, USER_FUNC]
        return lib.Darcy2_L2error_u(self.obj, f)

    def L2error_pressure(self, f):
        lib.Darcy2_L2error_p.restype    = ct.c_double
        lib.Darcy2_L2error_p.argtypes   = [ct.c_void_p, USER_FUNC]
        return lib.Darcy2_L2error_p(self.obj, f)


    def solve_umfpack(self):
        lib.Darcy2_get_nz.restype  = ct.c_int
        lib.Darcy2_get_nz.argtypes = [ct.c_void_p]
        n_val = lib.Darcy2_get_nz(self.obj)
        lib.Darcy2_get_size_A.restype  = ct.c_int
        lib.Darcy2_get_size_A.argtypes = [ct.c_void_p]
        n_rhs = lib.Darcy2_get_size_A(self.obj)

        row  = numpy.empty([n_rhs+1], dtype=numpy.int32)
        col  = numpy.empty([n_val], dtype=numpy.int32)
        data = numpy.empty([n_val], dtype=numpy.float64)
        rhs  = numpy.empty([n_rhs], dtype=numpy.float64)

        lib.Darcy2_get_CSR_data.restype  = None
        lib.Darcy2_get_CSR_data.argtypes = [ct.c_void_p,
        numpy.ctypeslib.ndpointer(dtype=numpy.int32)  ,
        numpy.ctypeslib.ndpointer(dtype=numpy.int32)  ,
        numpy.ctypeslib.ndpointer(dtype=numpy.float64),
        numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        lib.Darcy2_get_CSR_data(self.obj, row, col, data, rhs)

        A = csr_matrix((data, col, row), shape=(n_rhs, n_rhs))
        x = spsolve(A, rhs, use_umfpack=True)

        lib.Darcy2_give_back_sol.restype  = None
        lib.Darcy2_give_back_sol.argtypes = [ct.c_void_p,numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        lib.Darcy2_give_back_sol(self.obj, x)


def set_verbose(s) :
    lib.set_verbose.restype    = None
    lib.set_verbose.argtypes   = [ct.c_int]
    lib.set_verbose(s)




#
