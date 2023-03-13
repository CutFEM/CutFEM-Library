import numpy
import ctypes as ct
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

lib = ct.cdll.LoadLibrary('build/lib/lib_pyCutFEM.dylib')

USER_FUNC = ct.CFUNCTYPE(ct.c_double, ct.POINTER(
    ct.c_double), ct.c_int, ct.c_int)

USER_FUN_LS = ct.CFUNCTYPE(ct.c_double, ct.POINTER(
    ct.c_double))


def set_verbose(s):
    lib.set_verbose.restype = None
    lib.set_verbose.argtypes = [ct.c_int]
    lib.set_verbose(s)


class ProblemCutFEM(object):

    def build_mesh(self, nx, ny, orx, ory, lx, ly):
        lib.build_mesh.restypes = None
        lib.build_mesh.argtypes = [
            ct.c_void_p, ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_double, ct.c_double]
        lib.build_mesh(self.obj, nx, ny, orx, ory, lx, ly)

    def init_space(self, f, s):
        lib.init_space.restypes = None
        lib.init_space.argtypes = [
            ct.c_void_p, USER_FUN_LS, ct.c_char_p]
        lib.init_space(self.obj, f, s.encode('utf-8'))

    def add_bulk_integral(self, f):
        lib.add_bulk_integral.restype = None
        lib.add_bulk_integral.argtypes = [ct.c_void_p, USER_FUNC]
        lib.add_bulk_integral(self.obj, f)

    def add_interface_integral(self, f):
        lib.add_interface_integral.restype = None
        lib.add_interface_integral.argtypes = [ct.c_void_p, USER_FUNC]
        lib.add_interface_integral(self.obj, f)

    def add_lagrange_multiplier(self):
        lib.add_lagrange_multiplier.argtypes = [ct.c_void_p]
        lib.add_lagrange_multiplier.restype = None
        lib.add_lagrange_multiplier(self.obj)

    def set_stabilization_penalty(self, cu, cp):
        lib.set_stabilization_penalty.restype = None
        lib.set_stabilization_penalty.argtypes = [
            ct.c_void_p, ct.c_double, ct.c_double]
        lib.set_stabilization_penalty(self.obj, cu, cp)

    def add_full_stabilization(self, stab):
        lib.add_full_stabilization.restype = None
        lib.add_full_stabilization.argtypes = [ct.c_void_p, ct.c_int]
        lib.add_full_stabilization(self.obj, stab)

    def add_macro_stabilization(self, dlt, stab):
        lib.add_macro_stabilization.restype = None
        lib.add_macro_stabilization.argtypes = [
            ct.c_void_p, ct.c_double, ct.c_int]
        lib.add_macro_stabilization(self.obj, dlt, stab)

    def post_process_pressure(self, f):
        lib.post_processing_pressure.restype = None
        lib.post_processing_pressure.argtypes = [ct.c_void_p, USER_FUNC]
        lib.post_processing_pressure(self.obj, f)

    def write_vtk_file(self, s):
        lib.write_vtk_file.restype = None
        lib.write_vtk_file.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.write_vtk_file(self.obj, s.encode('utf-8'))

    def L2error_div(self, f):
        lib.L2error_div.restype = ct.c_double
        lib.L2error_div.argtypes = [ct.c_void_p, USER_FUNC]
        return lib.L2error_div(self.obj, f)

    def L2error_vel(self, f):
        lib.L2error_u.restype = ct.c_double
        lib.L2error_u.argtypes = [ct.c_void_p, USER_FUNC]
        return lib.L2error_u(self.obj, f)

    def L2error_pressure(self, f):
        lib.L2error_p.restype = ct.c_double
        lib.L2error_p.argtypes = [ct.c_void_p, USER_FUNC]
        return lib.L2error_p(self.obj, f)

    def solve(self):
        lib.get_nz.restype = ct.c_int
        lib.get_nz.argtypes = [ct.c_void_p]
        n_val = lib.get_nz(self.obj)
        lib.get_size_A.restype = ct.c_int
        lib.get_size_A.argtypes = [ct.c_void_p]
        n_rhs = lib.get_size_A(self.obj)

        row = numpy.empty([n_rhs+1], dtype=numpy.int32)
        col = numpy.empty([n_val], dtype=numpy.int32)
        data = numpy.empty([n_val], dtype=numpy.float64)
        rhs = numpy.empty([n_rhs], dtype=numpy.float64)

        lib.get_CSR_data.restype = None
        lib.get_CSR_data.argtypes = [ct.c_void_p,
                                     numpy.ctypeslib.ndpointer(
                                         dtype=numpy.int32),
                                     numpy.ctypeslib.ndpointer(
                                         dtype=numpy.int32),
                                     numpy.ctypeslib.ndpointer(
                                         dtype=numpy.float64),
                                     numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        lib.get_CSR_data(self.obj, row, col, data, rhs)

        A = csr_matrix((data, col, row), shape=(n_rhs, n_rhs))
        x = spsolve(A, rhs, use_umfpack=True)

        lib.give_back_sol.restype = None
        lib.give_back_sol.argtypes = [
            ct.c_void_p, numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        lib.give_back_sol(self.obj, x)
