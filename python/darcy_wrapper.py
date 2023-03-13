import numpy
import ctypes as ct

from cutfem_wrapper import *


# class DarcyCutFEM(object):
class DarcyCutFEM(ProblemCutFEM):
    def __init__(self):
        lib.Darcy2_new.restype = ct.c_void_p
        lib.Darcy2_new.argtypes = None
        self.obj = lib.Darcy2_new()

    def __del__(self):
        lib.delete_object.argtypes = [ct.c_void_p]
        lib.delete_object.restype = None
        lib.delete_object(self.obj)

    def add_natural_BC(self, f):
        lib.Darcy2_add_natural_BC.restype = None
        lib.Darcy2_add_natural_BC.argtypes = [ct.c_void_p, USER_FUNC]
        lib.Darcy2_add_natural_BC(self.obj, f)
