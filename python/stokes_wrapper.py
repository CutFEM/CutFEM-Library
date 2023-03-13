import numpy
import ctypes as ct

from cutfem_wrapper import *


class FictitiousStokesCutFEM(ProblemCutFEM):
    def __init__(self):
        lib.FictitiousStokes_new.restype = ct.c_void_p
        lib.FictitiousStokes_new.argtypes = None
        self.obj = lib.FictitiousStokes_new()

    def __del__(self):
        lib.delete_object.argtypes = [ct.c_void_p]
        lib.delete_object.restype = None
        lib.delete_object(self.obj)


class FictitiousStokesVorticityCutFEM(ProblemCutFEM):

    def __init__(self):
        lib.FictitiousStokesVorticity_new.restype = ct.c_void_p
        lib.FictitiousStokesVorticity_new.argtypes = None
        self.obj = lib.FictitiousStokesVorticity_new()

    def __del__(self):
        lib.delete_object.argtypes = [ct.c_void_p]
        lib.delete_object.restype = None
        lib.delete_object(self.obj)
