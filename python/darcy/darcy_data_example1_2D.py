import numpy

sq_SW = 0. + 1e-10
sq_LGTH = 1. + 2e-10


def func_level_set(P):
    shift = 0.5
    interfaceRad = 0.25
    return numpy.sqrt((P[0] - shift) * (P[0] - shift) + (P[1] - shift) * (P[1] - shift)) - interfaceRad


def func_div(P, c, dom):
    rad = 0.25
    r2 = rad*rad
    return -2./r2 if dom == 0 else -4./r2


def func_phat(P, c, dom):
    return 19./12


def func_neumann(P, c, dom):
    rad = 0.25
    shift = 0.5
    r2 = (P[0]-shift)*(P[0]-shift) + (P[1]-shift)*(P[1]-shift)
    radius2 = rad*rad
    return r2/(2*radius2)+3./2.


def func_velocity(P, c, dom):
    rad = 0.25
    shift = 0.5
    r2 = (P[0]-shift)*(P[0]-shift) + (P[1]-shift)*(P[1]-shift)
    rad2 = rad*rad
    if c == 0:
        return -1./rad2*(P[0]-shift) if dom == 0 else -2./rad2*(P[0]-shift)
    if c == 1:
        return -1./rad2*(P[1]-shift) if dom == 0 else -2./rad2*(P[1]-shift)


def func_pressure(P, c, dom):
    rad = 0.25
    shift = 0.5
    r2 = (P[0]-shift)*(P[0]-shift) + (P[1]-shift)*(P[1]-shift)
    rad2 = rad*rad
    if c == 0:
        return r2/(2*rad2)+3./2 if dom == 0 else r2/rad2


#
