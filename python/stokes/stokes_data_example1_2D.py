import numpy

sq_SW = 0. - 1e-10
sq_LGTH = 1. + 2e-10


def func_level_set(P):
    shift = 0.5
    interfaceRad = numpy.sqrt(0.25)
    return numpy.sqrt((P[0] - shift) * (P[0] - shift) + (P[1] - shift) * (P[1] - shift)) - interfaceRad


def func_div(P, c, dom):
    return 0


def func_rhs(P, c, dom):
    x = P[0]
    y = P[1]
    if c == 0:
        return 40 * x * x * x - 40 * x * y * y - 32 * y + 16
    if c == 1:
        return -40 * x * x * y + 32 * x + 40 * y * y * y - 16


def func_velocity(P, c, dom):
    x = P[0]
    y = P[1]
    if c == 0:
        return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1)
    if c == 1:
        return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1)


def func_pressure(P, c, dom):
    x = P[0]
    y = P[1]
    return 10 * (x * x - y * y) * (x * x - y * y)


#
