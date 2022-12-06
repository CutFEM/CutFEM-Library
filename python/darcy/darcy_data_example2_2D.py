import numpy as np

c0 = 0
c1 = 100.
c2 = 42.
c3 = 200.
a3 = 50.
c = 10.
mu_G = 2. * c / 1e2

sq_SW = 0 - 1e-10
sq_LGTH = 1. + 2e-10


def func_level_set(P):
    x = P[0]
    y = P[1]
    return (c1 * pow(x, 6) + 20 * x * x - c2 * y * y + c3 * pow(y, 4) - c)


def func_div(P, c, dom):
    x = P[0]
    y = P[1]
    val = (-12 * c3 * y * y + 2 * c2) / np.sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) + c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                                                pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)) - (30 * c1 * pow(x, 4) + 40) / np.sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) + c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) -
                                                                                                                                                                                c2 * y * y - c, 2) + pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)) - ((-4 * c3 * y * y * y + 2 * c2 * y) *
                                                                                                                                                                                                                                                  (2 * (-4 * c3 * y * y * y + 2 * c2 * y) * (-12 * c3 * y * y + 2 * c2) -
                                                                                                                                                                                                                                                   2 * c0 * (-4 * c3 * y * y * y + 2 * c2 * y) *
                                                                                                                                                                                                                                                   (c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c))) / (2 * pow(np.sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                                                                                                                                                                                                                                                                                                                                        c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                                                                                                                                                                                                                                                                                                                                        pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)), 3)) + ((2 * (6 * c1 * pow(x, 5) + 40 * x) * (30 * c1 * pow(x, 4) + 40) +
                                                                                                                                                                                                                                                                                                                                                                                           2 * c0 * (6 * c1 * pow(x, 5) + 40 * x) *
                                                                                                                                                                                                                                                                                                                                                                                           (c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c)) *
                                                                                                                                                                                                                                                                                                                                                                                          (6 * c1 * pow(x, 5) + 40 * x)) / (2 * pow(np.sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                            c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                            pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)), 3))
    return (dom == 0) * val


def func_phat(P, c, dom):
    return (-1. / 8 * mu_G + c / (2. * 1e2))


def fun_dxg(P):
    x = P[0]
    y = P[1]
    return (6 * c1 * pow(x, 5) + 20 * 2 * x)


def fun_dyg(P):
    x = P[0]
    y = P[1]
    return (-2 * c2 * y + 4 * c3 * y * y * y)


def fun_normgradg(P):
    x = P[0]
    y = P[1]
    g = func_level_set(P)
    return np.sqrt(pow(fun_dxg(P), 2) + pow(fun_dyg(P), 2) + c0 * g * g)


def fun_normal(P, c):
    x = P[0]
    y = P[1]
    dxg = fun_dxg(P)
    dyg = fun_dyg(P)
    normgradg = fun_normgradg(P)
    return -((c == 0) * dxg / normgradg + (c == 1) * dyg / normgradg)


def func_neumann(P, c, dom):
    return (dom == 0) * ((c == 0) * (fun_normal(P, 0)) +
                         (c == 1) * (fun_normal(P, 1)))


def func_velocity(P, c, dom):
    return (dom == 0) * ((c == 0) * (fun_normal(P, 0)) +
                         (c == 1) * (fun_normal(P, 1)))


def func_pressure(P, c, dom):
    return (dom == 0) * (func_level_set(P) + c) / 1e2


#
