/*
   We consider a time-dependent bulk problem on Omega2.
   We consider problems with both Neumann and Dirichlet boundary conditions,
   and we consider scheme II, III and the Reynold scheme.

   // Problem:
   Find u in Omega2 such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    on Omega2.

    // Numerical method:s
    A space-time Cutfem, using the level-set method,
    which allows for both dg and cg.

    Classical : Integration by parts on full convection term, no term added to
   make anti-symmetric Conservative: Reynold's transport theorem is used to make
   the bilinear form fulfill a conservation law
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
#include <iostream>
#include <experimental/filesystem>
#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#include "finiteElement.hpp"
#include "../num/gnuplot.hpp"
// #include "GenericMapping.hpp"
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
// #include "projection.hpp"
// #include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"

#include "../algoim/quadrature_general.hpp"

using namespace globalVariable;

namespace Example1 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return -(sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17) - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return -(sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17) - Epsilon;
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// Initial solution surface
R fun_init_surfactant(double *P, const int i) {
    double x = P[0], y = P[1];

    return 10 + 0.5 + 0.4 * cos(pi * x) * cos(pi * y) -
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) -
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

// Exact solution surface
R fun_sol_surfactant(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

    return 10 + 0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) -
           pi / 250 * sin(pi * x) * cos(pi * y) * cos(2 * pi * t) * (x - xc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) -
           pi / 250 * cos(pi * x) * sin(pi * y) * cos(2 * pi * t) * (y - yc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc));
}

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return (pi *
            (31250 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) + 31250 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             31250 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             15625 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             15625 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             31250 * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             31250 * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             3125000 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3125000 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             62500 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             52500 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             1750000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             62500 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4900 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4375 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             31250 * x * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * x * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * (x * x) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * (x * x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * y * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             31250 * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             93750 * (y * y) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * (y * y * y) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4900 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) +
             2450 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             7350 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4375 * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             1372 * pi * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             31250 * x * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * (x * x) * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             187500 * (x * x) * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * (x * x * x) * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             31250 * y * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             187500 * (y * y) * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             62500 * (y * y) * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * (y * y * y) * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             12500000 * (x * x) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             12500000 * (y * y) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             7350 * pi * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) -
             2450 * pi * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) -
             14700 * pi * (cos(t * pi) * cos(t * pi)) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             4900 * pi * (cos(t * pi) * cos(t * pi)) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             4375 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             8750 * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             1372 * pi * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) * sin(y * pi) +
             2744 * pi * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             980000 * (cos(t * pi) * cos(t * pi)) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * pi * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(y * pi) -
             14700 * pi * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(x * pi) +
             4375 * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             2744 * pi * cos(y * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(x * pi) -
             980000 * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             2450 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             686 * (pi * pi) * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             31250 * (x * x) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * (x * x) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * (x * x * x) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * (y * y) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * (y * y) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * (y * y * y) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             2450 * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) +
             7350 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             2450 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             686 * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) -
             1372 * (pi * pi) * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) *
                 sin(y * pi) +
             93750 * (x * x) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             125000 * (x * x * x) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             62500 * (x * x * x * x) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             93750 * (y * y) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             125000 * (y * y * y) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * (y * y * y * y) * (pi * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             2450 * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) +
             7350 * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) -
             2450 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             1372 * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) *
                 sin(x * pi) +
             686 * (pi * pi) * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             3125000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             2450 * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) +
             686 * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) -
             62500 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             62500 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             1562500 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             1562500 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             17500 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             62500 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             31250 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             3125000 * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             8750 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             187500 * x * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             62500 * x * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             125000 * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             187500 * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             125000 * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             62500 * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             6250000 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             12500000 * x * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             12500000 * y * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             52500 * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             35000 * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             17500 * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             17500 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             8750 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1750000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * x * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             26250 * (x * x) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             1372 * x * (pi * pi) * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) *
                 cos(y * pi) -
             17500 * (x * x * x) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4900 * y * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             8750 * (y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * x * (y * y) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * (x * x) * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             4900 * x * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) +
             4900 * x * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * (x * x) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * (x * x) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             4900 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) -
             14700 * y * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * (y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             26250 * (y * y) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 (sin(t * pi) * sin(t * pi) * sin(t * pi)) -
             17500 * (y * y * y) * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             686 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) +
             686 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             14700 * x * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) -
             4900 * x * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * (x * x) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             52500 * (x * x) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             4900 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) +
             14700 * y * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             78750 * (y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * (y * y) * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1372 * y * (pi * pi) * (cos(t * pi) * cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) *
                 sin(y * pi) -
             52500 * (y * y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1372 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) *
                 sin(y * pi) -
             1372 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) *
                 sin(x * pi) -
             14700 * x * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) +
             78750 * (x * x) * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             1372 * x * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi) * sin(t * pi)) * sin(x * pi) *
                 sin(y * pi) -
             52500 * (x * x * x) * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             4900 * y * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) -
             8750 * (y * y) * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             1750000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3125000 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3125000 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             686 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) +
             686 * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             35000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             62500 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             875000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             875000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1750000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             17500 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             35000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             375000 * x * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             125000 * x * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             125000 * x * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             1750000 * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             875000 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             875000 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             105000 * x * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * x * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * x * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             52500 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             17500 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             105000 * y * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * y * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             17500 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             1750000 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             7000000 * y * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             35000 * x * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             105000 * x * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             105000 * y * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             35000 * y * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             35000 * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             7000000 * x * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * (x * x) * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             29400 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             14700 * (y * y) * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             6250000 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             14700 * (x * x) * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) -
             4900 * (y * y) * (pi * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) * sin(y * pi) +
             490000 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             9375000 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3125000 * (x * x) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * (x * x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3125000 * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9375000 * (y * y) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * (y * y * y) * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             17500 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             8750 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             62500 * x * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             62500 * x * y * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * (x * x) * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             490000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             245000 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * x * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * x * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * (x * x) * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             14700 * y * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             52500 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             17500 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             52500 * (y * y) * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             125000 * x * (y * y) * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * (x * x) * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             245000 * pi * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             14700 * x * pi * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) +
             9800 * x * pi * (cos(t * pi) * cos(t * pi)) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             8750 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             52500 * x * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             35000 * (x * x) * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             52500 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * y * pi * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) +
             29400 * y * pi * (cos(t * pi) * cos(t * pi)) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * (pi * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             105000 * (y * y) * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             17500 * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             1372 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) +
             1372 * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             29400 * x * pi * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             105000 * (x * x) * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             9800 * y * pi * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(y * pi) +
             8750 * y * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * (y * y) * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             2744 * pi * cos(t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(2 * t * pi) * sin(y * pi) -
             2744 * pi * (cos(t * pi) * cos(t * pi)) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             19600 * x * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             19600 * y * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             6250000 * x * (y * y) * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * (x * x) * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             17500 * x * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             490000 * x * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * (y * y) * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             35000 * x * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * x * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             490000 * x * pi * cos(2 * t * pi) * cos(x * pi) * (sin(t * pi) * sin(t * pi)) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * (x * x) * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * cos(2 * t * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9800 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             9800 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             17500 * x * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * x * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * x * y * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             9800 * x * y * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             17500 * x * (y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             9800 * x * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (sin(t * pi) * sin(t * pi)) -
             17500 * (x * x) * y * (pi * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * x * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 (sin(t * pi) * sin(t * pi)) +
             9800 * (x * x) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sin(t * pi) -
             9800 * (y * y) * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * (x * x) * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * (y * y) * (pi * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             1372 * x * (pi * pi) * (cos(t * pi) * cos(t * pi)) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) *
                 sin(y * pi) -
             1372 * y * (pi * pi) * cos(t * pi) * cos(2 * t * pi) * (sin(t * pi) * sin(t * pi)) * sin(x * pi) *
                 sin(y * pi) +
             1750000 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1750000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             12500000 * x * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             35000 * x * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             3500000 * x * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1750000 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1750000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             70000 * x * y * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             3500000 * y * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9800 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             9800 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             70000 * x * y * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             980000 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * x * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                       pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))))) /
           (12500 *
            sqrt((pow((y + (7 * cos(pi * t)) / 25 - 1. / 2), 2) + pow(((7 * sin(pi * t)) / 25 - x + 1. / 2), 2))) *
            (-1250 * x - 1250 * y - 350 * cos(t * pi) + 350 * sin(t * pi) + 700 * y * cos(t * pi) -
             700 * x * sin(t * pi) + 1250 * (x * x) + 1250 * (y * y) + 98 * (cos(t * pi) * cos(t * pi)) +
             98 * (sin(t * pi) * sin(t * pi)) + 625));
}
} // namespace Example1

namespace Shi1 {
/* An Eulerian Formulation for Solving Partial Differential Equations
Along a Moving Interface – Jian-Jun Xu, Hong-Kai Zhao. */

R fun_init_surfactant(double *P, const int i) { return P[1] / sqrt(P[0] * P[0] + P[1] * P[1]) + 2.; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-t / 4) * (y / sqrt((x - t) * (x - t) + y * y)) + 2.;
    // return exp(-4*t) * (y / sqrt((x - t) * (x - t) + y * y)) + 2.;
}
R fun_one(double *P, const int cc, const R t) { return 1.; }

R fun_velocity(double *P, const int i) { return (i == 0) ? 1. : 0.; }
R fun_levelSet(double *P, const int i, const R t) { return sqrt((P[0] - t) * (P[0] - t) + P[1] * P[1]) - 2 - Epsilon; }
R fun_levelSet(double *P, const int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 2 - Epsilon; }

// R fun_levelSet(double *P, const int i, const R t) { return (P[0] - t) * (P[0] - t) + P[1] * P[1] - 2 - Epsilon; }
// R fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - 2 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];

    return -1. / 4 * exp(-t / 4) *
           (y * ((x - t) * (x - t) + y * y - 4.) / (pow(((x - t) * (x - t) + y * y), 3. / 2))); // u = exp(-t/4)...

    // return -(y * exp(-4 * t) * (4 * t * t - 8 * t * x + 4 * x * x + 4 * y * y - 1)) /
    //        pow((t * t - 2 * t * x + x * x + y * y), 3. / 2);        // u = exp(-4t)...
}
} // namespace Shi1

namespace Deckelnick {
// "Stability and error analysis for a diffuse interface approach to an advection-diffusion
// equation on a moving surface" – Example 1

R fun_init_surfactant(double *P, const int i) { return P[0] * P[1]; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-4 * t) * (x * y * cos(pi * t) + 0.5 * (x * x - y * y) * sin(pi * t));
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 0.5 * pi * P[1] : -0.5 * pi * P[0]; }
R fun_levelSet(double *P, const int i, const R t) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 1 - Epsilon; }
// R fun_levelSet(double *P, const int i, const R t) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 2 - Epsilon; }
R fun_levelSet(double *P, const int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 1 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return 0.;
}
} // namespace Deckelnick

namespace Deckelnick2 {
// "Stability and error analysis for a diffuse interface approach to an advection-diffusion
// equation on a moving surface" – Example 2

R fun_init_surfactant(double *P, const int i) { return (P[0] + 0.5) * P[1] + 2.; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-4 * t) * (x + 0.5 - 2 * t) * y + 2.; //! ORIGINAL
    // return exp(-t/4) * (x + 0.5 - 2*t)*y + 0.;    //! SLOWER
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 2. : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1]) - 1 - Epsilon;
}
R fun_levelSet(double *P, const int i) { return sqrt((P[0] + 0.5) * (P[0] + 0.5) + P[1] * P[1]) - 1 - Epsilon; }

// R fun_levelSet(double *P, const int i, const R t) { return (P[0]+0.5-2*t)*(P[0]+0.5-2*t) + P[1] * P[1] - 1 - Epsilon;
// } R fun_levelSet(double *P, const int i) { return (P[0]+0.5)*(P[0]+0.5) + P[1] * P[1] - 1 - Epsilon; }

R fun_one(double *P, const int cc, const R t) { return 1.; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    // return 0.;

    return -(2 * y * exp(-4 * t) * (2 * x - 4 * t + 1) *
             (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y - 3)) /
           (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1); //! ORIGINAL

    // return -(y * exp(-t / 4) * (2 * x - 4 * t + 1) *
    //          (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y - 63)) /
    //        (8 * (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1)); //! SLOWER
}
} // namespace Deckelnick2

// Setup two-dimensional class types
const int d = 2;
typedef GFESpace<MeshQuad2> FESpace;
typedef CutFESpace<MeshQuad2> CutSpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<MeshQuad2> Fun_h;
typedef FunFEM<Mesh2> Fun_h_triag;
// typedef Mesh2 Mesh;
// typedef GFESpace<Mesh> FESpace;
// typedef CutFESpace<Mesh> CutSpace;
// typedef TestFunction<d> FunTest;
// typedef FunFEM<Mesh2> Fun_h;

R fun_levelSet_(double *P, const int i) {
    double x = P[0], y = P[1];
    // return sqrt((x - 1) * (x - 1) + (y - 0.65) * (y - 0.65)) - 0.45;
    return x * x + y * y - .5;
}

R fun_one(double *P, const int cc) { return 1.; }

// namespace Deckelnick2 {

//     R fun_init_surfactant(double *P, const int i) { return (P[0]+0.5)*P[1] + 0.;}
//     R fun_sol_surfactant(double *P, const int i, const R t) {
//         R x = P[0], y = P[1];
//         return exp(-4*t) * (x + 0.5 - 2*t)*y + 0.;    
        
//     }
//     R fun_velocity(double *P, const int i) { return (i == 0) ? 2. : 0.; }

//     R fun_levelSet(double *P, const int i, const R t) { return (P[0]+0.5-2*t)*(P[0]+0.5-2*t) + P[1] * P[1] - 1 - Epsilon; }
//     R fun_levelSet(double *P, const int i) { return (P[0]+0.5)*(P[0]+0.5) + P[1] * P[1] - 1 - Epsilon; }

//     R fun_one(double *P, const int cc, const R t) { return 1.; }

//     R fun_rhs(double *P, const int cc, const R t) {
//         R x = P[0], y = P[1];
//         //return 0.;

//         return -(2 * y * exp(-4 * t) * (2 * x - 4 * t + 1) *
//                  (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y - 3)) /
//                (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1);    
//     }

// }   

template <int N> struct Levelsetsaye {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const { 
    
        if constexpr (N == 2) return (x(0)+0.5-2*t) * (x(0)+0.5-2*t) + x(1) * x(1) - 1.0;

        else
            return (x(0)-x(2)) * (x(0)-x(2)) + x(1) * x(1) - 0.5;
    }

    // gradient
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {
        if constexpr (N == 2) return algoim::uvector<T, N>(2.0 * (x(0)+0.5-2*t), 2.0 * x(1));

        else return algoim::uvector<T, N>(2.0 * (x(0)-x(2)), 2.0 * x(1), 2 * (x(2)-x(0)));
    }
};



Levelsetsaye<2> phi;

template <typename M>
double integralSaye(FunFEM<M> &fh, const Interface<M> &interface, int cu) {

    typedef M Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    typedef typename QFB::QuadraturePoint QuadraturePoint;

    double t = -1;
    
    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element();
            iface += interface.next_element()) {
            const int kb = interface.idxElementOfFace(iface); // idx on backMesh
            const R meas = interface.measure(iface);

            const auto &T(interface.get_element(kb));
            const auto &V0(T.at(0)); // vertex 0
            const auto &V2(T.at(2)); // vertex 2   diagonally opposed
            
            algoim::uvector<double, 2> xymin{V0[0], V0[1]};
            algoim::uvector<double, 2> xymax{V2[0], V2[1]};

            algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, 7);
                
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                const R weight = q.nodes.at(ipq).w;

                val += weight * fh.evalOnBackMesh(kb, 0, mip, t, cu, 0, 0);
            }
    }

    double val_receive = 0;
    #ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_SUM);
    #else
    val_receive = val;
    #endif

    return val_receive;
}


template <typename M>
double integralSaye(const ActiveMesh<M> &Th, const FunFEM<M> &fh, int c0) {
   int nb_dom = Th.get_nb_domain();
   
   double val = 0.;
   std::shared_ptr<const ExpressionVirtual> ui =
       std::make_shared<const ExpressionFunFEM<M>>(fh, c0, op_id);
   for (int i = 0; i < nb_dom; ++i) {
      val += integralSaye(Th, ui, i, 0);
   }
   return val;
}

template <typename M>
double integralSaye(const ActiveMesh<M> &Th,
                const std::shared_ptr<const ExpressionVirtual> &fh, int domain,
                int itq) {
    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element();
        k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        if (Th.isInactive(k, itq))
            continue;

        const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
        int kb = Th.idxElementInBackMesh(k);

        const auto &T(Th.Th[kb]);   // get element in background mesh
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed
        
        //for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

        // Coordinates of current quadrilateral element
        algoim::uvector<double, 2> xymin{V0[0], V0[1]};
        algoim::uvector<double, 2> xymax{V2[0], V2[1]};
        // the rectangle is [V0[0], V2[0]] x [V0[1], V2[1]
        algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, 7);
            
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            val += weight * fh->evalOnBackMesh(kb, domain, mip);
        }

        //}
    }
    double val_receive = 0;

#ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_SUM);
#else
    val_receive = val;
#endif
    return val_receive;
}

#define timenot
#define use_h
#define use_tnot

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    const std::string path_output_data = "../output_files/saye_test/data/";
    const std::string path_figures     = "../output_files/saye_test/paraview/";

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_figures);
    }

    algoim::uvector<double, 3> xmin{-1.1, -1.1, .0};
    algoim::uvector<double, 3> xmax{1.1, 1.1, 1.};

    Levelsetsaye<3> phi3D;
    algoim::QuadratureRule<3> q = algoim::quadGen<3>(phi3D, algoim::HyperRectangle<double, 3>(xmin, xmax), 2, 0, 4);
    std::cout << std::setprecision(16);
    std::cout << "Area = " << q.sumWeights() << "\n";
    return 0;
#ifdef saye
    int nx = 10, ny = 10;
    double x0 = -3.0, y0 = -3.0, lx = 6.0, ly = 6.0;
    MeshQuad2 Th(nx, ny, x0, y0, lx, ly);

    algoim::uvector<double, 2> xmin{-1.1, -1.1};
    algoim::uvector<double, 2> xmax{1.1, 1.1};
    
    algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(-1.1, 1.1), 2, -1, 4);
    std::ofstream f(path_figures + "scheme.vtp");
    algoim::outputQuadratureRuleAsVtpXML(q, f);
    std::cout << std::setprecision(16);
    std::cout << "Area = " << q.sumWeights() << "\n";

    FESpace Vh(Th, DataFE<MeshQuad2>::P1);  // continuous basis functions

    FESpace Lh(Th, DataFE<MeshQuad2>::P1);
    Fun_h levelSet(Lh, fun_levelSet_);
    InterfaceLevelSet<MeshQuad2> interface(Th, levelSet);

    ActiveMesh<MeshQuad2> Thi(Th);
    Thi.truncate(interface, 1);
    Paraview<MeshQuad2> writer(Th, path_figures + "Th.vtk");
    writer.add(levelSet, "levelSet.vtk", 0, 1);
    writer.writeActiveMesh(Thi, path_figures + "ActiveMesh" + ".vtk");

    double area = 0.0;
    for (int k = Thi.first_element(); k < Thi.last_element(); k += Thi.next_element()) {

        // const auto &T(Kh_i[k]);     // element k
        const auto &T(Th[Thi.idxElementInBackMesh(k)]);
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed
        
        algoim::uvector<double, 2> xymin{V0[0], V0[1]};
        algoim::uvector<double, 2> xymax{V2[0], V2[1]};

        area += algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, 7).sumWeights();
    }
    std::cout << std::setprecision(16);
    std::cout << area << "\n";

    Fun_h funone(Vh, fun_one);
    double intGamma = integral(funone, interface, 0);
    double intGammaSaye = integralSaye(funone, interface, 0);

    double intInterior = integral(Thi, funone, 0);
    double intInteriorSaye = integralSaye(Thi, funone, 0);

    std::cout << "intGamma = " << intGamma << "\n";
    std::cout << "intGammaSaye = " << intGammaSaye << "\n";
    std::cout << "intInterior = " << intInterior << "\n";
    std::cout << "intInteriorSaye = " << intInteriorSaye << "\n";

#elif defined(time)

    int nx = 10, ny = 10;
    double x0 = -3.0, y0 = -3.0, lx = 6.0, ly = 6.0;
    MeshQuad2 Th(nx, ny, x0, y0, lx, ly);

    algoim::uvector<double, 3> xmin{-1.1, -1.1, 0.};
    algoim::uvector<double, 3> xmax{1.1, 1.1, 1.};
    
    algoim::QuadratureRule<3> q = algoim::quadGen<3>(phi, algoim::HyperRectangle<double, 3>(-1.1, 1.1), 2, -1, 4);
    std::ofstream f(path_figures + "scheme.vtp");
    algoim::outputQuadratureRuleAsVtpXML(q, f);
    std::cout << std::setprecision(16);
    std::cout << "Area = " << q.sumWeights() << "\n";

  
#elif defined(quadratic)
    int nx = 10, ny = 10;
    algoim::uvector<double, 2> xrange_global{-1.1, 1.1};
    algoim::uvector<double, 2> yrange_global{-1.1, 1.1};
    MeshQuad2 Th(nx, ny, xrange_global(0), yrange_global(0), xrange_global(1), yrange_global(1));

    FESpace Lh(Th, DataFE<MeshQuad2>::P1);
    Fun_h levelSet(Lh, fun_levelSet_);
    InterfaceLevelSet<MeshQuad2> interface(Th, levelSet);

    Levelsetsaye<2> phi;
    auto q       = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(-1.1, 1.1), -1, -1, 4);
    double area0 = q([](const auto &x) { return 1.0; });
    std::cout << "area0 = " << area0 << "\n";

    algoim::HyperRectangle<double, 2> rect(xrange_global, yrange_global);

    algoim::QuadratureRule<2> qq =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xrange_global, yrange_global), -1, -1, 4);
    std::ofstream f(path_figures + "scheme.vtp");
    algoim::outputQuadratureRuleAsVtpXML(qq, f);

    int i = 1;
    std::cout << "rect.side " << std::to_string(i) << ": " << rect.side(i) << "\n";

    // const auto &T(Th[0]);     // element k
    // const auto &V0(T.at(0));   // vertex 0
    // const auto &V1(T.at(1));   // vertex 1
    // const auto &V2(T.at(2));   // vertex 2
    // const auto &V3(T.at(3));   // vertex 3
    // std::cout << "Vertex 0: " << V0 << "\n";
    // std::cout << "Vertex 1: " << V1 << "\n";
    // std::cout << "Vertex 2: " << V2 << "\n";
    // std::cout << "Vertex 3: " << V3 << "\n";

    ActiveMesh<MeshQuad2> Kh_i(Th);
    Kh_i.truncate(interface, 1);
    Paraview<MeshQuad2> writer(Th, path_figures + "Th.vtk");
    writer.add(levelSet, "levelSet.vtk", 0, 1);
    writer.writeActiveMesh(Kh_i, path_figures + "ActiveMesh" + ".vtk");

    double area = 0.0;
    for (int k = Kh_i.first_element(); k < Kh_i.last_element(); k += Kh_i.next_element()) {

        // const auto &T(Kh_i[k]);     // element k
        const auto &T(Th[Kh_i.idxElementInBackMesh(k)]);
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed
        algoim::uvector<double, 2> xrange{V0[0], V2[0]};
        algoim::uvector<double, 2> yrange{V0[1], V2[1]};

        algoim::uvector<double, 2> xmin{V0[0], V0[1]};
        algoim::uvector<double, 2> xmax{V2[0], V2[1]};

        area += algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xrange, yrange), -1, -1, 4).sumWeights();
    }
    std::cout << std::setprecision(16);
    std::cout << area << "\n";
    return 0;

#elif defined(triangulation)
    Mesh2 Th_(3, 3, 0, 0, 2, 2);
    GFESpace<Mesh2> Lh_(Th_, DataFE<Mesh2>::P1);
    Fun_h_triag levelSet_(Lh_, fun_levelset_triangulation);
    InterfaceLevelSet<Mesh2> interface_(Th_, levelSet_);

    const auto &T(Th_[3]); // first triangle and the unit triangle
    std::cout << Th_(T) << std::endl;
    std::cout << T.at(0) << std::endl;

    ActiveMesh<Mesh2> Kh_i_(Th_);
    Kh_i_.truncate(interface_, 1);
    Paraview<Mesh2> writer_(Th_, pathOutputFigures + "Th_triangulation.vtk");
    writer_.add(levelSet_, "levelSet.vtk", 0, 1);
    writer_.writeActiveMesh(Kh_i_, pathOutputFigures + "ActiveMesh_triangulation" + ".vtk");

    Th_.info();

#elif defined(cutfem)
    using namespace Deckelnick2;
    // General structure idea of a time-dependent implementation

    // Construct background mesh
    int nx, ny;
    double h  = 0.05; // starting mesh size
    double x0 = -3.0, y0 = -3.0, lx = 6.0, ly = 6.0;
    nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
    MeshQuad2 Th(nx, ny, x0, y0, lx, ly);

    // Set time parameters
    const double t0 = 0.;
    const double tfinal = 1.5; // Final time
    
    const int divisionMeshSize = 3;
    double dT              = h / divisionMeshSize;
    int total_number_iteration = int(tfinal / dT);
    dT        = tfinal / total_number_iteration;
    double time_step = dT;
    double final_time = total_number_iteration * time_step;

    // Problem parameters
    double D = 1.;
    double tau0 = 0, tau1 = 1., tau2 = .1;

    // Set up FE space in time
    Mesh1 Qh(total_number_iteration + 1, t0, final_time);
    // 1D Time space
    FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
    // Quadrature data
    const QuadratureFormular1d &qTime(*Lobatto(3));
    const Uint nbTime       = qTime.n;
    const Uint ndfTime      = Ih[0].NbDoF();
    const Uint lastQuadTime = nbTime - 1;

    // Set up FE space in space

    // ! Skip velocity for now
    // Lagrange2 FEvelocity(0);    //! Need to define this for Mesh2Quad
    // FESpace VelVh(Th, FEvelocity);
    //Fun_h vel(VelVh, fun_velocity);

    FESpace Lh(Th, DataFE<MeshQuad2>::P1);
    double dt_levelSet = dT / (nbTime - 1);
    std::vector<Fun_h> ls(nbTime);
    for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);

    FESpace Vh(Th, DataFE<MeshQuad2>::P1);  // continuous basis functions

    TimeInterface<MeshQuad2> interface(qTime);
    
    InterfaceLevelSet<MeshQuad2> interface(Th, levelSet);

    CutFEM<MeshQuad2> surfactant(qTime);

    ActiveMesh<MeshQuad2> Thi(Th);
    Thi.truncate(interface, 1);     // remove positive part of level set

#endif

    return 0;
}
