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

R fun_one(double *P, const int cc, const R t) { return 1.; }

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

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) -
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) -
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

// Exact solution surface
R fun_sol_surfactant(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) -
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

namespace Shi1Simplified {
/* An Eulerian Formulation for Solving Partial Differential Equations
Along a Moving Interface – Jian-Jun Xu, Hong-Kai Zhao. */

R fun_init_surfactant(double *P, const int i) { return P[1] / 2 + 2.; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-4 * t) * y / 2 + 2.;
    // return exp(-t / 4) * y/2 + 2.;
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 1. : 0.; }
R fun_levelSet(double *P, const int i, const R t) { return sqrt((P[0] - t) * (P[0] - t) + P[1] * P[1]) - 2 - Epsilon; }
R fun_levelSet(double *P, const int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 2 - Epsilon; }
// R fun_levelSet(double *P, const int i, const R t) { return (P[0] - t) * (P[0] - t) + P[1] * P[1] - 2 - Epsilon; }
// R fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - 2 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return (y * exp(-4 * t)) / (2 * ((x - t) * (x - t) + y * y)) - 2 * y * exp(-4 * t); // u = exp(-4t)y/2 + 2
    // return (y * exp(-t / 4)) / (2 * ((x - t) * (x - t) + y * y)) - (y * exp(-t / 4)) / 8;       // u = exp(-t/4)y/2 +
    // 2
}
} // namespace Shi1Simplified

namespace Shi2 {

R fun_init_surfactant(double *P, const int i) { return P[1] / sqrt(P[0] * P[0] + P[1] * P[1]) + 2; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-t / 4) * (y / sqrt((x - 0.1 * t) * (x - 0.1 * t) + y * y)) + 2;
}
R fun_velocity(double *P, const int i, const R t) { return (i == 0) ? 0.1 : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] - 0.1 * t) * (P[0] - 0.1 * t) + P[1] * P[1]) - 2 - Epsilon;
}
R fun_levelSet(double *P, const int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 2 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return -(5 * y * exp(-t / 4) * (t * t - 20 * t * x + 100 * x * x + 100 * y * y - 400)) /
           (2 * pow(t * t - 20 * t * x + 100 * x * x + 100 * y * y, 1.5));
}

} // namespace Shi2

namespace Shi3 {

R fun_init_surfactant(double *P, const int i) { return P[1] / sqrt(P[0] * P[0] + P[1] * P[1]) + 2; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-t / 4) * (y / sqrt((x - 0.9 * t) * (x - 0.9 * t) + y * y)) + 2;
}
R fun_velocity(double *P, const int i, const R t) { return (i == 0) ? 0.9 : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] - 0.9 * t) * (P[0] - 0.9 * t) + P[1] * P[1]) - 2 - Epsilon;
}

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return -(5 * y * exp(-t / 4) * (81 * (t * t) - 180 * t * x + 100 * (x * x) + 100 * (y * y) - 400)) /
           (2 * (81 * (t * t) - 180 * t * x + 100 * (x * x) + 100 * (y * y)) *
            sqrt((81 * (t * t) - 180 * t * x + 100 * (x * x) + 100 * (y * y))));
}

} // namespace Shi3

namespace Shi4 {

R fun_init_surfactant(double *P, const int i) { return P[1] / sqrt(P[0] * P[0] + P[1] * P[1]) + 2; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-t / 4) * (y / sqrt((x - 0.5 * t) * (x - 0.5 * t) + y * y)) + 2;
}
R fun_velocity(double *P, const int i, const R t) { return (i == 0) ? 0.5 : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] - 0.5 * t) * (P[0] - 0.5 * t) + P[1] * P[1]) - 2 - Epsilon;
}

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return -(y * exp(-t / 4) * ((t * t) - 4 * t * x + 4 * (x * x) + 4 * (y * y) - 16)) /
           (2 * ((t * t) - 4 * t * x + 4 * (x * x) + 4 * (y * y)) *
            sqrt(((t * t) - 4 * t * x + 4 * (x * x) + 4 * (y * y))));
}

} // namespace Shi4

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

namespace Deckelnick2ToShi1 {
// Here, we solve using the Deckelnick2 form of solution, but on the
// Shi1 geometry

R fun_init_surfactant(double *P, const int i) { return P[0] * P[1] + 0.; }
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-4 * t) * (x - t) * y + 0.;
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 1. : 0.; }
R fun_levelSet(double *P, const int i, const R t) { return sqrt((P[0] - t) * (P[0] - t) + P[1] * P[1]) - 2 - Epsilon; }
R fun_levelSet(double *P, const int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - 2 - Epsilon; }

R fun_one(double *P, const int cc, const R t) { return 1.; }

// R fun_levelSet(double *P, const int i) { return sqrt((P[0]+0.5)*(P[0]+0.5) + P[1] * P[1]) - 1 - Epsilon; }
// R fun_levelSet(double *P, const int i) { return (P[0]+0.5)*(P[0]+0.5) + P[1] * P[1] - 1 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return (4.0 * y * exp(-4 * t) * (t - x) * (t * t - 2 * t * x + x * x + y * y - 1)) /
           (t * t - 2.0 * t * x + x * x + y * y);
}
} // namespace Deckelnick2ToShi1

namespace Deckelnick2Normalized {
// "Stability and error analysis for a diffuse interface approach to an advection-diffusion
// equation on a moving surface" – Example 2

R fun_init_surfactant(double *P, const int i) {
    R x = P[0], y = P[1];
    return y / sqrt((x + 0.5) * (x + 0.5) + y * y);
}
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-4 * t) * y / sqrt((x + 0.5 - 2 * t) * (x + 0.5 - 2 * t) + y * y);
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 2. : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1]) - 1 - Epsilon;
}
R fun_levelSet(double *P, const int i) { return sqrt((P[0] + 0.5) * (P[0] + 0.5) + P[1] * P[1]) - 1 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return -(32 * y * exp(-4 * t) * (4 * t * t - 4 * t * x - 2 * t + x * x + x + y * y)) /
           ((16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1) *
            sqrt((16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1)));
}
} // namespace Deckelnick2Normalized

namespace Deckelnick2NormalizedSlower {
// "Stability and error analysis for a diffuse interface approach to an advection-diffusion
// equation on a moving surface" – Example 2

R fun_init_surfactant(double *P, const int i) {
    R x = P[0], y = P[1];
    return y / sqrt((x + 0.5) * (x + 0.5) + y * y) + 4.;
}
R fun_sol_surfactant(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return exp(-t / 4) * y / sqrt((x + 0.5 - 2 * t) * (x + 0.5 - 2 * t) + y * y) + 4.;
}

R fun_velocity(double *P, const int i) { return (i == 0) ? 2. : 0.; }
R fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1]) - 1 - Epsilon;
}
R fun_levelSet(double *P, const int i) { return sqrt((P[0] + 0.5) * (P[0] + 0.5) + P[1] * P[1]) - 1 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];
    return -(y * exp(-t / 4) * (16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y - 15)) /
           (2 * ((16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1) *
                 sqrt((16 * t * t - 16 * t * x - 8 * t + 4 * x * x + 4 * x + 4 * y * y + 1))));
}
} // namespace Deckelnick2NormalizedSlower

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// Note: standard for this program is to solve the equations on the inner domain
// Omega 2.

// Choose Discontinuous or Continuous Galerkin method (options: "dg", "cg")
#define cg
// Set numerical example (options: "example1", "shi1", "shi2", "deckelnick", "deckelnick2")
#define deckelnick2
// Set scheme for the dg method (options: "conservative", "classical" see
// thesis. Irrelevant if "cg" is defined instead of "dg")
#define conservative
// Set stabilization method (options: "fullstab", "macro")
#define fullstab
// Decide whether to solve for level set function, or to use exact (options:
// "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h
#define use_tnot

#if defined(frachon1)
using namespace NumericSurfactantEllipse2D;
#elif defined(zahedi1)
using namespace Zahedi1;
#elif defined(zahedi2)
using namespace Zahedi2;
#elif defined(shi1)
using namespace Shi1;
#elif defined(shi2)
using namespace Shi2;
#elif defined(shi3)
using namespace Shi3;
#elif defined(example1)
using namespace Example1;
#elif defined(deckelnick)
using namespace Deckelnick;
#elif defined(deckelnick2)
using namespace Deckelnick2;
#elif defined(deckelnick2toshi1)
using namespace Deckelnick2ToShi1;
#endif

int main(int argc, char **argv) {

    verbose = 0; // 2 for more info

    MPIcf cfMPI(argc, argv);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to
                                 // run only once and plot to paraview)
    int nx = 20, ny = 15;        // starting mesh size (only apply if use_n is defined)
    double h  = 0.1;             // starting mesh size
    // double h  = 0.00625/sqrt(2);             // starting mesh size
    double dT = 0.0625;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#if defined(frachon1)
    // Paths to store data
    const std::string path_output_data = "../output_files/surface/frachon1/data/";
    const std::string path_figures     = "../output_files/surface/frachon1/paraview/";
#elif defined(zahedi1)
    const std::string path_output_data = "../output_files/surface/zahedi1/data/";
    const std::string path_figures     = "../output_files/surface/zahedi1/paraview/";
#elif defined(zahedi2)
    const std::string path_output_data = "../output_files/surface/zahedi2/data/";
    const std::string path_figures     = "../output_files/surface/zahedi2/paraview/";
#elif defined(shi1)
    const std::string path_output_data = "../output_files/surface/shi1/data/";
    const std::string path_figures     = "../output_files/surface/shi1/paraview/";
#elif defined(shi2)
    const std::string path_output_data = "../output_files/surface/shi2/data/";
    const std::string path_figures     = "../output_files/surface/shi2/paraview/";
#elif defined(shi3)
    const std::string path_output_data = "../output_files/surface/shi3/data/";
    const std::string path_figures     = "../output_files/surface/shi3/paraview/";
#elif defined(example1)
    const std::string path_output_data = "../output_files/surface/example1/data/";
    const std::string path_figures     = "../output_files/surface/example1/paraview/";
#elif defined(deckelnick)
    const std::string path_output_data = "../output_files/surface/deckelnick/data/";
    const std::string path_figures     = "../output_files/surface/deckelnick/paraview/";
#elif defined(deckelnick2)
    const std::string path_output_data = "../output_files/surface/deckelnick2/data/";
    const std::string path_figures     = "../output_files/surface/deckelnick2/paraview/";
#elif defined(deckelnick2toshi1)
    const std::string path_output_data = "../output_files/surface/deckelnick2_to_shi1/data/";
    const std::string path_figures     = "../output_files/surface/deckelnick2_to_shi1/paraview/";
#endif

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_figures);
    }

    // Data file to hold problem data
    std::ofstream output_data(path_output_data + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors;             // array to hold bulk errors
    std::array<double, iterations> errors_max;         // array to hold bulk errors
    std::array<int, iterations> numb_stabilized_edges; // array to count stabilized edges
    std::array<double, iterations> reynold_error;
    std::array<double, iterations> gamma_length_h;
    std::array<double, iterations> nxs; // array to hold mesh sizes
    std::array<double, iterations> nys; // array to hold mesh sizes
    std::array<double, iterations> hs;  // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

#if defined(frachon1) or defined(zahedi1)
        const double lx = 4., ly = 4.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -2, -2, lx, ly);
#elif defined(zahedi2) || defined(example1)
        const double lx = 1., ly = 1.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0., 0., lx, ly);
#elif defined(shi1) || defined(shi2) || defined(shi3) || defined(deckelnick2toshi1)
        const double lx = 8., ly = 6.;

#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3, -3, lx, ly);

        // std::string f = "../mesh/square_seb_"+std::to_string(j+1)+".msh";
        // Mesh Th(f.c_str());
#elif defined(deckelnick) || defined(deckelnick2)
        const double lx = 4.8, ly = 4.8;
        // const double lx = 8., ly = 6.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif

        Mesh Th(nx, ny, -2.4, -2.4, lx, ly);
        // Mesh Th(nx, ny, -3, -3, lx, ly);
#endif

        // Parameters
        double tfinal = .5; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        int divisionMeshSize = 8;

        // int divisionMeshSize = 2*3*pi;
        // int divisionMeshSize = 18;

        double dT = h / divisionMeshSize;

        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "Iteration " << j + 1 << "/" << iterations << "\n";
        }
        std::cout << "h = " << h << "\n";
        std::cout << "nx = " << nx << "\n";
        std::cout << "ny = " << ny << "\n";
        std::cout << "dT = " << dT << "\n";

        double D = 1.;

        // CG stabilization parameters
        double tau0 = 0, tau1 = 1., tau2 = .1;

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space

        FESpace2 Vh(Th, DataFE<Mesh>::P1);  // continuous basis functions
        FESpace2 Vh2(Th, DataFE<Mesh>::P2); // higher order space
        FESpace2 Vh3(Th, DataFE<Mesh>::P3); // higher order space

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(3));

        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

// Velocity field
#if defined(shi1) || defined(shi2) || defined(shi3) || defined(deckelnick2)
        Lagrange2 FEvelocity(0);
#elif defined(example1) || defined(deckelnick)
        Lagrange2 FEvelocity(1);
#else
        Lagrange2 FEvelocity(2);
#endif
        FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        FESpace2 Lh2(Th, DataFE<Mesh>::P2);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime), ls2(nbTime);
#if defined(levelsetexact)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);
#elif defined(levelsetsolve)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet);
#endif
        // projection(ls_2[0], ls[nbTime - 1]);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);
        // KN<Mapping2 *> mapping(qTime);

        // Convection-Diffusion Problem Object
        CutFEM<Mesh> surfactant(qTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << "\n";

        int iter = 0;
        // double q0_0, q0_1, qp_0, qp_1;
        double q_init0, q_init1, qp1;
        double intF = 0, intG = 0, intGamma = 0; // hold integrals of rhs and Neumann bcs
        double errL2 = 0., err_max = 0.;

        std::vector<double> gamma_length, error_t;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double tid            = iter * time_step;

            const TimeSlab &In(Ih[iter]);

            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << "\n";
            std::cout << " Time      \t : \t" << current_iteration * time_step << "\n";

            swap(ls[0], ls[lastQuadTime]);
            // computation of the interface in each of the three quadrature
            // points
            for (int i = 0; i < nbTime; ++i) {

#if defined(levelsetexact) && not defined(example2)
                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);
#endif
                interface.init(i, Th, ls[i]);

#if defined(levelsetsolve) || defined(example2)
                // We solve for the level-set using Crank-Nicholson in time
                if (i < lastQuadTime) {
                    LevelSet::move_2D(ls[i], vel, vel, dt_levelSet, ls[i + 1]);
                }
#endif
            }

            // Create active meshes
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            CutSpace Wh(ThGamma, Vh);

            // Data for initial solution
            surfactant.initSpace(Wh, In);

            Rn datau0(surfactant.get_nb_dof(), 0.);
            surfactant.initialSolution(datau0); //! Q: does this only allocate memory? and if so, of what?
            KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(), 0)));
            if (iter == 0)
                interpolate(Wh, datas0, fun_init_surfactant);
            // Rn uh(datau0);
            Fun_h u0(Wh, datau0);

            // Objects needed for the weak form
            Normal n;
            Conormal t;

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega2

//// Assembling linear and bilinear forms

// Partial derivative in time, bulk

// reynold scheme
#ifdef conservative

            surfactant.addBilinear(-innerProduct(u, dt(v)), interface, In);

            surfactant.addBilinear(+innerProduct(u, v), *interface(lastQuadTime), In, (int)lastQuadTime);

            surfactant.addLinear(+innerProduct(u0.expr(), v), *interface(0), In, 0);

            // Added time penalty
            // surfactant.addBilinear(
            //     + innerProduct(u, v)
            //     , *interface(0)
            //     , In
            //     , 0
            // );

            // surfactant.addLinear(
            //     + innerProduct(u0.expr(), v)
            //     , *interface(0)
            //     , In
            //     , 0
            // );

// classical scheme
#elif defined(classical)

            surfactant.addBilinear(+innerProduct(dt(u), v), interface, In);
            surfactant.addBilinear(innerProduct(u, v), *interface(0), In, 0);
            surfactant.addLinear(innerProduct(u0.expr(), v), *interface(0), In, 0);

#endif

            // Scheme for diffusion

            surfactant.addBilinear(+innerProduct(D * gradS(u), gradS(v)), interface, In
                                   //, mapping
            );

#ifdef dg
            // integral on inner edges (E_h)
            surfactant.addBilinear(-innerProduct(D * average(gradS(u) * t, 0.5, -0.5), jump(v)) -
                                       innerProduct(D * jump(u), average(gradS(v) * t, 0.5, -0.5)) +
                                       innerProduct(lambdaA * jump(u), jump(v)),
                                   interface, INTEGRAL_INNER_NODE_2D, In);
#endif

            // Schemes for convection

#if defined(dg) && (defined(conservative))

            // Convection terms

            for (int i = 0; i < nbTime; ++i) {

                // Convection term
                surfactant.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v)))
                                       //- innerProduct(u*dxS(vel[i]) + u*dyS(vel[i]), v)
                                       ,
                                       interface, In, i);

                // Flux terms
                surfactant.addBilinear(+innerProduct(average(vel[i] * t * u, 0.5, -0.5), jump(v)) +
                                           innerProduct(lambdaB * jump(fabs(vel[i] * t) * u), jump(v))
                                       // + innerProduct(average((vel.at(i) * t) * u), jump(v))
                                       // - innerProduct(jump(u), average((vel.at(i) * t) * v))
                                       // + innerProduct(1 * jump(u * fabs(vel.at(i) * n)),
                                       // jump(v))
                                       ,
                                       interface, INTEGRAL_INNER_NODE_2D, In, i);
            }

#elif defined(cg) || defined(classical)

#if defined(classical)

            surfactant.addBilinear(+innerProduct((vel.exprList() * grad(u)), v) + innerProduct(u * divS(vel), v),
                                       interface, In);

#elif defined(conservative)
            surfactant.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), interface, In);
#endif

#endif

            // Stabilization
            double stab_surf_face      = tau1;
            double stab_surf_interface = h * h * tau2;
            double stab_mass           = 0.; // tau1 * h;
            double stab_dt             = 0.; // tau1 *h;

#if defined(fullstab)

            surfactant.addFaceStabilization(+innerProduct(stab_surf_face * jump(grad(u) * n), jump(grad(v) * n)),
                                            ThGamma, In);

            // stabilize in last quadrature point
            double ccend = 1. / In.T.mesure() * 1. / qTime[lastQuadTime].a;
            double ccmid = 1. / In.T.mesure() * 1. / qTime[lastQuadTime - 1].a;
            surfactant.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), ccend * jump(grad(v) * n)),
                                            ThGamma, In, lastQuadTime); //! note: which quadrature point?

            // stabilize time derivative term
            surfactant.addFaceStabilization(-innerProduct(stab_dt * jump(grad(u) * n), jump(grad(dt(v)) * n)), ThGamma,
                                            In);

#elif defined(macro)
            // double gamma = 0.1 + 2.95 * dT;
            TimeMacroElementSurface<Mesh> TimeMacro(ThGamma, interface, qTime, 0.25);
            // gammas.at(j) = gamma;

            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, path_figures + "Th" + std::to_string(iter + 1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);
                // std::cout << writerMacro.get_nb_stab_elems(ThGamma, 0) <<
                // "\n";
                //  domain = 0,

                writerMacro.writeFaceStab(ThGamma, 0,
                                              path_figures + "FullStabilization" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeActiveMesh(ThGamma, path_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0, path_figures + "macro" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroInnerEdge(TimeMacro, 0,
                                                    path_figures + "macro_inner_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroOutterEdge(TimeMacro, 0,
                                                     path_figures + "macro_outer_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeSmallElements(TimeMacro, 0,
                                                   path_figures + "small_element" + std::to_string(iter + 1) + ".vtk");
            }

            surfactant.addFaceStabilization(+innerProduct(tau0 / h / h * jump(u), jump(v)) +
                                                    innerProduct(tau1 * jump(grad(u) * n), jump(grad(v) * n)),
                                                ThGamma, In, TimeMacro);

#endif

            // numb_stabilized_edges.at(j) = surfactant.number_of_stabilized_edges;

            surfactant.addBilinear(+innerProduct(stab_surf_interface * grad(u) * n, grad(v) * n), interface, In);

            //     surfactant.addLinear(
            //         + innerProduct(funrhs.expr(), tau1*h*(dx(v)*vx +
            //         dy(v)*vy)) , interface , In , i
            //     );
            // }

            Fun_h funrhs(Vh, In, fun_rhs);

            // Add RHS on surface
            surfactant.addLinear(+innerProduct(funrhs.expr(), v), interface, In);

            // surfactant.addLinear(fun_rhs, +innerProduct(1., v), interface, In);

            // Compute integrals //! PUT BACK
            // intF = integral(funrhs, In, interface, 0);

            //! PUT BACK
            Fun_h funone(Vh, In, fun_one);
            intGamma = integral(funone, interface(1), 0);
            // std::cout << std::setprecision(16);
            // gamma_length.push_back(intGamma);

            // #ifdef conservative
            // #ifdef fullstab
            // #ifdef use_t
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_conservative_fullstab_h" + std::to_string(h) +
            //                 "_"
            //                                                     +
            //                                                     std::to_string(j
            //                                                     + 1) + ".dat");
            // #else
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_,
            //                                path_output_data +
            //                                "mat_conservative_fullstab_j" +
            //                                std::to_string(j + 1) + ".dat");
            // #endif
            // #elif defined(macro)
            // #ifdef use_t
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_conservative_macro_h" + std::to_string(h) + "_"
            //                                                     +
            //                                                     std::to_string(j
            //                                                     + 1) + ".dat");
            // #else
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_conservative_macro_j" + std::to_string(j + 1) +
            //                 ".dat");
            // #endif
            // #endif
            // #elif defined(classical)
            // #ifdef fullstab
            // #ifdef use_t
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_classical_fullstab_h" + std::to_string(h) + "_"
            //                                                     +
            //                                                     std::to_string(j
            //                                                     + 1) + ".dat");
            // #else
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_classical_fullstab_j" + std::to_string(j + 1) +
            //                 ".dat");
            // #endif
            // #elif defined(macro)
            // #ifdef use_t
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_,
            //                                path_output_data + "mat_classical_macro_h"
            //                                + std::to_string(h) + "_" +
            //                                std::to_string(j + 1) + ".dat");
            // #else
            //             if (iter == total_number_iteration - 1)
            //                 matlab::Export(surfactant.mat_, path_output_data +
            //                 "mat_classical_macro_j" + std::to_string(j + 1) +
            //                 ".dat");
            // #endif
            // #endif
            // #endif

#ifndef USE_MPI
            if ((iter == total_number_iteration - 1) && MPIcf::IamMaster()) {
                matlab::Export(surfactant.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
#elif defined(USE_MPI)
            if (iter == total_number_iteration - 1) {
                matlab::Export(surfactant.mat_[0], path_output_data + "mat_rank_" + std::to_string(MPIcf::my_rank()) +
                                                           "_" + std::to_string(j + 1) + ".dat");
            }
#endif

                // Solve linear system
                surfactant.solve("mumps");

                KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
                datau0 = dw;
                surfactant.saveSolution(datau0);

                // Compute error in Reynold relation
                {

                    //     CutFEM<Mesh> reynold(qTime);
                    //     reynold.initSpace(Wh, In);

                    //     reynold.addBilinear(innerProduct(u, v), ThGamma, (int)lastQuadTime, In);
                    //     reynold.addLinear(innerProduct(u0.expr(), v), ThGamma, 0, In);
                    //     for (int i = 0; i < nbTime; ++i) {
                    //         ExpressionFunFEM<Mesh2> vx(vel[i], 0, op_id);
                    //         ExpressionFunFEM<Mesh2> vy(vel[i], 1, op_id);
                    //         reynold.addBilinear(-innerProduct(dt(u), v) - innerProduct(u, dt(v)) -
                    //                                 innerProduct(dx(u) * vx + dy(u) * vy, v) -
                    //                                 innerProduct(u, dx(v) * vx + dy(v) * vy) -
                    //                                 innerProduct(u * dxS(vel[i]) + u * dyS(vel[i]), v),
                    //                             //- innerProduct(u*divS(vel[i]), v)
                    //                             interface, In, i);
                    //     }

                    //     int N = Wh.NbDoF();
                    //     Rn lhs(2 * N);
                    //     multiply(2 * N, 2 * N, reynold.mat_, uh, lhs);

                    //     lhs -= reynold.rhs_;

                    //     reynold_error.at(j) = lhs.linfty();
                    //     std::cout << " e_r^n = " << reynold_error.at(j) << "\n";
                    //
                }

                {
                    // Compute error
                    Rn sol(Wh.get_nb_dof(), 0.);
                    sol += datau0(SubArray(Wh.get_nb_dof(), 0));
                    sol += datau0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));

                    Fun_h funuh_0(Wh, datau0);
                    Fun_h funuh(Wh, sol);

                    errL2 = L2normSurf(funuh_0, fun_sol_surfactant, *interface(0), tid, 0, 1);
                    std::cout << " t_n -> || u-uex||_2 = " << errL2 << "\n";
                    errL2 = L2normSurf(funuh, fun_sol_surfactant, *interface(lastQuadTime), tid + dT, 0, 1);
                    // err_max = max_norm_surface(funuh, fun_sol_surfactant, *interface(lastQuadTime), tid + dT, 0, 1);
                    std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << "\n";
                    // std::cout << " t_{n+1} -> || u-uex||_max = " << err_max << "\n";

                    // Conservation error

                    //! PUT BACK
                    double q0 = integral(funuh_0, interface(0), 0);
                    double q1 = integral(funuh, interface(lastQuadTime), 0);
                    if (iter == 0) {
                        q_init0 = q0;
                        q_init1 = q1;
                        qp1     = q1;
                        q_init1 = integral(u0, interface(0), 0);
                    }

                    output_data << std::setprecision(10);
                    output_data << tid << "," << (q1 - qp1) << "," << intF << "," << ((q1 - qp1) - intF) << "\n";
                    qp1 = q1;
                }

                error_t.push_back(errL2);
                errors.at(j)         = errL2;
                errors_max.at(j)     = err_max;
                gamma_length_h.at(j) = intGamma;

                if (iterations == 1) {
                    Fun_h sol(Wh, datau0);

                    Paraview<Mesh> writer(ThGamma, path_figures + "surfactant_" + std::to_string(iter + 1) + ".vtk");

                    Fun_h uS_ex(Wh, fun_sol_surfactant, tid);
                    writer.add(u0, "surfactant", 0, 1);
                    writer.add(uS_ex, "surfactant_exact", 0, 1);
                    writer.add(fabs(u0.expr() - uS_ex.expr()), "surfactant_error");
                    writer.add(ls[0], "levelSet", 0, 1);
                    // writer.add(ls[2], "levelSet2", 0, 1);
                }

                if (iterations > 1 && iter == total_number_iteration - 1)
                    output_data << h << "," << dT << "," << errL2 << "\n";

                iter++;
            }

            //! PUT BACK
            // std::cout << "\n";
            // std::cout << "Length of Gamma(t) = [";
            // for (auto & length : gamma_length) {

            //     std::cout << length;

            //     std::cout << ", ";

            // }
            // std::cout << "]"
            // << "\n";

            // std::cout << "\n";
            // std::cout << "Error t = [";
            // for (auto & err : error_t) {

            //     std::cout << err;

            //     std::cout << ", ";

            // }
            // std::cout << "]"
            // << "\n";

            // Refine mesh

#ifdef use_n
            nx *= 2;
            ny *= 2;
#elif defined(use_t)
        dT *= 0.5;
#elif defined(use_h)
        // h *= sqrt(0.5);     //! CHANGE BACK
        h *= 0.5;
#endif
        }

        std::cout << std::setprecision(16);
        std::cout << "\n";
        std::cout << "Errors = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << errors.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";
        // std::cout << "Max Errors = [";
        // for (int i = 0; i < iterations; i++) {

        //     std::cout << errors_max.at(i);
        //     if (i < iterations - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]"
        //           << "\n";

        // std::cout << "\n";
        // std::cout << "Number of stabilized edges = [";
        // for (int i = 0; i < iterations; i++) {

        //     std::cout << numb_stabilized_edges.at(i);
        //     if (i < iterations - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]"
        //           << "\n";

        std::cout << "\n";
        std::cout << "h = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << hs.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";
        std::cout << "nx = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << nxs.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";

        std::cout << "ny = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << nys.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";

        std::cout << "\n";
        std::cout << "dT = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << dts.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";

        std::cout << "\n";
        std::cout << "length of Gamma(t) vs h = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << gamma_length_h.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]"
                  << "\n";

        std::cout << "\n";
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
                  << "[s]"
                  << "\n";

        return 0;
    }
