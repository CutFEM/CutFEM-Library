/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/

/**
 * @brief Time-dependent convection diffusion equation.
 * @note We consider a time-dependent coupled bulk-surface problem.

 *  Numerical method:
    A space-time Cutfem, using the level-set method.

 *  Classical scheme: Integration by parts on diffusion term.
 *  Conservative scheme: Reynold's transport theorem is used to make
    the bilinear form fulfill a conservation law.
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
#include <iostream>
#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#include "finiteElement.hpp"
// #include "../num/gnuplot.hpp"
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../num/matlab.hpp"
#include "../num/redirectOutput.hpp"
#include "paraview.hpp"

using namespace globalVariable; // to access some globally defined constants

// Numerical examples
namespace Example1 {
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17 - Epsilon;
    // return -sqrt((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.22)*(P[1]-0.22)) - 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(pi * y) * sin(pi * x) * (2 * cos(pi * t) * cos(pi * t) - 1) * ((7 * sin(pi * t)) / 25 - x + 0.5)) /
               (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2))) -
           (pi * cos(pi * x) * sin(pi * y) * (2 * cos(pi * t) * cos(pi * t) - 1) * (y + (7 * cos(pi * t)) / 25 - 0.5)) /
               (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2)));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// // Normal x-direction
// R n1(double *P, const R t) {
//    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
//    return (P[0] - xc) /
//           (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// // Normal y-direction
// R n2(double *P, const R t) {
//    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
//    return (P[1] - yc) /
//           (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// Initial solution surface
R fun_uSurfInit(double *P, const int i) {
    double x = P[0], y = P[1];
    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) +
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) +
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double r0 = 1., x = P[0], y = P[1];
    return (2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 -
           (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (y + (7 * cos(pi * t)) / 25 - 1. / 2)) /
               (250 *
                sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow(((7 * sin(t * pi)) / 25 - x + 1. / 2), 2)))) +
           (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * ((7 * sin(pi * t)) / 25 - x + 1. / 2)) /
               (250 *
                sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow(((7 * sin(t * pi)) / 25 - x + 1. / 2), 2)))) +
           1. / 2;
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * pi * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 125 -
           (4 * pi * cos(pi * x) * cos(pi * y) * sin(2 * pi * t)) / 5 -
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (x - 1. / 2)) / 5 +
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (y - 1. / 2)) / 5;

    // return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
    //        (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
    //        (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 0.5)) / 5 +
    //        (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 0.5)) / 5;
}

// RHS fB bulk
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi *
            (46875 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) + 46875 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             2450 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             7350 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             1372 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
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
             31250 * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             125000 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             125000 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             35000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             8750 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             8750 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             35000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             31250 * x * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * y * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             7350 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             2450 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             1372 * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             52500 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             52500 * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             1750000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             17500 * x * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             52500 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             17500 * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             62500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             4900 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             4900 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             4900 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             31250 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             31250 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             93750 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4900 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             2450 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             7350 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4375 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             1372 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             31250 * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * x * x * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             187500 * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * x * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             31250 * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             187500 * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             62500 * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * y * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             12500000 * x * x * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             12500000 * y * y * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             7350 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
             2450 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             14700 * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             4900 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             8750 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             1372 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             2744 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             62500 * x * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             980000 * cos(t * pi) * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             14700 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             4375 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             2744 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             4900 * x * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             17500 * x * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             14700 * y * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * y * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             980000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             14700 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             52500 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             4900 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1372 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             1372 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             31250 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             7350 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             686 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) -
             1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             125000 * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             62500 * x * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             93750 * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             125000 * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * y * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             7350 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             1372 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             3125000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             2450 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             686 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
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
             8750 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             62500 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * x * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             1750000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             26250 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) *
                 cos(y * pi) -
             17500 * x * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4900 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             8750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * x * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             4900 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
             14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             26250 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                 sin(t * pi) -
             17500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             14700 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             52500 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             78750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) *
                 sin(y * pi) -
             52500 * y * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1372 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
             1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             14700 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             78750 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             1372 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sin(y * pi) -
             52500 * x * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             4900 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
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
             686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
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
             35000 * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             7000000 * x * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * x * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             29400 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             14700 * y * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             6250000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             35000 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             14700 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             4900 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             490000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             9375000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3125000 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3125000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9375000 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             9800 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             9800 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             8750 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             62500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             62500 * x * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             490000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * x * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             14700 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             52500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             52500 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             125000 * x * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * x * x * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             245000 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             14700 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             9800 * x * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             8750 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             52500 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             35000 * x * x * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             52500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             29400 * y * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             105000 * y * y * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             17500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             1372 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             1372 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             29400 * x * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             105000 * x * x * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             9800 * y * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             8750 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * y * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             2744 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             2744 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             19600 * x * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             19600 * y * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             6250000 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             490000 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             35000 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             490000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * x * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             9800 * x * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             17500 * x * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             9800 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
             17500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                 sin(t * pi) +
             9800 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sin(t * pi) -
             9800 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * x * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) *
                 sin(y * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
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
             700 * x * sin(t * pi) + 1250 * x * x + 1250 * y * y + 98 * cos(t * pi) * cos(t * pi) +
             98 * sin(t * pi) * sin(t * pi) + 625));
}

R fun_neumann_left(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_bottom(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}

R fun_neumann_right(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_top(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}

} // namespace Example1

namespace Example1_new {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17 - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(pi * y) * sin(pi * x) * (2 * cos(pi * t) * cos(pi * t) - 1) * ((7 * sin(pi * t)) / 25 - x + 0.5)) /
               (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2))) -
           (pi * cos(pi * x) * sin(pi * y) * (2 * cos(pi * t) * cos(pi * t) - 1) * (y + (7 * cos(pi * t)) / 25 - 0.5)) /
               (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2)));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// Initial solution surface
R fun_uSurfInit(double *P, const int i) {
    double x = P[0], y = P[1];

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) +
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) +
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

R fun_uSurfInitT(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) +
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) +
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) +
           pi / 250 * sin(pi * x) * cos(pi * y) * cos(2 * pi * t) * (x - xc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) +
           pi / 250 * cos(pi * x) * sin(pi * y) * cos(2 * pi * t) * (y - yc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc));
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * pi * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 125 -
           (4 * pi * cos(pi * x) * cos(pi * y) * sin(2 * pi * t)) / 5 -
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (x - 1. / 2)) / 5 +
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (y - 1. / 2)) / 5;
}

// RHS fS surface
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return -(pi *
             (46875 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              46875 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              2450 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              7350 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              1372 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              31250 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              15625 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              15625 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              93750 * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              31250 * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              31250 * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
              3125000 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3125000 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              31250 * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              125000 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              125000 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              31250 * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              35000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              8750 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              8750 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              35000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              31250 * x * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              93750 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * x * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              93750 * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              31250 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * y * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              7350 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              2450 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              1372 * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              52500 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              17500 * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
              35000 * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
              17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              52500 * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              17500 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              1750000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              17500 * x * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              52500 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              17500 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              17500 * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              62500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              62500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              4900 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              4900 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              4900 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              31250 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              31250 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              93750 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              31250 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              31250 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              31250 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              31250 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              93750 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              62500 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              4900 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
              2450 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              7350 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              4375 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              1372 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              31250 * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              62500 * x * x * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              187500 * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              125000 * x * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              31250 * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              187500 * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              62500 * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              125000 * y * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              12500000 * x * x * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              12500000 * y * y * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              7350 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
              2450 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
              14700 * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              4900 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              8750 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              8750 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              1372 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              2744 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              62500 * x * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * x * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              980000 * cos(t * pi) * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              4900 * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
              14700 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
              4375 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              2744 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              4900 * x * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              17500 * x * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              14700 * y * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              52500 * y * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              980000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              14700 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              52500 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              4900 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              17500 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              1372 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              1372 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              31250 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              93750 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * x * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              93750 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              31250 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
              7350 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              686 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) -
              1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              93750 * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              125000 * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              62500 * x * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              93750 * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              125000 * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              62500 * y * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              7350 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
              2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              1372 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              3125000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              2450 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              686 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              62500 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              62500 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              1562500 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              1562500 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              17500 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              62500 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              31250 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              31250 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              62500 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              3125000 * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              8750 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              17500 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              17500 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              187500 * x * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              62500 * x * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              125000 * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
              187500 * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              125000 * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              62500 * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              6250000 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              12500000 * x * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              12500000 * y * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              52500 * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              35000 * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              17500 * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              17500 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              8750 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              62500 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              62500 * x * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              1750000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3500000 * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              26250 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) *
                  cos(y * pi) -
              17500 * x * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              4900 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              8750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              62500 * x * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              62500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              4900 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
              4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              17500 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              8750 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
              14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              52500 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              26250 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              1372 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                  sin(t * pi) -
              17500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
              686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              14700 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
              4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              8750 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              52500 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              78750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              17500 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) *
                  sin(y * pi) -
              52500 * y * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              1372 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
              1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              14700 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              78750 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              1372 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                  sin(y * pi) -
              52500 * x * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              4900 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              8750 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              1750000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3125000 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3125000 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              6250000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              35000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              62500 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              62500 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              875000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              875000 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              1750000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * x * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              17500 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              17500 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              35000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              17500 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              52500 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              375000 * x * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              125000 * x * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              125000 * x * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              1750000 * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              875000 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              875000 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              105000 * x * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              35000 * x * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              35000 * x * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
              52500 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              17500 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              105000 * y * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              35000 * y * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
              17500 * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              17500 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              1750000 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              7000000 * y * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              4900 * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              4900 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              35000 * x * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
              105000 * x * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              105000 * y * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              35000 * y * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
              35000 * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              35000 * x * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              7000000 * x * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              4900 * x * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              29400 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
              9800 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
              9800 * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              14700 * y * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              6250000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              6250000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              35000 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              14700 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              4900 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              490000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              9375000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3125000 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3125000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              9375000 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              6250000 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              9800 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              9800 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              8750 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              62500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              62500 * x * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              62500 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
              62500 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              490000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              4900 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              8750 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              17500 * x * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              14700 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              52500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              17500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
              17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              52500 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              125000 * x * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              125000 * x * x * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              245000 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              245000 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              14700 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
              9800 * x * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
              8750 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              17500 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              52500 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              35000 * x * x * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              52500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              4900 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              29400 * y * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
              35000 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              17500 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              105000 * y * y * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
              17500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              1372 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
              1372 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
              4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              29400 * x * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              35000 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              105000 * x * x * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
              9800 * y * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
              8750 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              35000 * y * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
              2744 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
              2744 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
              19600 * x * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
              19600 * y * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
              6250000 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              6250000 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              490000 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              490000 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3500000 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              35000 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
              17500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              490000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              3500000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              490000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
              35000 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
              17500 * x * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              9800 * x * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
              17500 * x * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
              9800 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
              17500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              1372 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                  sin(t * pi) +
              9800 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
              1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sin(t * pi) -
              9800 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
              17500 * x * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
              17500 * x * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
              3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              6250000 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              6250000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) *
                  sin(y * pi) -
              1372 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                  sin(y * pi) -
              1750000 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              1750000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              12500000 * x * y * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              35000 * x * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
              3500000 * x * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              1750000 * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              1750000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              70000 * x * y * pi * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) -
              35000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              3500000 * y * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
              9800 * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
              9800 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
              70000 * x * y * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
              980000 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * x * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
              3500000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                  sqrt((pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) +
                        pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))))) /
           (12500 *
            sqrt((pow((y + (7 * cos(pi * t)) / 25 - 1. / 2), 2) + pow(((7 * sin(pi * t)) / 25 - x + 1. / 2), 2))) *
            (-1250 * x - 1250 * y - 350 * cos(t * pi) + 350 * sin(t * pi) + 700 * y * cos(t * pi) -
             700 * x * sin(t * pi) + 1250 * x * x + 1250 * y * y + 98 * cos(t * pi) * cos(t * pi) +
             98 * sin(t * pi) * sin(t * pi) + 625));
}

R fun_neumann_left(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_bottom(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}

R fun_neumann_right(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_top(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}
} // namespace Example1_new

namespace Example1_omega2 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return -(sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17) - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return -(sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17) - Epsilon;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // In Omega_2(t), there is no Neumann boundary condition, since the only boundary is Gamma(t)
    // to which only the coupling conditions apply
    return 0;
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// Initial solution surface
R fun_uSurfInit(double *P, const int i) {
    double x = P[0], y = P[1];

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) -
           pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) -
           pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
}

// R fun_uSurfInitT(double *P, const int i, const R t) {
//     double x = P[0], y = P[1];

//     return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) -
//         pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) -
//         pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5));
// }

// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

    return 0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) -
           pi / 250 * sin(pi * x) * cos(pi * y) * cos(2 * pi * t) * (x - xc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) -
           pi / 250 * cos(pi * x) * sin(pi * y) * cos(2 * pi * t) * (y - yc) /
               sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc));
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * pi * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 125 -
           (4 * pi * cos(pi * x) * cos(pi * y) * sin(2 * pi * t)) / 5 -
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (x - 1. / 2)) / 5 +
           (2 * pi * pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (y - 1. / 2)) / 5;
}

// RHS fS surface
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // new version
    // return (pi*(15625*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 15625*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 2450*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) -
    // 7350*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) -
    // 1372*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(x*pi) + 31250*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 15625*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) - 15625*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 93750*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 31250*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 31250*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) + 3125000*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 3125000*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25
    // + 1./2),2))) - 6250000*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 31250*x*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 31250*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)
    // + 8750*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 8750*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 31250*(x*x)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) - 93750*(x*x)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 62500*(x*x*x)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 93750*(y*y)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 31250*(y*y)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 62500*(y*y*y)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 7350*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 2450*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 1372*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 52500*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) - 17500*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) -
    // 35000*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) - 17500*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*x*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 52500*y*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 1750000*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 3500000*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) + 17500*x*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 52500*x*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) + 17500*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 17500*y*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) + 62500*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 62500*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) + 4900*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 4900*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 4900*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 4375*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) + 31250*x*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 31250*x*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 93750*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 31250*(x*x)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 62500*(x*x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 31250*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 31250*y*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 31250*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) - 93750*(y*y)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 62500*(y*y*y)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 4900*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))
    // + 2450*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 7350*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 8750*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 8750*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 4375*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 1372*pi*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 31250*x*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 62500*(x*x)*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 187500*(x*x)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) + 125000*(x*x*x)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) +
    // 31250*y*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 187500*(y*y)*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 62500*(y*y)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) + 125000*(y*y*y)*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 12500000*(x*x)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 12500000*(y*y)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 7350*pi*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) -
    // 2450*pi*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) -
    // 14700*pi*(cos(t*pi)*cos(t*pi))*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 4900*pi*(cos(t*pi)*cos(t*pi))*cos(y*pi)*sin(2*t*pi)*sin(x*pi) +
    // 4375*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 8750*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 8750*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 1372*pi*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(y*pi) +
    // 2744*pi*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(x*pi)*sin(2*t*pi)*sin(y*pi) +
    // 62500*x*(y*y)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 62500*(x*x)*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 980000*(cos(t*pi)*cos(t*pi))*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((-
    // x + (7*sin(t*pi))/25 + 1./2),2))) - 4900*pi*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(y*pi) -
    // 14700*pi*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(x*pi) +
    // 4375*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 2744*pi*cos(y*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(x*pi) +
    // 4900*x*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 17500*(x*x)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 14700*y*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 52500*(y*y)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 980000*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((-
    // x + (7*sin(t*pi))/25 + 1./2),2))) + 14700*x*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) -
    // 52500*(x*x)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 4900*y*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) -
    // 17500*(y*y)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 1372*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) -
    // 1372*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 2450*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 686*(pi*pi)*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 31250*(x*x)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 93750*(x*x)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 62500*(x*x*x)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 93750*(y*y)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 31250*(y*y)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 62500*(y*y*y)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 2450*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) +
    // 7350*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 2450*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 686*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi)) -
    // 1372*(pi*pi)*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 93750*(x*x)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 125000*(x*x*x)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 62500*(x*x*x*x)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 93750*(y*y)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 125000*(y*y*y)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 62500*(y*y*y*y)*(pi*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)
    // - 2450*(pi*pi)*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) +
    // 7350*(pi*pi)*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) -
    // 2450*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 1372*(pi*pi)*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(x*pi) +
    // 686*(pi*pi)*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 3125000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 2450*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) +
    // 686*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) -
    // 62500*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) - 62500*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 1562500*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 1562500*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) - 17500*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 62500*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 31250*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 31250*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) + 62500*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 3125000*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 8750*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 17500*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 17500*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 187500*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) + 62500*x*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) +
    // 125000*x*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) - 187500*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 125000*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) + 62500*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) -
    // 6250000*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 12500000*x*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) - 6250000*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 12500000*y*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 52500*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 35000*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) + 17500*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) +
    // 17500*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) - 8750*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 62500*x*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) - 62500*x*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 1750000*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 3500000*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 4900*x*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 26250*(x*x)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 1372*x*(pi*pi)*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 17500*(x*x*x)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 4900*y*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 8750*(y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 62500*x*(y*y)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 62500*(x*x)*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 4900*x*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) +
    // 4900*x*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 17500*(x*x)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 8750*(x*x)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 4900*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) -
    // 14700*y*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 52500*(y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 26250*(y*y)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 1372*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi)) -
    // 17500*(y*y*y)*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 686*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) +
    // 686*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 14700*x*(pi*pi)*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) -
    // 4900*x*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 8750*(x*x)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 52500*(x*x)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 4900*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) +
    // 14700*y*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 78750*(y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 17500*(y*y)*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 1372*y*(pi*pi)*(cos(t*pi)*cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 52500*(y*y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 6250000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 6250000*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 1372*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) -
    // 1372*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 14700*x*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) +
    // 78750*(x*x)*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 1372*x*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) -
    // 52500*(x*x*x)*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) +
    // 4900*y*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) -
    // 8750*(y*y)*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 1750000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 6250000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) + 3125000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y
    // + (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 3125000*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 6250000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 686*(pi*pi)*cos(t*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) +
    // 686*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) +
    // 35000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) + 62500*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 62500*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 875000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 875000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 1750000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 6250000*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 6250000*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) + 17500*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 35000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 17500*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) - 52500*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 375000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 125000*x*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 125000*x*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) - 1750000*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y
    // + (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 875000*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 875000*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 105000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) - 35000*x*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 35000*x*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) - 52500*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 17500*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) - 105000*y*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 35000*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) - 17500*y*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 17500*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 1750000*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 7000000*y*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 4900*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 4900*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 35000*x*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) + 105000*x*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) -
    // 105000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) + 35000*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) +
    // 35000*y*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) + 35000*x*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 7000000*x*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 4900*(x*x)*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 29400*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) +
    // 9800*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) +
    // 9800*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) -
    // 14700*(y*y)*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 6250000*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 6250000*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) - 35000*x*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 14700*(x*x)*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) -
    // 4900*(y*y)*(pi*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) +
    // 490000*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) + 9375000*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 3125000*(x*x)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 6250000*(x*x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 3125000*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 9375000*(y*y)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 6250000*(y*y*y)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 9800*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 9800*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 17500*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 8750*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) - 62500*x*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 62500*x*y*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) - 62500*x*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 62500*(x*x)*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 490000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 245000*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 245000*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) - 4900*x*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 8750*x*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 17500*(x*x)*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 14700*y*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 52500*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) -
    // 17500*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) -
    // 17500*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 52500*(y*y)*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) + 125000*x*(y*y)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) +
    // 125000*(x*x)*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) +
    // 245000*pi*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 245000*pi*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) - 14700*x*pi*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi) +
    // 9800*x*pi*(cos(t*pi)*cos(t*pi))*cos(y*pi)*sin(2*t*pi)*sin(x*pi) +
    // 8750*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 17500*x*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 52500*x*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 35000*(x*x)*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) +
    // 52500*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 4900*y*pi*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) +
    // 29400*y*pi*(cos(t*pi)*cos(t*pi))*cos(x*pi)*sin(2*t*pi)*sin(y*pi) -
    // 35000*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 17500*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 17500*y*(pi*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 105000*(y*y)*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) +
    // 17500*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 1372*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi) +
    // 1372*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 4900*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) +
    // 4900*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 29400*x*pi*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(x*pi) -
    // 35000*x*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 105000*(x*x)*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) +
    // 9800*y*pi*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(y*pi) +
    // 8750*y*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 35000*(y*y)*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) +
    // 2744*pi*cos(t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(2*t*pi)*sin(y*pi) -
    // 2744*pi*(cos(t*pi)*cos(t*pi))*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) -
    // 19600*x*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) -
    // 19600*y*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) -
    // 6250000*x*(y*y)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 6250000*(x*x)*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 17500*x*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 490000*x*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 490000*y*pi*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 3500000*(y*y)*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x
    // + (7*sin(t*pi))/25 + 1./2),2))) + 35000*x*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) +
    // 17500*x*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 490000*x*pi*cos(2*t*pi)*cos(x*pi)*(sin(t*pi)*sin(t*pi))*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) +
    // pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 3500000*(x*x)*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x
    // + (7*sin(t*pi))/25 + 1./2),2))) + 490000*y*pi*cos(2*t*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sqrt((pow((y
    // + (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 9800*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 9800*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 17500*x*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) -
    // 35000*x*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 9800*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 9800*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 17500*x*y*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 9800*x*y*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi) -
    // 17500*x*(y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) +
    // 9800*x*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) -
    // 17500*(x*x)*y*(pi*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 1372*x*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*(sin(t*pi)*sin(t*pi)) +
    // 9800*(x*x)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 1372*y*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) -
    // 9800*(y*y)*(pi*pi)*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) +
    // 17500*(x*x)*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) +
    // 17500*x*(y*y)*(pi*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) +
    // 3500000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 6250000*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25
    // - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 6250000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 1372*x*(pi*pi)*(cos(t*pi)*cos(t*pi))*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) -
    // 1372*y*(pi*pi)*cos(t*pi)*cos(2*t*pi)*(sin(t*pi)*sin(t*pi))*sin(x*pi)*sin(y*pi) +
    // 1750000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 3500000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 1750000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 3500000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 12500000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 35000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) +
    // 3500000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 3500000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 1750000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) + 1750000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) +
    // 70000*x*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) - 35000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) -
    // 3500000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2))) - 9800*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) +
    // 9800*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) -
    // 70000*x*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) -
    // 980000*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((-
    // x + (7*sin(t*pi))/25 + 1./2),2))) - 3500000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y +
    // (7*cos(t*pi))/25 - 1./2),2) + pow((- x + (7*sin(t*pi))/25 + 1./2),2))) -
    // 3500000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y + (7*cos(t*pi))/25 - 1./2),2) + pow((- x +
    // (7*sin(t*pi))/25 + 1./2),2)))))/(12500*sqrt((pow((y + (7*cos(pi*t))/25 - 1./2),2) + pow(((7*sin(pi*t))/25 - x
    // + 1./2),2)))*(- 1250*x - 1250*y - 350*cos(t*pi) + 350*sin(t*pi) + 700*y*cos(t*pi) - 700*x*sin(t*pi) + 1250*(x*x)
    // + 1250*(y*y) + 98*(cos(t*pi)*cos(t*pi)) + 98*(sin(t*pi)*sin(t*pi)) + 625));

    return (pi *
            (46875 * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) + 46875 * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             2450 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             7350 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             1372 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
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
             31250 * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             125000 * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             125000 * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             35000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             8750 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             8750 * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             35000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             31250 * x * x * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * y * y * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             7350 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             2450 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             1372 * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             52500 * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * pi * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * x * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             52500 * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             17500 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             1750000 * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3500000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             17500 * x * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             52500 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             17500 * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             62500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             4900 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             4900 * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             4900 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             31250 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             31250 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             31250 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             93750 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4900 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             2450 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             7350 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             4375 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             1372 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             31250 * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * x * x * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             187500 * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * x * x * x * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             31250 * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             187500 * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             62500 * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * y * y * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             12500000 * x * x * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             12500000 * y * y * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             7350 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
             2450 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             14700 * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             4900 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             4375 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             8750 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             1372 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             2744 * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             62500 * x * y * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             980000 * cos(t * pi) * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             14700 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             4375 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             2744 * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             4900 * x * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             17500 * x * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             14700 * y * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * y * y * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             980000 * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             14700 * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             52500 * x * x * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             4900 * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * y * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1372 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             1372 * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             31250 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             93750 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             31250 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             7350 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             686 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) -
             1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             93750 * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             125000 * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             62500 * x * x * x * x * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             93750 * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             125000 * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             62500 * y * y * y * y * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             2450 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             7350 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             2450 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             1372 * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             3125000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             2450 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             686 * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
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
             8750 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             62500 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * x * y * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             1750000 * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * cos(t * pi) * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             26250 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) *
                 cos(y * pi) -
             17500 * x * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             4900 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             8750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             62500 * x * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             62500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             4900 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             8750 * x * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
             14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             52500 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             26250 * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                 sin(t * pi) -
             17500 * y * y * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             14700 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) -
             4900 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             52500 * x * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             4900 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             14700 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             78750 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             17500 * y * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) *
                 sin(y * pi) -
             52500 * y * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             6250000 * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             1372 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) -
             1372 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             14700 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             78750 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             1372 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sin(y * pi) -
             52500 * x * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             4900 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             8750 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
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
             686 * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             686 * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
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
             35000 * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             7000000 * x * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             4900 * x * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             29400 * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             9800 * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             14700 * y * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             6250000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             35000 * x * y * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             14700 * x * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             4900 * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             490000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             9375000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             3125000 * x * x * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3125000 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9375000 * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * y * y * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             9800 * x * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             9800 * y * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             8750 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             62500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             62500 * x * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             62500 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) +
             62500 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             490000 * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             4900 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             8750 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * x * x * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             14700 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             52500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) -
             17500 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) -
             17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             52500 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             125000 * x * y * y * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             125000 * x * x * y * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             245000 * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             245000 * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             14700 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) +
             9800 * x * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) +
             8750 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             52500 * x * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             35000 * x * x * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             52500 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             29400 * y * pi * cos(t * pi) * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) -
             35000 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             17500 * y * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             105000 * y * y * pi * cos(t * pi) * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) +
             17500 * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             1372 * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) +
             1372 * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) +
             4900 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) +
             29400 * x * pi * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             35000 * x * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             105000 * x * x * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             9800 * y * pi * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) +
             8750 * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * y * y * pi * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) +
             2744 * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             2744 * pi * cos(t * pi) * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             19600 * x * pi * cos(t * pi) * cos(x * pi) * sin(t * pi) * sin(2 * t * pi) * sin(y * pi) -
             19600 * y * pi * cos(t * pi) * cos(y * pi) * sin(t * pi) * sin(2 * t * pi) * sin(x * pi) -
             6250000 * x * y * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * x * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             490000 * x * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * y * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             35000 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) +
             17500 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             490000 * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             3500000 * x * x * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             490000 * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             17500 * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) -
             35000 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi) -
             9800 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(t * pi) * sin(x * pi) -
             17500 * x * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) -
             9800 * x * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) -
             17500 * x * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) +
             9800 * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) * sin(t * pi) -
             17500 * x * x * y * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) *
                 sin(t * pi) +
             9800 * x * x * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sin(t * pi) -
             9800 * y * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(t * pi) +
             17500 * x * x * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) +
             17500 * x * y * y * pi * pi * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) * sin(y * pi) +
             3500000 * y * pi * cos(t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) +
             6250000 * x * y * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             6250000 * x * y * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 sqrt(
                     (pow((y + (7 * cos(t * pi)) / 25 - 1. / 2), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 1. / 2), 2))) -
             1372 * x * pi * pi * cos(t * pi) * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(x * pi) *
                 sin(y * pi) -
             1372 * y * pi * pi * cos(t * pi) * cos(2 * t * pi) * sin(t * pi) * sin(t * pi) * sin(x * pi) *
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
             700 * x * sin(t * pi) + 1250 * x * x + 1250 * y * y + 98 * cos(t * pi) * cos(t * pi) +
             98 * sin(t * pi) * sin(t * pi) + 625));
}

R fun_neumann_left(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_bottom(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}

R fun_neumann_right(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
}

R fun_neumann_top(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
}
} // namespace Example1_omega2

namespace Example2 {
R fun_levelSet(double *P) { return +(sqrt((P[0] - 0.1) * (P[0] - 0.1) + (P[1] - 0.0) * (P[1] - 0.0)) - 0.3) + Epsilon; }

// R fun_levelSet(double *P, const int i, const R t) {
//     return +(sqrt((P[0]-0.1)*(P[0]-0.1) + (P[1]-0.0)*(P[1]-0.0)) - 0.3) + Epsilon;
// }

R fun_rhsBulk(double *P, const int i, const R t) { return 0.; }
R fun_rhsSurf(double *P, const int i, const R t) { return 0.; }
R fun_neumann_Gamma(double *P, const int i) { return 0.; }
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return -0.5 * (1 + cos(M_PI * P[0])) * sin(M_PI * P[1]);
    else
        return 0.5 * (1 + cos(M_PI * P[1])) * sin(M_PI * P[0]);
}
//    R fun_velocity(double *P, const int i){
//        return 0.0;
//    }
R fun_uSurfInit(double *P, int i) { return 0.; }
R fun_w(double r) { return 0.5 * (1 - cos((r - 0.3) * M_PI / (0.5 * 0.3))); }
R fun_uBulkInit(double *P, int i) {
    double r = sqrt((P[0] - 0.1) * (P[0] - 0.1) + P[1] * P[1]);
    if (r > 1.5 * 0.3)
        return 0.5 * (1 - P[0] * P[0]) * (1 - P[0] * P[0]);
    else if (0.3 < r && r <= 1.5 * 0.3)
        return 0.5 * (1 - P[0] * P[0]) * (1 - P[0] * P[0]) * fun_w(r);
    else
        return 0.;
}

// below are just so that program does not crash, it is not known
R fun_uBulkD(double *P, const int i, const int d, const R t) { return 0.; }

R fun_uBulk(double *P, const int i, const R t) { return 0.; }

R fun_uSurf(double *P, const int i, const R t) { return 0.; }
} // namespace Example2

namespace Kex {

    // Level-set function
    double fun_levelSet(double *P, const int i, const R t) {
        R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
        return sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17 - Epsilon;
    }

    // Level-set function initial
    double fun_levelSet(double *P, const int i) {
        return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon;
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(pi * y) * sin(pi * x) * (2 * cos(pi * t) * cos(pi * t) - 1) * ((7 * sin(pi * t)) / 25 - x + 0.5)) /
                (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2))) -
            (pi * cos(pi * x) * sin(pi * y) * (2 * cos(pi * t) * cos(pi * t) - 1) * (y + (7 * cos(pi * t)) / 25 - 0.5)) /
                (250 * sqrt(pow(y + (7 * cos(t * pi)) / 25 - 0.5, 2) + pow((7 * sin(t * pi)) / 25 - x + 0.5, 2)));
    }

    // Velocity field
    R fun_velocity(double *P, const int i) {
        if (i == 0)
            return M_PI * (0.5 - P[1]);
        else
            return M_PI * (P[0] - 0.5);
    }

    // Initial solution bulk
    R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

    // Exact solution bulk
    R fun_uBulk(double *P, const int i, const R t) {
        return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
    }

    // Initial solution surface
    R fun_uSurfInit(double *P, const int i) {
        double x = P[0], y = P[1];

        return (0.5 + 0.4 * cos(pi * x) * cos(pi * y) +
            pi / 250 * sin(pi * x) * cos(pi * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) +
            pi / 250 * cos(pi * x) * sin(pi * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)))
            /(1.5 + 0.4 * cos(pi * x) * cos(pi * y));
    }

    // Exact solution surface
    R fun_uSurf(double *P, const int i, const R t) {
        double x = P[0], y = P[1];

        R xc = 0.5 + 0.28 * sin(pi * t), yc = 0.5 - 0.28 * cos(pi * t);

        return (0.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t) +
            pi / 250 * sin(pi * x) * cos(pi * y) * cos(2 * pi * t) * (x - xc) /
                sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) +
            pi / 250 * cos(pi * x) * sin(pi * y) * cos(2 * pi * t) * (y - yc) /
                sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)))
                /(1.5 + 0.4 * cos(pi * x) * cos(pi * y) * cos(2 * pi * t));
    }

    R fun_uBulkD(double *P, const int i, const int d, const R t) {
        return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
    }

    // RHS fB bulk
    R fun_rhsBulk(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * pi * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 125 -
            (4 * pi * cos(pi * x) * cos(pi * y) * sin(2 * pi * t)) / 5 -
            (2 * pi * pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (x - 1. / 2)) / 5 +
            (2 * pi * pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (y - 1. / 2)) / 5;
    }

    /*
    Errors Bulk = [0.0211878, 0.00658872]
    Errors Surface = [0.00419287, 0.00178343]
    */

    // RHS fS surface
    R fun_rhsSurf(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        // return -(pi*(17578125*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 17578125*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 8268750*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 24806250*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 4630500*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)+ 70312500*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 4687500*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 4687500*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 35156250*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 35156250*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 37500000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5000000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 210937500*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 70312500*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 70312500*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 4687500000*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 9375000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 105468750*x*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 70312500*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 105468750*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 19687500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 29531250*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 29531250*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 19687500*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 6250000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1000000*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 6250000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 1000000*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 4687500*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 4687500*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 105468750*x*x*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 316406250*x*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*pow(x,3)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 316406250*y*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 105468750*y*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*pow(y,3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 24806250*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 8268750*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 4630500*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 10500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 2800000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 118125000*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 39375000*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 2500000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2500000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 59062500*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 59062500*x*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 177187500*y*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 59062500*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 2625000000*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 10500000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 9375000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 59062500*x*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 177187500*x*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 59062500*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 59062500*y*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 140625000*x*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 140625000*y*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 16537500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 16537500*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 253125000*x*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 168750000*pow(x,3)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 35000000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 4000000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 84375000*x*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 84375000*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 22500000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 22500000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2000000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 2000000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 253125000*y*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 168750000*pow(y,3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 35000000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 4000000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 11025000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 9843750*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 70312500*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*x*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 70312500*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 210937500*y*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 140625000*pow(y,3)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 6615000*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 6300000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 560000*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 19845000*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 3704400*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 9800000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1120000*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 75000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 11025000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 5512500*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 16537500*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 19687500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 19687500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 9843750*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 3087000*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 75000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 70312500*x*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 140625000*x*x*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 421875000*x*x*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*pow(x,3)*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 70312500*y*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 421875000*y*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 140625000*y*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*pow(y,3)*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 19845000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 3704400*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 9800000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 1120000*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 6615000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 6300000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 560000*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi)+ 18750000000*x*x*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 18750000000*y*y*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5880000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 700000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 18750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 75000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 37500000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 10000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 16537500*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5512500*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 33075000*pi*pow(cos(t*pi),2)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11025000*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 9843750*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 19687500*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 19687500*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 3087000*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),3)*sin(y*pi)+ 6174000*pi*pow(cos(t*pi),3)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 18750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 75000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 10000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 210937500*x*y*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*x*x*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 1470000000*pow(cos(t*pi),2)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 5880000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5880000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 2625000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 784000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 2800000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 700000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 11025000*pi*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 33075000*pi*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi)+ 9843750*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 6174000*pi*cos(y*pi)*pow(sin(t*pi),3)*sin(2*t*pi)*sin(x*pi)+ 1250000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1250000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2500000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 16537500*x*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 59062500*x*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 49612500*y*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 177187500*y*y*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 1470000000*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5880000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 2940000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 2625000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 2800000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3) 
        // - 1400000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 49612500*x*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 177187500*x*x*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 16537500*y*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 59062500*y*y*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 4630500*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 4630500*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 67500000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 45000000*pow(x,3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 6000000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 4000000*pow(x,3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 22500000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 22500000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 2000000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 2000000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 67500000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 45000000*pow(y,3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 6000000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 4000000*pow(y,3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 1543500*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 70312500*x*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 210937500*x*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(x,3)*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 210937500*y*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 70312500*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*pow(y,3)*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 1764000*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 156800*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi) 
        // - 5292000*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 987840*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 470400*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 87808*pow(cos(t*pi),3)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 5512500*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 16537500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 1543500*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3) 
        // - 3087000*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 210937500*x*x*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 281250000*pow(x,3)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 140625000*pow(x,4)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 210937500*y*y*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 281250000*pow(y,3)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 140625000*pow(y,4)*pi*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 18750000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 18750000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 5292000*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 987840*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 470400*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 87808*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 1764000*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 156800*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 2940000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 823200*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 392000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 109760*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 28125000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 28125000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 18750000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10000000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 5512500*pi*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 16537500*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 5512500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 3087000*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)+ 1543500*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 28125000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 28125000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 18750000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 18750000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10000000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 4687500000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2205000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 735000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 2940000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 411600*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 823200*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3) 
        // - 784000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 392000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 109760*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),3) 
        // - 1250000000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5512500*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 1543500*pi*pi*cos(2*t*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 140625000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 140625000*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2343750000*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2343750000*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 93750000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 84375000*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 84375000*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 93750000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 735000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 2205000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 411600*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 784000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 39375000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 140625000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 70312500*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 70312500*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 140625000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 4687500000*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 23625000*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 26250000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 10000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 19687500*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 39375000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 75000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 421875000*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 140625000*x*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 281250000*x*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 421875000*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 281250000*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 140625000*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 26250000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 23625000*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 9375000000*x*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*x*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*y*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 118125000*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 78750000*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 39375000*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 39375000*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 19687500*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 37500000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 75000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 210937500*x*y*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 210937500*x*y*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 2625000000*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 45000000*x*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 4000000*x*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 45000000*x*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 4000000*x*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 59062500*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 3087000*x*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 39375000*pow(x,3)*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 11025000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 19687500*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 140625000*x*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 140625000*x*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 8820000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),2)*sin(x*pi)*pow(sin(y*pi),2)+ 1646400*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),3)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 5880000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi)+ 625000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 625000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 3528000*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 313600*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 12600000*x*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1120000*x*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 10584000*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 37800000*y*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 940800*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi)+ 3360000*y*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 11025000*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 19687500*x*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 11025000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2) 
        // - 33075000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 118125000*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 59062500*y*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 3087000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3) 
        // - 39375000*pow(y,3)*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 37500000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 75000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 75000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 37500000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10584000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 37800000*x*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 940800*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 3360000*x*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 12600000*y*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 1120000*y*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 3528000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 313600*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi)+ 5880000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 31500000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1646400*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 21000000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 4200000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 219520*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 2800000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 1543500*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 1543500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 5880000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 10500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 784000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 1400000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 33075000*x*pi*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 11025000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 19687500*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 118125000*x*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 56250000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 56250000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 56250000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 56250000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 37500000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 37500000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 20000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 20000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 11025000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 33075000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 177187500*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 39375000*y*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 3087000*y*pi*pi*pow(cos(t*pi),3)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 118125000*pow(y,3)*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 21000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 10500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 9375000000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 987840*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 87808*pow(cos(t*pi),2)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 987840*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 87808*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5000000000*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5000000000*y*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 4410000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 1470000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 5880000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 15750000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 10500000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 823200*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 10500000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 1568000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 784000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 5600000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3)+ 1400000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 3087000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 3087000*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 4410000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 1470000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 5880000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 15750000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 10500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 31500000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 1646400*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3) 
        // - 21000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 784000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 4200000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 219520*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),3) 
        // - 2800000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 2500000000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 177187500*x*x*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 3087000*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 118125000*pow(x,3)*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 11025000*y*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 19687500*y*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 7500000*x*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 10500000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 7500000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 7500000*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 10000000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 5000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2625000000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 21000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 9375000000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 7500000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 4687500000*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 168750000*x*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 823200*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 823200*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 109760*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 109760*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 168750000*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 700000000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1470000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 4410000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 10500000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 15750000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1543500*pi*pi*cos(t*pi)*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 1543500*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 1470000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 4410000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 15750000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 823200*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 10500000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 5600000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 196000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 784000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 78750000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 588000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 109760*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2500000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 140625000*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 140625000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 1312500000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1312500000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2625000000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*x*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 47250000*x*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 47250000*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 47250000*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 141750000*y*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 700000000*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 411600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 411600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 439040*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(t*pi),2)*pow(sin(y*pi),3) 
        // - 439040*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 39375000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 78750000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 588000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 109760*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*pow(sin(x*pi),2)+ 700000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 39375000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 118125000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 196000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 784000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 5000000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 843750000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 281250000*x*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 281250000*x*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 2625000000*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1312500000*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1312500000*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5000000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 141750000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 47250000*x*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 47250000*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 47250000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 42000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 5600000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 236250000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*x*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*x*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 118125000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 1400000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 700000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 236250000*y*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 39375000*y*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 2625000000*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 13230000*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 13230000*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 10500000000*y*cos(t*pi)*cos(x*pi)*cos(y*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 11025000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 21000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 42000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 78750000*x*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 236250000*x*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1400000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 236250000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 78750000*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 78750000*y*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 37500000*x*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 10000000*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 56250000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 56250000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 112500000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 112500000*y*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 75000000*pow(x,3)*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 20000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 20000000*y*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 37500000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5000000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 118125000*x*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 10500000000*x*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 21000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 11025000*x*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 66150000*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 22050000*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 22050000*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 33075000*y*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 8820000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 15750000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1646400*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5600000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 2940000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 2800000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 1400000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 9375000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 75000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 75000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 9375000000*y*y*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 5000000000*x*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 118125000*x*y*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 5000000000*y*y*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 33075000*x*x*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 11025000*y*y*pi*pi*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 2940000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 21000000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 2800000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 735000000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 8820000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 5250000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 15750000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 1646400*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),3)*sin(y*pi)+ 5600000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 2800000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 14062500000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 4687500000*x*x*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 4687500000*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 14062500000*y*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*pow(y,3)*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*x*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 392000000*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 33075000*y*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 168750000*x*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 45000000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 4000000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 168750000*x*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 45000000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 4000000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 19687500*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2500000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 140625000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 140625000*x*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 140625000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)+ 140625000*x*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 5250000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 735000000*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 367500000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 367500000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 392000000*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 13230000*x*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 12600000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1120000*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 47250000*x*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 12600000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 12600000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 1120000*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 1120000*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 39690000*y*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 141750000*y*y*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 37800000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 3360000*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(y*pi) 
        // - 11025000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 19687500*x*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 39375000*x*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 196000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 33075000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 118125000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 39375000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 39375000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 118125000*y*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi) 
        // - 5000000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 281250000*x*y*y*pi*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 281250000*x*x*y*pi*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 367500000*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 367500000*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5000000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 39690000*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 141750000*x*x*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 37800000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 3360000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi) 
        // - 47250000*y*y*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 12600000*x*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 12600000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 1120000*x*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi)+ 1120000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 13230000*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 12600000*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 1120000*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 21000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 2800000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 10500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1400000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 33075000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 22050000*x*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(x*pi)+ 19687500*x*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 39375000*x*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 118125000*x*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 78750000*x*x*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 118125000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 56250000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 56250000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 392000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 11025000*y*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)+ 66150000*y*pi*pow(cos(t*pi),2)*cos(x*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 78750000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)+ 39375000*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 39375000*y*pi*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 236250000*y*y*pi*cos(t*pi)*cos(x*pi)*sin(2*t*pi)*sin(y*pi)+ 39375000*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 3704400*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 3528000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 313600*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 3704400*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 3528000*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 313600*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 15750000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 5250000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 10500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(sin(y*pi),3) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 1400000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 3087000*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)+ 3087000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 11025000*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 15750000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 10500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 21000000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 5600000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3) 
        // - 2800000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 66150000*x*pi*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 78750000*x*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 236250000*x*x*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 392000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 22050000*y*pi*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi)+ 19687500*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*y*y*pi*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 37500000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5000000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 10000000*x*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 5000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 10000000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 5000000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 20000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 15000000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 20000000*y*y*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 15000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 10000000*pow(x,3)*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 10000000*pow(y,3)*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 18750000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 37500000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5000000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 56250000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 10500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 15750000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 6174000*pi*cos(t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 6174000*pi*pow(cos(t*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 5250000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 15750000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3) 
        // - 4410000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1568000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 1176000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 700000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 2800000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 219520*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 112500000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 75000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 4410000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 823200*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 784000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 392000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 1400000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 75000000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 112500000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1470000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 1470000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3) 
        // - 5880000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 4410000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 823200*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(x*pi) 
        // - 784000*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 392000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi)+ 1400000*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 8820000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 4410000*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 1646400*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1568000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 1176000*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 700000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 2800000*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 219520*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(y*pi) 
        // - 44100000*x*pi*cos(t*pi)*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 44100000*y*pi*cos(t*pi)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 31500000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*x*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 2800000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 5880000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 31500000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 17640000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 31500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 63000000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 2800000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 11200000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 8400000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 150000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 31500000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 5600000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 2800000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 98000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 150000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 98000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 31500000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 5600000*x*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 2800000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 42000000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 31500000*x*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 5880000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 42000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 31500000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 63000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 11200000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 8400000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 2800000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 63000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 31500000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 21000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 2800000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 9375000000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 9375000000*x*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 56250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 75000000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 37500000*pow(x,4)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 56250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 37500000*pow(y,4)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 39375000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)+ 1646400*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 8820000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 784000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 63000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 1646400*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 8820000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 784000*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 735000000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 42000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 735000000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*y*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 1470000*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 25200000*x*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2240000*x*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(x*pi)+ 78750000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)+ 39375000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 11760000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 11760000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 735000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 5250000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 735000000*y*pi*cos(2*t*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1470000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 25200000*x*y*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 2240000*x*y*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 22050000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 21000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)+ 2800000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)+ 22050000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 39375000*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 78750000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi) 
        // - 7056000*x*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 7056000*y*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 627200*x*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),4)*pow(cos(y*pi),3)*sin(t*pi)*sin(y*pi) 
        // - 627200*y*cos(t*pi)*pow(cos(2*t*pi),4)*pow(cos(x*pi),3)*pow(cos(y*pi),4)*sin(t*pi)*sin(x*pi)+ 1250000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 22050000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 31500000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2) 
        // - 21000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 21000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 11200000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*pow(sin(x*pi),3)+ 2800000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 22050000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 39375000*x*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 75000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 112500000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 15000000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 10000000*x*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 37500000*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 37500000*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 112500000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 15000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 20000000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*x*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 75000000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 112500000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 112500000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 11760000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 1568000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 11760000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 1568000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 21000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)+ 31500000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 11200000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 8820000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 1400000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2800000*x*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 10500000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 31500000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 784000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 4200000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5600000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 2352000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 2800000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 5600000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 8400000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 21000000*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 150000000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 31500000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 4200000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 150000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 31500000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 5880000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 3136000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*sin(t*pi)*pow(sin(y*pi),3)+ 5880000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3136000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),3)+ 11760000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)+ 31500000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 4200000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 21000000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 42000000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 31500000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 2352000*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 2800000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5600000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 784000*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 4200000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 8400000*x*x*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 10500000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 17640000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 8820000*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 63000000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(x*pi),2)*sin(y*pi)+ 1400000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5600000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 2800000*y*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 31500000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 22050000*x*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 39375000*x*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi) 
        // - 2469600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 219520*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*sin(x*pi) 
        // - 392000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(sin(t*pi),2)*sin(x*pi)*pow(sin(y*pi),2)+ 63000000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 2469600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)+ 219520*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 392000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 1568000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 42000000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 11760000*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 1250000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1250000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 22050000*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2) 
        // - 39375000*x*x*y*pi*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 3292800*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2)+ 3292800*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 75000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 350000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 75000000*x*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 350000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 3087000*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)+ 22050000*x*x*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 11760000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 21000000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2) 
        // - 1568000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 2800000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3) 
        // - 3087000*y*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi) 
        // - 22050000*y*y*pi*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)+ 39375000*x*x*y*pi*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi) 
        // - 21000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 63000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 21000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 350000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 350000000*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 8820000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 2940000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)+ 11760000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2) 
        // - 31500000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(y*pi),2)+ 21000000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2) 
        // - 21000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)+ 1568000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2) 
        // - 2800000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi)+ 39375000*x*y*y*pi*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 21000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 112500000*x*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 75000000*pow(x,3)*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 15000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 20000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 15000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 10000000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 63000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 21000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 37500000*x*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 37500000*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 15000000*x*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 15000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 20000000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 10000000*x*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 10000000*pow(x,3)*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5250000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 112500000*y*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 75000000*pow(y,3)*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 9375000000*x*y*pi*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 9375000000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 18750000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 2800000000*y*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1646400*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)+ 11760000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 219520*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*pow(sin(t*pi),2)+ 1568000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 1646400*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 11760000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi) 
        // - 219520*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 1568000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),3)*sin(t*pi) 
        // - 3087000*x*pi*pi*pow(cos(t*pi),2)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 2940000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 8820000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 21000000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 31500000*x*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3087000*y*pi*pi*cos(t*pi)*cos(2*t*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 392000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 1568000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 5880000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 2940000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 1176000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 392000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 4200000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 5600000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 2800000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2) 
        // - 219520*x*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2800000*pow(x,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 5880000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 8820000*pi*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 1646400*pi*pow(cos(t*pi),3)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 2625000000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*x*pi*cos(2*t*pi)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1176000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 4200000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 2625000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 18750000000*x*y*pi*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2625000*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 2800000000*x*cos(2*t*pi)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 94500000*x*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 823200*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 5880000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 823200*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 5880000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1176000*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2)+ 4200000*x*x*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 78750000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi) 
        // - 8820000*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 1646400*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),3)*sin(2*t*pi)*sin(x*pi)+ 392000*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 1176000*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 2800000*x*x*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 4200000*y*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 5600000*y*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 219520*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),3)*pow(sin(x*pi),2) 
        // - 2800000*pow(y,3)*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 1400000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 5250000000*x*pi*cos(t*pi)*cos(2*t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*x*pi*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2625000000*x*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 392000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 1568000*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi) 
        // - 2625000000*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2625000*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 94500000*x*y*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 109760*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 439040*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 109760*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 439040*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(y*pi)+ 1400000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*pow(cos(y*pi),2)*sin(t*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 157500000*x*y*pi*cos(t*pi)*cos(y*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 78750000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi) 
        // - 2800000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(y*pi),2)*pow(sin(x*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000000*y*pi*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1250000000*x*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 3750000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 3750000000*y*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*pow(x,3)*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*pow(y,3)*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1250000000*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 26460000*x*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 26460000*y*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi) 
        // - 22050000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)+ 22050000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(t*pi)*sin(y*pi)+ 2800000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*sin(t*pi)*pow(sin(y*pi),2)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 157500000*x*y*pi*cos(x*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 1470000000*pi*cos(t*pi)*cos(2*t*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 112500000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)+ 75000000*x*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 10000000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi)+ 98000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 112500000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 75000000*x*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)+ 10000000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi) 
        // - 98000000*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2940000*x*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 17640000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 8820000*y*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 2500000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*x*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 8820000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 700000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2940000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 700000000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi) 
        // - 700000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 700000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 75000000*x*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 75000000*x*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11760000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 5880000*x*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 784000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 5600000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 21000000*x*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 2352000*x*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 8400000*x*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2)+ 11760000*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)+ 17640000*y*pi*pow(cos(t*pi),2)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 63000000*y*y*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 5250000000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 5250000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 21000000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 17640000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 63000000*x*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 2352000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(t*pi),2)*pow(sin(x*pi),2) 
        // - 8400000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 21000000*y*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 784000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2)+ 5600000*x*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2)+ 5880000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi)+ 5250000000*x*y*pi*cos(2*t*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 5250000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 1646400*pi*pow(cos(t*pi),2)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 219520*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(t*pi),2)*pow(sin(y*pi),2) 
        // - 219520*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1568000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 1568000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 1646400*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(2*t*pi)*sin(y*pi) 
        // - 2500000000*x*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 2500000000*x*y*y*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 196000000*x*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*y*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*y*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 63000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(x*pi) 
        // - 5600000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(y*pi)+ 196000000*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*pow(sin(t*pi),2)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 196000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*pow(sin(t*pi),2)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 1400000000*x*x*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 84000000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 5600000*x*y*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 63000000*x*y*pi*pi*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 1568000*x*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi) 
        // - 17640000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(y*pi) 
        // - 17640000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)+ 1568000*y*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi) 
        // - 84000000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 2940000*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 5250000*x*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 8820000*y*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 47250000*y*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*y*pi*pi*pow(cos(t*pi),3)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 31500000*pow(y,3)*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 23520000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(y*pi)*sin(t*pi)*pow(sin(x*pi),2)*sin(y*pi) 
        // - 23520000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*sin(t*pi)*sin(x*pi)*pow(sin(y*pi),2) 
        // - 8820000*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 47250000*x*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*x*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),3)*sin(x*pi)*sin(y*pi) 
        // - 31500000*pow(x,3)*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 2940000*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi) 
        // - 5250000*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 411600*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 2500000000*x*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 2500000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 75000000*x*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 75000000*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 225000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 700000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 21000000*x*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 5600000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*pow(sin(x*pi),2)+ 11200000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(x*pi) 
        // - 21000000*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 21000000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi)+ 8400000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*pow(sin(y*pi),2) 
        // - 63000000*y*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(2*t*pi)*sin(y*pi)+ 63000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 700000000*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 63000000*x*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 8400000*x*y*pi*pi*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2)+ 21000000*x*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 21000000*y*pi*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 5600000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 11200000*x*y*pi*pi*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 21000000*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 63000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)+ 5880000*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi) 
        // - 1568000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*cos(y*pi)*sin(t*pi)*pow(sin(y*pi),2) 
        // - 3136000*x*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),3)*pow(cos(y*pi),2)*sin(t*pi)*sin(y*pi)+ 1568000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*cos(x*pi)*pow(cos(y*pi),3)*sin(t*pi)*pow(sin(x*pi),2) 
        // - 3136000*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),3)*pow(cos(x*pi),2)*pow(cos(y*pi),3)*sin(t*pi)*sin(x*pi)+ 5880000*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi)+ 10500000*x*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)+ 10500000*x*y*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*x*pi*pi*pow(cos(t*pi),2)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi) 
        // - 823200*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*pow(sin(t*pi),2)*sin(x*pi)*sin(y*pi)+ 5000000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 42000000*x*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(2*t*pi)*sin(x*pi) 
        // - 1400000000*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 42000000*x*y*pi*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 392000000*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2))) 
        // - 11760000*x*pi*cos(t*pi)*cos(2*t*pi)*pow(cos(x*pi),2)*cos(y*pi)*sin(t*pi)*sin(2*t*pi)*sin(y*pi) 
        // - 11760000*y*pi*cos(t*pi)*cos(2*t*pi)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(2*t*pi)*sin(x*pi)+ 1400000000*x*y*pi*cos(t*pi)*pow(cos(2*t*pi),2)*pow(cos(x*pi),2)*cos(y*pi)*sin(y*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 
        // - x+0.5),2)))+ 1400000000*x*y*pi*pow(cos(2*t*pi),2)*cos(x*pi)*pow(cos(y*pi),2)*sin(t*pi)*sin(x*pi)*sqrt((pow((y+ (7*cos(t*pi))/25 
        // -0.5),2)+ pow(((7*sin(t*pi))/25 - x+0.5),2))) 
        // - 10500000*x*y*pi*pi*cos(t*pi)*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(x*pi)*sin(y*pi) 
        // - 10500000*x*y*pi*pi*pow(cos(2*t*pi),2)*cos(x*pi)*cos(y*pi)*sin(t*pi)*sin(x*pi)*sin(y*pi)))
        // /(12500*sqrt(pow((y+ (7*cos(pi*t))/25 -0.5),2)+ pow(((7*sin(pi*t))/25 - x+0.5),2))
        // *pow((4*cos(2*pi*t)*cos(pi*x)*cos(pi*y)+ 15),3)*(350*sin(t*pi) - 1250*y- 350*cos(t*pi) 
        // - 1250*x+ 700*y*cos(t*pi) - 700*x*sin(t*pi)+ 1250*x*x+ 1250*y*y+ 98
        // *pow(cos(t*pi),2)+ 98*pow(sin(t*pi),2)+ 625));


        return -(pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                 1) *
                   ((pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     1) *
                        (-((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                               (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                            pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                           2.5)) +
                           (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                            pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                           2.5)) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (2 * y + (14 * cos(t * pi)) / 25 - 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                         (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                         (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(x * pi) * sin(y * pi) *
                          sin(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                         (4 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                          (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)))) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
                    ((2 * y + (14 * cos(t * pi)) / 25 - 1) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     ((2 * y + (14 * cos(t * pi)) / 25 - 1) * pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2)) /
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2)) *
                        ((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                           (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                         (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
                    ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                    ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)))) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                    ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                            2)) -
               ((4 * pi * cos(x * pi) * cos(y * pi) * sin(2 * t * pi)) / 5 -
                (pi * pi * cos(y * pi) * sin(2 * t * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                    (125 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (7 * pi * pi * cos(t * pi) * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                    (6250 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (7 * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(t * pi) * sin(y * pi)) /
                    (6250 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                (pi * pi * cos(x * pi) * sin(2 * t * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                    (125 *
                     sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                 ((14 * pi * sin(t * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) / 25 -
                  (14 * pi * cos(t * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) / 25) *
                 (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                    (500 *
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         1.5)) +
                (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                 ((14 * pi * sin(t * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) / 25 -
                  (14 * pi * cos(t * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) / 25) *
                 (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                    (500 *
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         1.5))) /
                   ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
               (pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                    (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                1) *
                   ((pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                         (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                     1) *
                        (((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                              (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                          (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                           (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                          (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                           (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) -
                          (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                           pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          2.5)) +
                          (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                           (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                           (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                              (250 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                          pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                         1.5)) +
                          (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                           pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2)) /
                              (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          2.5))) /
                             ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                         (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                         (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(y * pi) * cos(y * pi) * sin(x * pi) *
                          sin(x * pi) *
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           0.5)) /
                             (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                         (4 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                          ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                            (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                            pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                           (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                            (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)) -
                           (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                            (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                               (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                           pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                          1.5)))) /
                             (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) -
                    (((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                       (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                          (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                     1.5)) -
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                       (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                     1.5))) /
                         ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                     (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       0.5)) /
                         (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                        ((-2 * x + (14 * sin(t * pi)) / 25 + 1) / (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) +
                                                                   pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                         (pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                             pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                  pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                 2)) -
                    (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                    ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)))) /
                        (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                    (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                            2)) +
               pi * (x - 0.5) *
                   ((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                        ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                    (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                     ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      0.5)) /
                        (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) +
               pi *
                   (((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                      (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                        ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                    (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                     ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                      (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                      (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                          (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                      0.5)) /
                        (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                   (y - 0.5) +
               ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                ((pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) /
                      (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                  1) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3))) -
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-2 * x + (14 * sin(t * pi)) / 25 + 1) * pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)), 2) +
                 ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                  (((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                        (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                    (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                     pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (1000 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2.5)) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (250 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) +
                    (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     pow((-2 * x + (14 * sin(t * pi)) / 25 + 1), 2)) /
                        (1000 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             2.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                   (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(y * pi) * cos(y * pi) * sin(x * pi) *
                    sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                   (4 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                      (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)))) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)))) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         2))) /
                   (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
               ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                ((pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2) /
                      (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                  1) *
                     (((2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) / 5 -
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi)) /
                           (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                       (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) -
                       (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                       (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                       pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                      1.5)) +
                       (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-x + (7 * sin(t * pi)) / 25 + 0.5) * (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5)) -
                       (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                        (-2 * x + (14 * sin(t * pi)) / 25 + 1) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                           (1000 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       2.5))) /
                          ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                      (2 * pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                       (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) -
                      (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                       ((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                         (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)) -
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                         (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (500 * pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                        pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                                       1.5)))) /
                          (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                      (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * sin(x * pi) *
                       sin(y * pi) *
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                        (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                        (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                            (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                         pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                        0.5)) /
                          (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3))) -
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) -
                 ((-x + (7 * sin(t * pi)) / 25 + 0.5) * (y + (7 * cos(t * pi)) / 25 - 0.5) *
                  (-((2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 -
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) /
                         (125 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * pi * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * pi * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (3 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                      pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (1000 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              2.5)) +
                     (3 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                      pow((2 * y + (14 * cos(t * pi)) / 25 - 1), 2) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (1000 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              2.5)) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)) +
                   (8 * pi * pi * cos(2 * t * pi) * cos(2 * t * pi) * cos(x * pi) * cos(x * pi) * sin(y * pi) *
                    sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (25 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 3)) +
                   (4 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    (-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)) -
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                      (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (500 *
                          pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                              1.5)))) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2)))) /
                     (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
                 ((((2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) / 5 -
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                     (-2 * x + (14 * sin(t * pi)) / 25 + 1)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (-2 * x + (14 * sin(t * pi)) / 25 + 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) -
                   (2 * pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (2 * y + (14 * cos(t * pi)) / 25 - 1) * pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)), 2) +
                 (((-(2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) / 5 +
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * sin(x * pi) * sin(y * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * pi * cos(2 * t * pi) * cos(x * pi) * cos(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                     pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                    (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5)) -
                    (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (2 * y + (14 * cos(t * pi)) / 25 - 1) *
                     (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                        (500 *
                         pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                             1.5))) /
                       ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 + 1.5) +
                   (2 * pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) *
                    ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                     (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                     (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                         (250 * sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) +
                                      pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                     0.5)) /
                       (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2))) *
                  (2 * y + (14 * cos(t * pi)) / 25 - 1) * (-x + (7 * sin(t * pi)) / 25 + 0.5) *
                  (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     pow((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)),
                         2))) /
                   (pow((y + (7 * cos(t * pi)) / 25 - 0.5), 2) + pow((-x + (7 * sin(t * pi)) / 25 + 0.5), 2)) +
               (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                   (250 *
                    sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
               (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                   (250 *
                    sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
               (4 * pi * cos(x * pi) * cos(y * pi) * sin(2 * t * pi) *
                ((2 * cos(2 * t * pi) * cos(x * pi) * cos(y * pi)) / 5 +
                 (pi * cos(2 * t * pi) * cos(x * pi) * sin(y * pi) * (y + (7 * cos(t * pi)) / 25 - 0.5)) /
                     (250 *
                      sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) -
                 (pi * cos(2 * t * pi) * cos(y * pi) * sin(x * pi) * (-x + (7 * sin(t * pi)) / 25 + 0.5)) /
                     (250 *
                      sqrt((pow(((7 * sin(pi * t)) / 25 - x + 0.5), 2) + pow((y + (7 * cos(pi * t)) / 25 - 0.5), 2)))) +
                 0.5)) /
                   (5 * pow(((2 * cos(2 * pi * t) * cos(pi * x) * cos(pi * y)) / 5 + 1.5), 2));
    }

    R fun_neumann_left(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
    }

    R fun_neumann_bottom(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
    }

    R fun_neumann_right(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x)) / 250;
    }

    R fun_neumann_top(double *P, const int i, const R t) {
        R x = P[0], y = P[1];

        return -(pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y)) / 250;
    }
} 

namespace Example1Omega2 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17);
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// Initial solution surface
R fun_uSurfInit(double *P, const int i) {
    double x = P[0], y = P[1];

    return (0.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) 
    - M_PI / 250 * sin(M_PI * x) * cos(M_PI * y) * (x - 0.5) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)) 
    - M_PI / 250 * cos(M_PI * x) * sin(M_PI * y) * (y - 0.22) / sqrt((y - 0.22) * (y - 0.22) + (x - 0.5) * (x - 0.5)))
    / (1.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y));

}


// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);

    return (0.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t)
    - M_PI / 250 * sin(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t) * (x - xc) / sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)) 
    - M_PI / 250 * cos(M_PI * x) * sin(M_PI * y) * cos(2 * M_PI * t) * (y - yc) / sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc)))
    / (1.5 + 0.4 * cos(M_PI * x) * cos(M_PI * y) * cos(2 * M_PI * t));
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // return M_PI*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*(-4.0/5.0)+((M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI))/1.25E+2-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*(x-1.0/2.0)*(2.0/5.0)+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*(y-1.0/2.0)*(2.0/5.0);

    return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
           (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 1. / 2)) / 5 +
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 1. / 2)) / 5;
}

// RHS fS surface
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (M_PI*1.0/sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.0/pow(cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*4.0+1.5E+1,3.0)*(cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.23046875E+8+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.23046875E+8+cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.26875E+6+cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.480625E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*4.6305E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*7.03125E+7+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.96875E+7+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.96875E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.515625E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.515625E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*3.75E+7+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*5.0E+6+M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8-M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.03125E+7-M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.03125E+7+cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*9.375E+6-x*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.0546875E+8-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.515625E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.515625E+8-y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.0546875E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*9.84375E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.953125E+7+cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.953125E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*9.84375E+7+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.625E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.0E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.625E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.0E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.6875E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.6875E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.5E+6+(x*x)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.0546875E+8+(x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.1640625E+8-(x*x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8+(y*y)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.1640625E+8+(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.0546875E+8-(y*y*y)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.480625E+7+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*8.26875E+6-pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*4.6305E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.05E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8-M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.9375E+7-M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.5E+6-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*9.375E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.5E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.5E+6+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.90625E+7+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.90625E+7+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.771875E+8+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.90625E+7+cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.05E+7-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*9.375E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*9.375E+6-x*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.90625E+7-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.771875E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.90625E+7-y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*5.90625E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.65375E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.65375E+7+(x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.53125E+8-(x*x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.5E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6+(x*x)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.4375E+7+(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.4375E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.25E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.25E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*2.0E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.53125E+8-(y*y*y)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.5E+7-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*9.84375E+6+x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7-(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8+(y*y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.615E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*6.3E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*5.6E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.9845E+7-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.54E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.5125E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.65375E+7-(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*9.84375E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.087E+6+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.03125E+7-(x*x)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.40625E+8-(x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.21875E+8+(x*x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.03125E+7-(y*y)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.21875E+8-(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.40625E+8+(y*y*y)*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.9845E+7+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*3.7044E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.54E+7+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*6.615E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*6.3E+6+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+5-(x*x)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-(y*y)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.25E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.0E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.5E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.0E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.0E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.65375E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.5125E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.3075E+7-M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*9.84375E+6-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.96875E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.96875E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*3.087E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.174E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.5E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.0E+7-x*(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-(x*x)*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-pow(cos(t*M_PI),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.94E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.88E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.25E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.625E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*5.25E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.8E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*7.0E+5-M_PI*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.1025E+7-M_PI*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.3075E+7+(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*9.84375E+6-M_PI*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.174E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.5E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.5E+6-x*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.65375E+7-(x*x)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*5.90625E+7-y*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*4.96125E+7-(y*y)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.771875E+8-cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.88E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.94E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.625E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*2.8E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*1.4E+6-x*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.96125E+7+(x*x)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.771875E+8-y*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.65375E+7+(y*y)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*5.90625E+7-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.6305E+6+pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*4.6305E+6+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*6.75E+7-(x*x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*6.0E+6-(x*x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.25E+7+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.25E+7+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.0E+6+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*6.75E+7-(y*y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*6.0E+6-(y*y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*5.5125E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.5435E+6-(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.764E+6+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.568E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.292E+6-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*9.8784E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.704E+5-pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*8.7808E+4+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*5.5125E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.65375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*1.5435E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.087E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.8125E+8+(x*x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.40625E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.109375E+8+(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.8125E+8-(y*y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.40625E+8-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.875E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.875E+7+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.292E+6+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*9.8784E+5+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.704E+5+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*8.7808E+4+pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.764E+6+pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.568E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.94E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*8.232E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0976E+5-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.8125E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.8125E+7+(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*1.0E+7-(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.65375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*3.087E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.8125E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.8125E+7+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.875E+7-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.875E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.0E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.205E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.35E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*2.94E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*8.232E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*7.84E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*1.0976E+5+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*5.5125E+6+(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.40625E+8+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.34375E+9-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.34375E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.4375E+8-x*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.4375E+7-y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.4375E+7-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.4375E+8-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.35E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.205E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*4.116E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*7.84E+5-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.03125E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.03125E+7+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.3625E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*6.825E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.96875E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*7.5E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.0E+7-x*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*4.21875E+8+x*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.40625E+8+x*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*4.21875E+8+y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+y*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.40625E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*6.825E+7+pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.3625E+7-x*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+y*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.5E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.96875E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.5E+7+x*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*2.109375E+8+x*y*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*2.109375E+8-cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7-x*(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6-(x*x)*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7-(x*x)*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*5.90625E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.087E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.1025E+7+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.96875E+7+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*8.82E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.6464E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*6.25E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*6.25E+8-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*3.528E+6-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*3.136E+5-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.26E+7-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0584E+7-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.78E+7-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*9.408E+5-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*3.36E+6-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.1025E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.96875E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.1025E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.3075E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.18125E+8+(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*5.90625E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*3.087E+6-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.0584E+7+(x*x)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.78E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*9.408E+5+(x*x)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.36E+6+(y*y)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.26E+7+(y*y)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.528E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.136E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*3.15E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.6464E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*4.2E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.1952E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*1.5435E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.5435E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*5.88E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.05E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*7.84E+5+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.4E+6-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.3075E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.18125E+8+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7-x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7-x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7-(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.75E+7+(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*3.75E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.0E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.0E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.1025E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.771875E+8-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.05E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*9.8784E+5+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*8.7808E+4-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*9.8784E+5-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.7808E+4+x*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+y*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.41E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.47E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.232E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.05E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.6E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.087E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.087E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.41E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.47E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*5.88E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.05E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*3.15E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*1.6464E+6-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*7.84E+5+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*4.2E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*2.1952E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.771875E+8-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.18125E+8+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.5E+6-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.05E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.5E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+6+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*8.232E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*8.232E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*1.0976E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.0976E+5+x*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.47E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.41E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.05E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.5435E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.47E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.41E+6+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*8.232E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.05E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*1.568E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*5.6E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.96E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*7.84E+5+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*7.875E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+5+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0976E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.725E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.725E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.725E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.4175E+8+pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.116E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),3.0)*4.3904E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*4.3904E+5+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*7.875E+7+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.88E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0976E+5-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.96E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*7.84E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*8.4375E+8-x*y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8-x*y*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.3125E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.4175E+8-x*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*4.725E+7-y*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.725E+7-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*4.725E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*4.2E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*5.6E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*2.3625E+8-x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7-x*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-y*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.3625E+8-y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.323E+7-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.323E+7-y*cos(t*M_PI)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.05E+10-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*4.2E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*5.6E+6-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.1025E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.1025E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.1E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*4.2E+7+x*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+x*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.3625E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-y*M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.3625E+8+y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+y*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.0E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.0E+7+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6-x*y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8+x*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.05E+10-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.1E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.615E+7+M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.205E+7+M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.205E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.82E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.575E+7+M_PI*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6464E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.875E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.94E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.4E+6+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.875E+7-(x*x)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*y*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8-(y*y)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+(x*x)*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.3075E+7-(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.1025E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.94E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.8E+6-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.4E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*5.25E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.575E+7+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*1.6464E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.40625E+10-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*4.6875E+9-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.40625E+10+(y*y*y)*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+x*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.3075E+7-pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+y*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.3075E+7-x*(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.6875E+8+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*4.5E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*4.0E+6-(x*x)*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.6875E+8+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*4.5E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*4.0E+6-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*1.96875E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8-x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.40625E+8+(x*x)*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.40625E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*5.25E+6+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*2.1E+7-cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8-x*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.323E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.26E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.12E+6-(x*x)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.725E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.26E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.26E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*1.12E+6+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.969E+7-(y*y)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.4175E+8+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.78E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(y*M_PI)*3.36E+6-x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.1025E+7-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*1.96875E+7+(x*x)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.3075E+7+y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*1.18125E+8-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*3.9375E+7-y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*1.18125E+8+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9+x*(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.8125E+8+(x*x)*y*M_PI*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.8125E+8+M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8-M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.675E+8+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.969E+7+(x*x)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.4175E+8-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.78E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.36E+6+(y*y)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.725E+7-x*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.26E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.26E+7-x*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*1.12E+6-y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.323E+7-y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.26E+7-y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.05E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.4E+6-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.3075E+7+x*M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.205E+7+x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7+x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.18125E+8+(x*x)*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.875E+7+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.18125E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*5.625E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.625E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*2.0E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*2.0E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.1025E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.615E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*7.875E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*3.9375E+7+(y*y)*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.3625E+8+(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.9375E+7+pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.528E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*3.136E+5-cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*3.7044E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*3.528E+6-cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*3.136E+5-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*5.25E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.05E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(sin(y*M_PI),3.0)*5.6E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.6E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.4E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.087E+6+M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.087E+6+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.1025E+7+(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*1.1025E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.575E+7+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.05E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*5.6E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+x*M_PI*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.615E+7-x*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.875E+7-(x*x)*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.3625E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8+y*M_PI*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.205E+7+y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.96875E+7-(y*y)*M_PI*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.875E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.75E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6-(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7-(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.0E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.875E+7-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*5.625E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.0E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.5E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.0E+7-(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.5E+7-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.875E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.625E+7-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.75E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.0E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6-(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.0E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*5.6E+6+M_PI*cos(t*M_PI)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.174E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.174E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.25E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.575E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*5.6E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*5.6E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*4.41E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.568E+6-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.176E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.0E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.1952E+5+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.125E+8-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.41E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.232E+5-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*7.84E+5+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.4E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.4E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*7.5E+7+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.125E+8-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*7.5E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.47E+6-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.47E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*1.568E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*1.568E+6-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*4.41E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*8.232E+5-M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*7.84E+5-M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*3.92E+5+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.4E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.4E+6+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.41E+6-(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.6464E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.568E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.176E+6-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*7.0E+5-(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(y*M_PI)*2.1952E+5-x*M_PI*cos(t*M_PI)*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.41E+7-y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.41E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.15E+7+(x*x)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.764E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*3.15E+7+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.8E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.12E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*8.4E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+8+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.8E+6-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.5E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.6E+6+x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*2.8E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.764E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.15E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*5.88E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.12E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*8.4E+6+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*2.8E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*6.3E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*3.15E+7+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6-x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(x*x)*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.5E+7+(x*x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.75E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.625E+7+(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*7.5E+7-(y*y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.75E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*1.6464E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*8.82E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*7.84E+5-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*6.3E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.6464E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*8.82E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*7.84E+5-x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.47E+6+(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.52E+7-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(x*M_PI)*2.24E+6+x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*7.875E+7+x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.176E+7-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.35E+8+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*1.47E+6+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+x*y*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.52E+7+x*y*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*2.24E+6-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6+y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*7.875E+7+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*7.056E+6+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*7.056E+6+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),4.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(y*M_PI)*6.272E+5+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),4.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),4.0)*sin(t*M_PI)*sin(x*M_PI)*6.272E+5-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.205E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.15E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),3.0)*1.12E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.0E+7+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7+x*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.75E+7+y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.75E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.0E+7-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.125E+8-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.125E+8-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.15E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*1.12E+7+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*8.82E+6+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.4E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.6E+6+(x*x)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.8E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.05E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*7.84E+5-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.352E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.8E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*5.6E+6+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*8.4E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+8-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*8.82E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*3.15E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.5E+8-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*sin(t*M_PI)*pow(sin(y*M_PI),3.0)*3.136E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.88E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),3.0)*3.136E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*8.82E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*4.2E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*3.15E+7-x*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*2.352E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+6+y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*7.84E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6+(x*x)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*8.4E+6-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.05E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.764E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*8.82E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*6.3E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.4E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*5.6E+6+(y*y)*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.8E+6+M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*2.205E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*3.9375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.4696E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*2.1952E+5-(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*1.764E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*6.3E+7+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*2.4696E+6+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*2.1952E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*3.92E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.568E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*1.176E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+x*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*2.205E+7-(x*x)*y*(M_PI*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.9375E+7-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*3.2928E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*3.2928E+6-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8-x*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*3.087E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*1.176E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*1.568E+6-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*2.8E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*3.087E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*2.205E+7+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*2.1E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.1E+7-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*2.1E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8+M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.5E+8-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.82E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.94E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*1.176E+7-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.15E+7+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.1E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*1.568E+6-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.8E+6+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.9375E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.1E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.125E+8+(x*x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7-x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*2.0E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.5E+7+x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7+(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*1.0E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.1E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*3.75E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*3.75E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.5E+7+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*2.0E+7-x*(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7-(x*x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.0E+7+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.1E+7-(y*y)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.125E+8+(y*y*y)*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7+x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.375E+9-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.875E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.875E+7-y*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*1.6464E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*2.1952E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.6464E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*1.176E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*2.1952E+5-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*1.568E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.94E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*8.82E+6+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.1E+7-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*3.15E+7-y*(M_PI*M_PI)*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*3.087E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*3.92E+5+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.568E+6+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.94E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*3.92E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.8E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.1952E+5-(x*x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*5.88E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*8.82E+6+M_PI*pow(cos(t*M_PI),3.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.6464E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*1.176E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*4.2E+6+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*y*M_PI*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.875E+10+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.625E+6+x*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*9.45E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*8.232E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.88E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.232E+5-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*5.88E+6-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.176E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*7.875E+7-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*8.82E+6-M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),3.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.6464E+6+x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.92E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*1.176E+6-(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*2.8E+6+(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*4.2E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6-y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),3.0)*pow(sin(x*M_PI),2.0)*2.1952E+5-(y*y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.8E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.94E+6+x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*3.92E+5+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*1.568E+6+y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.625E+9+(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.625E+6+x*y*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*9.45E+7+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.0976E+5-(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*4.3904E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*1.0976E+5+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*4.3904E+5-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*M_PI*cos(t*M_PI)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.575E+8-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*7.875E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(y*M_PI),2.0)*pow(sin(x*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-y*M_PI*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.75E+9-(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.75E+9-(x*x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(y*y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.25E+9+x*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.646E+7+y*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*2.646E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*2.205E+7+y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(y*M_PI)*2.205E+7-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.8E+9-x*y*M_PI*cos(x*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.575E+8-M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.47E+9-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*1.125E+8+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*7.5E+7-x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.0E+7-M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*1.125E+8-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*7.5E+7+x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*1.0E+7+M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*9.8E+7+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.764E+7-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*4.2E+7+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*4.2E+7+x*(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7+(x*x)*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.176E+7+x*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*5.88E+6+x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*7.84E+5+x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+(x*x)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*2.352E+6-x*(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+6+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.176E+7+y*M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.764E+7+(y*y)*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.3E+7-x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9+x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.1E+7+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.764E+7-(x*x)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.3E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(t*M_PI),2.0)*pow(sin(x*M_PI),2.0)*2.352E+6-(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.4E+6-(y*y)*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*7.84E+5+(x*x)*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*5.88E+6-x*y*M_PI*cos(t*M_PI*2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.25E+9-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.1E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6-M_PI*pow(cos(t*M_PI),2.0)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.6464E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*pow(sin(y*M_PI),2.0)*2.1952E+5-y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*2.1952E+5+(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.568E+6-(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.568E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.6464E+6+(x*x)*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*(y*y)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+(y*y)*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(x*M_PI)*6.3E+7-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(y*M_PI)*5.6E+6-x*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.96E+8+(x*x)*M_PI*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+7-x*y*M_PI*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*5.6E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*6.3E+7-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*1.568E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(y*M_PI)*1.764E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*1.764E+7+y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.568E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*8.4E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-(x*x)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6+y*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+(y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.725E+7-y*(M_PI*M_PI)*pow(cos(t*M_PI),3.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-(y*y*y)*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*sin(y*M_PI)*2.352E+7-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*sin(t*M_PI)*sin(x*M_PI)*pow(sin(y*M_PI),2.0)*2.352E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.82E+6+(x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.725E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),3.0)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-(x*x*x)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*3.15E+7+y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*2.94E+6-(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*5.25E+6+(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*4.116E+5+x*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*2.5E+9-x*y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*7.5E+7-x*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*7.5E+7+x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*2.25E+8+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*pow(sin(x*M_PI),2.0)*5.6E+6+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(x*M_PI)*1.12E+7-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7+x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*pow(sin(y*M_PI),2.0)*8.4E+6-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*6.3E+7+x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.3E+7-M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*7.0E+8+x*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*6.3E+7+x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*8.4E+6+x*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7+y*M_PI*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*2.1E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*5.6E+6-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*1.12E+7+y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*2.1E+7-y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*6.3E+7+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*5.88E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*cos(y*M_PI)*sin(t*M_PI)*pow(sin(y*M_PI),2.0)*1.568E+6-x*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),3.0)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(y*M_PI)*3.136E+6+y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*cos(x*M_PI)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*pow(sin(x*M_PI),2.0)*1.568E+6-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),3.0)*pow(cos(x*M_PI),2.0)*pow(cos(y*M_PI),3.0)*sin(t*M_PI)*sin(x*M_PI)*3.136E+6+M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*5.88E+6+(x*x)*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7+x*(y*y)*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7-x*(M_PI*M_PI)*pow(cos(t*M_PI),2.0)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*pow(sin(t*M_PI),2.0)*sin(x*M_PI)*sin(y*M_PI)*8.232E+5-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*5.0E+9-x*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9+x*y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI*2.0)*sin(x*M_PI)*4.2E+7+y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*M_PI*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*4.2E+7+M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*3.92E+8-x*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(y*M_PI)*1.176E+7-y*M_PI*cos(t*M_PI)*cos(t*M_PI*2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(t*M_PI*2.0)*sin(x*M_PI)*1.176E+7-x*y*M_PI*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*pow(cos(x*M_PI),2.0)*cos(y*M_PI)*sin(y*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*M_PI*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*pow(cos(y*M_PI),2.0)*sin(t*M_PI)*sin(x*M_PI)*sqrt(pow(y+cos(t*M_PI)*(7.0/2.5E+1)-1.0/2.0,2.0)+pow(-x+sin(t*M_PI)*(7.0/2.5E+1)+1.0/2.0,2.0))*1.4E+9-x*y*(M_PI*M_PI)*cos(t*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7-x*y*(M_PI*M_PI)*pow(cos(t*M_PI*2.0),2.0)*cos(x*M_PI)*cos(y*M_PI)*sin(t*M_PI)*sin(x*M_PI)*sin(y*M_PI)*1.05E+7))/(x*-1.5625E+7-y*1.5625E+7-cos(t*M_PI)*4.375E+6+sin(t*M_PI)*4.375E+6+y*cos(t*M_PI)*8.75E+6-x*sin(t*M_PI)*8.75E+6+(x*x)*1.5625E+7+(y*y)*1.5625E+7+pow(cos(t*M_PI),2.0)*1.225E+6+pow(sin(t*M_PI),2.0)*1.225E+6+7.8125E+6);;
}

R fun_one(double *P, const int i) { return 1.; }

} // namespace Example1Omega2


// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// Note: standard for this program is to solve the equations on the inner domain
// Omega 2. To solve on Omega 1, define the word "omega1" for the pre-processor
// (only works for example 1)

//* Set numerical example (options: "example1", "example2")
#define example1

//* Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominated

//* Choose domain to solve on (options: "omega1", "omega2")
#define omega2
#define neumann // boundary condition on outer domain (options: "dirichlet", "neumann")

//* Set scheme for the method (options: "classical", "conservative".)
#define classical

//* Set stabilization method (options: "fullstab", "macro")
#define fullstab

//* Decide whether to solve for level set function, or to use exact (options:
// "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h    // to set mesh size using the h parameter. Write use_n to decide
                 // using nx, ny.
#define use_tnot // write use_t to control dT manually. Otherwise it is set
                 // proportional to h.

// Do not touch
#ifdef example1
#ifdef omega1
//using namespace Example1_new;
using namespace Kex;
#else
//using namespace Example1_omega2; // on Omega 2
using namespace Example1Omega2;
#endif
#endif
#if defined(lehrenfeld)
using namespace Lehrenfeld_Convection_Dominated;
#elif defined(example2)
using namespace Example2;
#endif

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);
    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;            // starting mesh size
    double dT = 0.125;

    int total_number_iteration;
    double time_step;
    const double t0 = 0.;

    int thread_count = 1;

#ifdef example1
    // Paths to store data
    const std::string path_output_data    = "../output_files/coupled_langmuir/example1/data/";
    const std::string path_output_figures = "../output_files/coupled_langmuir/example1/paraview/";
#elif defined(example2)
    const std::string path_output_data    = "../output_files/coupled_langmuir/example2/data/";
    const std::string path_output_figures = "../output_files/coupled_langmuir/example2/paraview/";
#endif

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    // Data file to hold problem data
    std::ofstream output_data(path_output_data + "data.dat", std::ofstream::out);
    // Data file to hold DOF indices
    std::ofstream indices(path_output_data + "indices.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors_bulk;    // array to hold bulk errors
    std::array<double, iterations> errors_surface; // array to hold surface errors
    std::array<double, iterations> hs;             // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
#ifdef example1
        const double x0 = 0., y0 = 0., lx = 1., ly = 1.;
#elif defined(example2)
        const double x0 = -1., y0 = -1., lx = 2., ly = 2.;
#endif

#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h                          = lx / (nx - 1);
#endif

        const Mesh Th(nx, ny, x0, y0, lx, ly);

        // Parameters
        const double tfinal = .1; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 4;

        double dT              = h / divisionMeshSize;
        // double dT = 3*h;
        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "Iteration " << j + 1 << "/" << iterations << '\n';
        }

        std::cout << "h  = " << h << '\n';
        std::cout << "nx = " << nx << '\n';
        std::cout << "ny = " << ny << '\n';
        std::cout << "dT = " << dT << '\n';

        const double D      = 0.01; // Diffusion parameter bulk
        const double DGamma = 1.;   // Diffusion parameter surface

        const double lambda = 1.; // Nitsche's method penalty parameter for Dirichlet BCs

        // Stabilization parameters
        const double tau1 = 1., tau2 = 1.;
        double stab_bulk_face = h * tau1;
        double stab_surf_face = tau2;

        const FESpace2 Vh(Th, DataFE<Mesh>::P2); // continuous basis functions

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        const FESpace2 Vh2(Th,
                           DataFE<Mesh>::P2); // higher order space for interpolation
        const FESpace2 Vh3(Th,
                           DataFE<Mesh>::P3); // higher order space for interpolation

        // 1D Time mesh
        const double final_time = total_number_iteration * time_step;
        const Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        const FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly);
        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(3)); // specify order of quadrature in time
        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Matrix to hold non-linear part in Newton iteration
        std::vector<std::map<std::pair<int, int>, double>> jacobian(thread_count);

        // Velocity field
        const Lagrange2 FEvelocity(1);

        const FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        const FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        const double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

#if defined(levelsetexact) && not defined(example2)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);
#elif defined(levelsetsolve) || defined(example2)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet);
#endif

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double q0_0, q0_0_surface, q0_1, q0_1_surface, qp_0, qp_0_surface, qp_1,
            qp_1_surface;                                             // integral values to be computed
        double intF = 0, int_outflow = 0, intF_surface = 0, intG = 0; // hold integrals of rhs and Neumann bcs
        // double int_uB = 0, int_uB_gamma = 0, int_uS = 0, int_uB_init = 0, int_uS_init = 0;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * time_step;
            const TimeSlab &In(Ih[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';
            std::cout << "dT = " << dT << '\n';

            swap(ls[0], ls[lastQuadTime]);

            // computation of the interface in each of the three quadrature points
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
            ActiveMesh<Mesh> Thi(Th);

            //Thi.truncate(interface, -1); // remove part with negative sign of level set
            Thi.truncate(interface, 1); // remove part with negative sign of level set

            // Cut FE space
            CutSpace Wh(Thi, Vh);

            // Surface active mesh
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            // Surface FE space
            CutFESpace WhGamma(ThGamma, Vh);

            // Add time slab to cut space
            convdiff.initSpace(Wh, In);

            // Add time slab to surface space on top of the data
            convdiff.add(WhGamma, In);

            // Objects needed for the weak form
            const Normal n;
            const Tangent t;

            // Right hand side functions
            Fun_h f(Vh, In, fun_rhsBulk);
            Fun_h fS(Vh2, In, fun_rhsSurf);

#ifndef example2
            Fun_h g(Vh, In, fun_uBulk);                         // create an FE-function of the exact bulk solution
            // Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);         // computed Neumann BC
            // Fun_h g_Neumann_left(Wh, In, fun_neumann_left);     // label 1
            // Fun_h g_Neumann_bottom(Wh, In, fun_neumann_bottom); // label 4
            // Fun_h g_Neumann_right(Wh, In, fun_neumann_right);   // label 2
            // Fun_h g_Neumann_top(Wh, In, fun_neumann_top);       // label 3
#endif

            // Test and Trial functions
            FunTest uB(Wh, 1), vB(Wh, 1);           // bulk
            FunTest uS(WhGamma, 1), vS(WhGamma, 1); // bulk

            // Data for initial solution

            // solution array looks like [uB^0 uB^1 uS^0 uS^1]
            // so to access the DOFs corresponding to the first point
            // in the time slab, we want uB^0 for the bulk and uS^0
            // for the surface
            int idx_s0  = Wh.NbDoF() * In.NbDoF();      // index where uS^0 start in the solution array
            int idx_s1  = WhGamma.NbDoF() * In.NbDoF(); // degrees of freedom in surface x time space
            int tot_dof = (int)convdiff.get_nb_dof();   // total degrees of freedom

            // Data for initial solution
            Rn data_init(tot_dof, 0.);                                    // initial data total     
            convdiff.initialSolution(data_init);
            Rn data_all(data_init);

            KN_<double> data_uhB0(data_init(SubArray(Wh.NbDoF(), 0)));                      // data for uhB(t_{n-1}^-)
            KN_<double> data_uhB(data_all(SubArray(Wh.NbDoF() * In.NbDoF(), 0)));           // data for uhB
            KN_<double> data_uhS0(data_init(SubArray(WhGamma.NbDoF(), idx_s0)));            // data for uhS(t_{n-1}^-)
            KN_<double> data_uhS(data_all(SubArray(WhGamma.NbDoF() * In.NbDoF(), idx_s0))); // data for uhS

            FunFEM<Mesh2> uhB0(Wh, data_uhB0);
            FunFEM<Mesh2> uhB(Wh, data_uhB);
            FunFEM<Mesh2> uhS0(WhGamma, data_uhS0);
            FunFEM<Mesh2> uhS(WhGamma, data_uhS);

            if (iter == 0) {
                interpolate(Wh, data_uhB0, fun_uBulkInit);
                interpolate(WhGamma, data_uhS0, fun_uSurfInit);
            }

            int newton_iterations = 0;
            double tol            = 1e-10;
            double residual       = 1.;

            // Newton's method
            while (1) {
                std::cout << " Newton iteration \t : \t" << newton_iterations << '\n';

                // Test functions in Jacobian (DF(u)(w))
                FunTest wB(Wh, 1), wS(WhGamma, 1);

                // Assemble linear part of Jacobian DF
                if (newton_iterations == 0) {

#ifdef conservative

                     // Time derivative terms
                    convdiff.addBilinear(-innerProduct(wB, dt(vB)), Thi, In);
                    convdiff.addBilinear(-innerProduct(wS, dt(vS)), interface, In);

                    // Time penalty terms (only) LHS
                    convdiff.addBilinear(+innerProduct(wB, vB), Thi, lastQuadTime, In);
                    convdiff.addBilinear(+innerProduct(wS, vS), *interface(lastQuadTime), In, lastQuadTime);

                    // Convection
                    convdiff.addBilinear(-innerProduct(wB, (vel.exprList() * grad(vB))), Thi, In);
                    convdiff.addBilinear(-innerProduct(wS, (vel.exprList() * grad(vS))),
                                         interface, In);

                    // Outer boundary term 
                    convdiff.addBilinear(innerProduct(vel.exprList() * wB * n, vB), Thi, INTEGRAL_BOUNDARY, In);

#elif defined(classical)

                    // Time derivative terms
                    convdiff.addBilinear(+innerProduct(dt(wB), vB), Thi, In);
                    convdiff.addBilinear(+innerProduct(dt(wS), vS), interface, In);

                    // Time penalty terms (only) LHS
                    convdiff.addBilinear(+innerProduct(wB, vB), Thi, 0, In);
                    convdiff.addBilinear(+innerProduct(wS, vS), *interface(0), In, 0);

                    // Convection
                    convdiff.addBilinear(+innerProduct((vel.exprList() * grad(wB)), vB), Thi, In);
                    convdiff.addBilinear(+innerProduct((vel.exprList() * grad(wS)), vS) +innerProduct(wS * divS(vel), vS),
                                         interface, In);
#endif
                    // Diffusion
                    convdiff.addBilinear(+innerProduct(D * grad(wB), grad(vB)), Thi, In);
                    convdiff.addBilinear(+innerProduct(DGamma * gradS(wS), gradS(vS)), interface, In);

                    

                    // Linear part of coupling term
                    convdiff.addBilinear(innerProduct(wB - wS, vB - vS), interface, In);

                    // Stabilization
                    // convdiff.addFaceStabilization(
                    //     +innerProduct(stab_bulk_face * jump(grad(wB) * n), jump(grad(vB) * n)), Thi, In);
                    convdiff.addFaceStabilization(
                        + innerProduct(h * tau1 * jump(grad(wB) * n), jump(grad(vB) * n)) 
                        + innerProduct(h * h * h * tau1 * jump(grad(grad(wB) * n) * n), jump(grad(grad(vB) * n) * n)),
                    Thi, In);

                    // convdiff.addFaceStabilization(
                    //     +innerProduct(stab_surf_face * jump(grad(wS) * n), jump(grad(vS) * n)), ThGamma, In);
                    convdiff.addFaceStabilization(
                        + innerProduct(tau2 * jump(grad(wS) * n), jump(grad(vS) * n)) 
                        + innerProduct(h * h * tau2 * jump(grad(grad(wS) * n) * n), jump(grad(grad(vS) * n) * n)),
                        ThGamma, In);

                    convdiff.addBilinear(
                        + innerProduct(tau2 * grad(wS) * n, grad(vS) * n)
                        + innerProduct(tau2 * h * h * grad(grad(wS) * n) * n, grad(grad(vS) * n) * n)
                        , interface
                        , In);

                }

                // Multiply the linear part of the Jacobian by the initial guess to evaluate in the initial guess (which
                // is updated in each Newton iteration)
                convdiff.addMatMul(data_all);   // -> goes into the rhs (this does not change the matrix convdiff.mat_)

                // Save the Jacobian matrix
                jacobian[0] = convdiff.mat_[0];

                // Switch the matrix to the jacobian instead of convdiff.mat_
                convdiff.set_map(jacobian);

                // Assemble remaining non-linear part of Jacobian (evaluated in the initial guess so it becomes linear)
                convdiff.addBilinear(
                    - innerProduct(wB * uhS.expr(), vB - vS) 
                    - innerProduct(wS * uhB.expr(), vB - vS)
                    , interface
                    , In);

                
                // Add the terms of the residual remaining from the Jacobian (construct F) as rhs terms (addLinear)
                convdiff.addLinear(-innerProduct(uhB.expr() * uhS.expr(), vB - vS), interface, In);

                // Time penalty RHS
                convdiff.addLinear(-innerProduct(uhB0.expr(), vB), Thi, 0, In);
                convdiff.addLinear(-innerProduct(uhS0.expr(), vS), *interface(0), In, 0);

                // RHS force functions
                convdiff.addLinear(-innerProduct(f.expr(), vB), Thi, In);        // rhs in bulk
                convdiff.addLinear(-innerProduct(fS.expr(), vS), interface, In); // rhs on surface
                
                // Solve DF(u0)(w) = F(u0)
                convdiff.solve(jacobian[0], convdiff.rhs_);

                // Retrieve the bulk solution
                KN_<double> data_wB(convdiff.rhs_(SubArray(Wh.NbDoF() * In.NbDoF(), 0)));
                // Retrieve the surface solution
                KN_<double> data_wS(convdiff.rhs_(SubArray(WhGamma.NbDoF() * In.NbDoF(), idx_s0)));

                // Compute the residual norms
                double residual_norm_B = data_wB.l1() / data_wB.size();
                double residual_norm_S = data_wS.l1() / data_wS.size();
                residual               = std::max(residual_norm_B, residual_norm_S);

                std::cout << " Residual norm bulk \t : \t" << residual_norm_B << '\n';
                std::cout << " Residual norm surface \t : \t" << residual_norm_S << '\n';
                std::cout << " Residual \t : \t" << residual << '\n';

                // Get the total residual
                KN_<double> data_w(convdiff.rhs_(SubArray(convdiff.get_nb_dof())));
                // Update initial guess
                data_all -= data_w;

                // Reset the rhs vector
                convdiff.rhs_.resize(convdiff.get_nb_dof());
                convdiff.rhs_ = 0.0; // reset rhs

                // Update Newton iteration counter
                newton_iterations++;

                if (newton_iterations > 5 || residual < tol) {
                    convdiff.saveSolution(data_all);
                    convdiff.cleanBuildInMatrix();
                    convdiff.set_map();     // return matrix to mat_ instead of jacobian
                    break;
                }
            }

            // Compute conservation error
            if (iterations == 1) {

                // auto outflow = (vel * n) * funuh.expr();
                //  int_outflow =
                //      integral(outflow, In, interface, 0); // integral(vel.exprList() * funuh * n, In, interface, 0);
                //  int_uB = integral(Thi, In, b0h, 0, qTime);
                //  int_uB_gamma       = integral(b0h.expr(), In, interface, 0);
                //  int_uS       = integral(s0h.expr(), In, interface, 0);

                intF         = integral(Thi, In, f, 0, qTime);
                intF_surface = integral(fS, In, interface, 0);

                Fun_h funuh(Wh, data_all);

                //* bulk

                // std::cout << "Wh.NbDoF() = "<< Wh.NbDoF() << ", idx_s0 = " << idx_s0 << ", tot_dof = " << tot_dof <<
                // ", WhGamma.get_nb_dof() = " << WhGamma.get_nb_dof() << std::endl; getchar();

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_all(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Thi, funsol, 0, 0); // integral of solution in bulk in t_{n-1}
                sol2 += data_all(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Thi, funsol, 0, lastQuadTime); // integral of solution in bulk in t_{n}

                //* surface

                Rn solsurf(WhGamma.NbDoF(), 0.);
                Fun_h funsolsurf(WhGamma, solsurf);
                solsurf += data_all(SubArray(WhGamma.NbDoF(), idx_s0));
                double q_0_surf = integral(funsolsurf, interface(0), 0);
                solsurf += data_all(SubArray(WhGamma.NbDoF(), idx_s0 + WhGamma.NbDoF()));
                double q_1_surf = integral(funsolsurf, interface(lastQuadTime), 0);

                Fun_h funuh_surf(WhGamma, solsurf);
                if (iter == 0) {
                    q0_0 = q_0;
                    q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Thi, uhB0, 0, 0);

                    q0_0_surface = q_0_surf;
                    q0_1_surface = q_1_surf;
                    qp_1_surface = q_1_surf;
                    q0_1_surface = integral(uhS0, interface(0), 0);

                    // std::cout << "total mass: " << (q0_1 + q0_1_surface) << "\n";
                    // return 0;
                }

                output_data << std::setprecision(10);
                output_data << current_time << "\t\t,"

#ifndef example2
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface)) << "\t\t,"
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface) - intF - intF_surface) << "\n";
#else
                            << (q_1 + q_1_surf) << "\t\t," << (q0_1 + q0_1_surface) << "\t\t,"
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface)) << "\n";
#endif
                

                qp_1_surface = q_1_surf;
                qp_1         = q_1;
            }

#ifndef example2
            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_all(SubArray(Wh.get_nb_dof(), 0));               // add first time dof
            sol += data_all(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof())); // add second time dof
            Fun_h funuh(Wh, sol);
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            Rn solsurf(WhGamma.get_nb_dof(), 0.);
            solsurf += data_all(SubArray(WhGamma.get_nb_dof(), idx_s0));
            solsurf += data_all(SubArray(WhGamma.get_nb_dof(), idx_s0 + WhGamma.get_nb_dof()));
            Fun_h funuh_surf(WhGamma, solsurf);
            double errSurf = L2normSurf(funuh_surf, fun_uSurf, *interface(lastQuadTime), current_time + dT, 0, 1);
            std::cout << " t_n -> || uS-uSex||_2 = " << errSurf << '\n';

            errors_bulk.at(j)    = errBulk;
            errors_surface.at(j) = errSurf;

            if (iterations > 1 && iter == total_number_iteration - 1)
                output_data << h << "," << dT << "," << errBulk << '\n';

#endif

            // #ifdef USE_MPI
            if ((iterations == 1) && MPIcf::IamMaster()) {
                // #else
                // Fun_h sol(Wh, data_all);
                Fun_h sol(Wh, data_all);
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk" + std::to_string(iter + 1) + ".vtk");
                writer.add(uhB0, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                // Expression2 uuh(sol, 0, op_id);
                // Expression2 uuex(uBex, 0, op_id);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(fabs(uhB0.expr() - uBex.expr()), "bulk_error");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);

                Paraview<Mesh> writer_surface(ThGamma,
                                              path_output_figures + "surface" + std::to_string(iter + 1) + ".vtk");
                writer_surface.add(uhS0, "surface", 0, 1);
                Fun_h uSex(WhGamma, fun_uSurf, current_time);
                Fun_h fS(WhGamma, fun_rhsSurf, current_time);
                // Expression2 uuhsurf(sol, 0, op_id);
                // Expression2 uuex(uBex, 0, op_id);
                writer_surface.add(fabs(uhS0.expr() - uSex.expr()), "surface_error");
                writer_surface.add(uSex, "surface_exact", 0, 1);
                writer_surface.add(fS, "surface_rhs", 0, 1);
                writer_surface.add(ls[0], "levelSet0", 0, 1);
                writer_surface.add(ls[1], "levelSet1", 0, 1);
                writer_surface.add(ls[2], "levelSet2", 0, 1);
            }

            iter++;
        }

// Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_t)
        dT *= 0.5;
#elif defined(use_h)
        h *= 0.5;
#endif
    }

    std::cout << '\n';
    std::cout << "Errors Bulk = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_bulk.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "Errors Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_surface.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << "dT = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << dts.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}
