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
    R fun_levelSet(double *P) {
        return +(sqrt((P[0]-0.1)*(P[0]-0.1) + (P[1]-0.0)*(P[1]-0.0)) - 0.3) + Epsilon;
    }

    // R fun_levelSet(double *P, const int i, const R t) {
    //     return +(sqrt((P[0]-0.1)*(P[0]-0.1) + (P[1]-0.0)*(P[1]-0.0)) - 0.3) + Epsilon;
    // }
    
    R fun_rhsBulk(double *P, const int i, const R t) {return 0.;}
    R fun_rhsSurf(double *P,const int i, const R t) {return 0.;}
    R fun_neumann_Gamma(double *P, const int i) {return 0.;}
    R fun_velocity(double *P, const int i){
        if(i == 0) return -0.5*(1+cos(M_PI*P[0]))*sin(M_PI*P[1]);
        else       return  0.5*(1+cos(M_PI*P[1]))*sin(M_PI*P[0]);
    }
//    R fun_velocity(double *P, const int i){
//        return 0.0;
//    }
    R fun_uSurfInit(double *P, int i) { return 0.;}
    R fun_w(double r) {
        return 0.5*(1-cos((r-0.3)*M_PI/(0.5*0.3)));
    }
    R fun_uBulkInit(double *P, int i) {
        double r = sqrt((P[0]-0.1)*(P[0]-0.1) + P[1]*P[1]);
        if(r > 1.5*0.3) return 0.5*(1-P[0]*P[0])*(1-P[0]*P[0]);
        else if(0.3 < r && r <= 1.5*0.3) return 0.5*(1-P[0]*P[0])*(1-P[0]*P[0])*fun_w(r);
        else return 0.;}

    // below are just so that program does not crash, it is not known
    R fun_uBulkD(double *P, const int i, const int d, const R t) {
        return 0.;
    }

    R fun_uBulk(double *P, const int i, const R t) {
        return 0.;
    }

    R fun_uSurf(double *P, const int i, const R t) {
        return 0.;
    }
}

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
#define example2

//* Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominated

//* Choose domain to solve on (options: "omega1", "omega2")
#define omega1
#define neumann // boundary condition on outer domain (options: "dirichlet", "neumann")

//* Set scheme for the method (options: "classical", "conservative".)
#define conservative

//* Set stabilization method (options: "fullstab", "macro")
#define fullstab

//* Decide whether to solve for level set function, or to use exact (options:
// "levelsetsolve", "levelsetexact")
#define levelsetsolve

#define use_h       // to set mesh size using the h parameter. Write use_n to decide
                    // using nx, ny.
#define use_tnot    // write use_t to control dT manually. Otherwise it is set
                    // proportional to h.

// Do not touch
#ifdef example1
#ifdef omega1
using namespace Example1_new;
#else
using namespace Example1_omega2; // on Omega 2
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
    const size_t iterations = 1;    // number of mesh refinements   (set to 1 to run
                                    // only once and plot to paraview)
    int nx = 15, ny = 15;           // starting mesh size
    double h  = 0.05;              // starting mesh size
    double dT = 0.125;

    int total_number_iteration;
    double time_step;
    const double t0 = 0.;

#ifdef example1
    // Paths to store data
    const std::string path_output_data    = "../output_files/coupled/example1/data/";
    const std::string path_output_figures = "../output_files/coupled/example1/paraview/";
#elif defined(example2)
    const std::string path_output_data    = "../output_files/coupled/example2/data/";
    const std::string path_output_figures = "../output_files/coupled/example2/paraview/";
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
        h = lx / (nx - 1);
#endif

        const Mesh Th(nx, ny, x0, y0, lx, ly);

        // Parameters
        const double tfinal = 2.0; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 3;

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

        const FESpace2 Vh(Th, DataFE<Mesh>::P1); // continuous basis functions

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        const FESpace2 Vh2(Th,
                           DataFE<Mesh>::P2); // higher order space for interpolation

        // 1D Time mesh
        const double final_time = total_number_iteration * time_step;
        const Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        const FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(3)); // specify order of quadrature in time
        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

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
        //double int_uB = 0, int_uB_gamma = 0, int_uS = 0, int_uB_init = 0, int_uS_init = 0;
        
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

            Thi.truncate(interface, -1); // remove part with negative sign of level set

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
            Fun_h g(Vh, In, fun_uBulk); // create an FE-function of the exact bulk solution
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);         // computed Neumann BC
            Fun_h g_Neumann_left(Wh, In, fun_neumann_left);     // label 1
            Fun_h g_Neumann_bottom(Wh, In, fun_neumann_bottom); // label 4
            Fun_h g_Neumann_right(Wh, In, fun_neumann_right);   // label 2
            Fun_h g_Neumann_top(Wh, In, fun_neumann_top);       // label 3
        #endif

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);             // bulk
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
            Rn data_u0(convdiff.get_nb_dof(),0.);                                    // initial data total
            //convdiff.initialSolution(data_u0);                          // allocate memory        //! DANGEROUS: THIS RUNS BUT GIVES INCORRECT SOLUTION
            convdiff.initialSolution(data_u0);                          // allocate memory

            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));           // initial data bulk
            KN_<R> data_S0(data_u0(SubArray(WhGamma.NbDoF(), idx_s0))); // initial data surface

            if (iter == 0) {
                // Interpolate initial solution
                interpolate(Wh, data_B0, fun_uBulkInit);
                interpolate(WhGamma, data_S0, fun_uSurfInit);

                Fun_h uB_init(Vh, fun_uBulkInit);
                Fun_h uS_init(Vh, fun_uSurfInit);
            }

            // Make function objects to use in bilinear forms
            Fun_h b0h(Wh, data_B0);
            Fun_h s0h(WhGamma, data_S0);

            //** Assembling linear and bilinear forms

            //* Time terms

#ifdef conservative

            // Bulk
            convdiff.addBilinear(-innerProduct(u, dt(v)), Thi, In);
            // Surface
            convdiff.addBilinear(-innerProduct(uS, dt(vS)), interface, In);

            // Time penalty term bulk LHS
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);

            // Time penalty term surface LHS
            convdiff.addBilinear(+innerProduct(uS, vS), *interface(lastQuadTime), In, (int)lastQuadTime);

            // Time penalty term surface RHS

            // if (iter == 0) convdiff.addLinear(fun_uSurfInitT, +innerProduct(1., vS), *interface(0), In, 0);
            // else
            convdiff.addLinear(+innerProduct(s0h.expr(), vS), *interface(0), In, 0);

#else // classical scheme

            // Bulk
            convdiff.addBilinear(+innerProduct(dt(u), v), Thi, In);
            // Surface
            convdiff.addBilinear(+innerProduct(dt(uS), vS), interface, In);

            // Time penalty term bulk LHS
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);

            // Time penalty term surface LHS
            convdiff.addBilinear(+innerProduct(uS, vS), *interface(0), In, 0);

            // Time penalty term surface RHS
            convdiff.addLinear(+innerProduct(s0h.expr(), vS), *interface(0), In, 0);

#endif

            //* Diffusion terms

            // Bulk
            convdiff.addBilinear(+innerProduct(D * grad(u), grad(v)), Thi, In);

            // Surface
            convdiff.addBilinear(+innerProduct(DGamma * gradS(uS), gradS(vS)), interface, In);

            //* Convection terms

#if defined(classical)
            // Bulk
            convdiff.addBilinear(+innerProduct((vel.exprList() * grad(u)), v), Thi, In);

            // Surface
            convdiff.addBilinear(+innerProduct((vel.exprList() * grad(uS)), vS) + innerProduct(uS * divS(vel), vS),
                                 interface, In);

#elif defined(conservative)
            convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Thi, In);
            convdiff.addBilinear(-innerProduct(uS, (vel.exprList() * grad(vS))), interface, In);
#endif

            //* Stabilization
            // Stabilization parameters
            const double tau1 = 1., tau2 = 1., tau3 = 0.;

            double stab_bulk_face      = h * tau1;
            double stab_surf_face      = tau2;
            double stab_surf_interface = 0.; //h *tau3;

            // Stabilization along the interface
            // convdiff.addBilinear(+innerProduct(stab_surf_interface * grad(uS) * n, grad(vS) * n), interface, In);

#if defined(fullstab)

            convdiff.addFaceStabilization(+innerProduct(stab_bulk_face * jump(grad(u) * n), jump(grad(v) * n)), Thi,
                                          In);
            convdiff.addFaceStabilization(+innerProduct(stab_surf_face * jump(grad(uS) * n), jump(grad(vS) * n)),
                                          ThGamma, In);

#elif defined(macro)

            TimeMacroElement<Mesh> TimeMacro(Thi, qTime, 0.25);
            TimeMacroElementSurface<Mesh> TimeMacroGamma(ThGamma, interface, qTime, 0.25);

            // Visualize macro elements
            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, path_output_figures + "Th" + std::to_string(iter + 1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

                // domain = 0,

                writerMacro.writeFaceStab(
                    Thi, 0, path_output_figures + "FullStabilization" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeActiveMesh(Thi,
                                            path_output_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0,
                                              path_output_figures + "macro" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroInnerEdge(
                    TimeMacro, 0, path_output_figures + "macro_inner_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroOutterEdge(
                    TimeMacro, 0, path_output_figures + "macro_outer_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeSmallElements(
                    TimeMacro, 0, path_output_figures + "small_element" + std::to_string(iter + 1) + ".vtk");

                Paraview<Mesh> writerMacroSurf(ThGamma,
                                               path_output_figures + "ThGamma" + std::to_string(iter + 1) + ".vtk");
                writerMacroSurf.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacroSurf.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacroSurf.add(ls[2], "levelSet2.vtk", 0, 1);

                writerMacroSurf.writeFaceStab(
                    ThGamma, 0, path_output_figures + "FullStabilizationSurface" + std::to_string(iter + 1) + ".vtk");
                writerMacroSurf.writeActiveMesh(ThGamma, path_output_figures + "ActiveMeshSurface" +
                                                             std::to_string(iter + 1) + ".vtk");
                writerMacroSurf.writeMacroInnerEdge(TimeMacroGamma, 0,
                                                    path_output_figures + "macro_inner_edge_surface" +
                                                        std::to_string(iter + 1) + ".vtk");
                writerMacroSurf.writeMacroOutterEdge(TimeMacroGamma, 0,
                                                     path_output_figures + "macro_outer_edge_surface" +
                                                         std::to_string(iter + 1) + ".vtk");
                writerMacroSurf.writeSmallElements(TimeMacroGamma, 0,
                                                   path_output_figures + "small_element_surface" +
                                                       std::to_string(iter + 1) + ".vtk");
            }

            // Stabilization of the bulk
            // convdiff.mat_.clear();
            convdiff.addFaceStabilization(+innerProduct(stab_bulk_face * jump(grad(u) * n), jump(grad(v) * n)), Thi, In,
                                          TimeMacro);
            convdiff.addFaceStabilization(+innerProduct(stab_surf_face * jump(grad(uS) * n), jump(grad(vS) * n)),
                                          ThGamma, In, TimeMacroGamma);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // getchar();

#endif

            //* Boundary conditions

            // Coupling term
            // convdiff.addBilinear(innerProduct(jump(u,uS), jump(v,vS)), interface, In);
            convdiff.addBilinear(innerProduct(u - uS, v - vS), interface, In);

#if defined(omega1)
#if defined(dirichlet)
            // Dirichlet outer
            convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v), // added penalty
                                 Thi, INTEGRAL_BOUNDARY, In);

            convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v), Thi,
                               INTEGRAL_BOUNDARY, In);

#if defined(conservative)
            // Set inflow and outflow conditions

            // convdiff.addLinear(-innerProduct(g.expr(), vel*n*v), Thi, INTEGRAL_BOUNDARY, In);
            convdiff.addLinear(+innerProduct(g.expr(), 0.5 * fabs(vel * n) * v) -
                                   innerProduct(g.expr(), 0.5 * (vel * n) * v),
                               Thi, INTEGRAL_BOUNDARY, In);

            convdiff.addBilinear(+innerProduct(u, 0.5 * fabs(vel * n) * v) + innerProduct(u, 0.5 * (vel * n) * v), Thi,
                                 INTEGRAL_BOUNDARY, In);
#endif

#elif defined(neumann)
            // add nothing from Neumann BC, since it's zero for this example

#if defined(conservative)
            convdiff.addBilinear(innerProduct(vel.exprList() * u * n, v), Thi, INTEGRAL_BOUNDARY, In);
#endif

#endif
#endif

            convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);         // rhs in bulk
            convdiff.addLinear(+innerProduct(fS.expr(), vS), interface, In); // rhs on surface
            // convdiff.addLinear(fun_rhsSurf, +innerProduct(1., vS), interface, In);

            // Write matrix and indices of DOFS to file
            if (iter == total_number_iteration - 1) {

                matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");

                indices << idx_s0 << ", " << idx_s1
                        << "\n"; // export indices needed for computing the scaled condition number
            }

            // Solve linear system
            convdiff.solve("mumps");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

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

                Fun_h funuh(Wh, data_u0);

                //* bulk

                // std::cout << "Wh.NbDoF() = "<< Wh.NbDoF() << ", idx_s0 = " << idx_s0 << ", tot_dof = " << tot_dof <<
                // ", WhGamma.get_nb_dof() = " << WhGamma.get_nb_dof() << std::endl; getchar();

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_u0(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Thi, funsol, 0, 0); // integral of solution in bulk in t_{n-1}
                sol2 += data_u0(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Thi, funsol, 0, lastQuadTime); // integral of solution in bulk in t_{n}

                //* surface

                Rn solsurf(WhGamma.get_nb_dof(), 0.);
                Fun_h funsolsurf(WhGamma, solsurf);
                solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0));
                double q_0_surf = integral(funsolsurf, interface(0), 0);
                solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0 + WhGamma.get_nb_dof()));
                double q_1_surf = integral(funsolsurf, interface(lastQuadTime), 0);

                Fun_h funuh_surf(WhGamma, solsurf);
                if (iter == 0) {
                    q0_0 = q_0;
                    q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Thi, b0h, 0, 0);

                    q0_0_surface = q_0_surf;
                    q0_1_surface = q_1_surf;
                    qp_1_surface = q_1_surf;
                    q0_1_surface = integral(s0h, interface(0), 0);

                    // std::cout << "total mass: " << (q0_1 + q0_1_surface) << "\n";
                    // return 0;
                    
                }

                                output_data << std::setprecision(10);
                output_data << current_time
                            << "\t\t,"
                            //<< ((q_1 - qp_1) - intF - int_outflow - intG) << "\t\t," << ((q_1_surf - qp_1_surface) -
                            // intF_surface) << "\t,"
                            //<< ((q_1 + q_1_surf) - (q0_0 + q0_0_surface) - intF - intF_surface) << '\t\t'
                        #ifndef example2
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface)) << "\t\t,"
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface) - intF - intF_surface) << "\n";
                        #else
                            << (q_1 + q_1_surf) << "\t\t,"
                            << (q0_1 + q0_1_surface) << "\t\t,"
                            << ((q_1 + q_1_surf) - (qp_1 + qp_1_surface)) << "\n";
                        #endif
                //<< ((q_1 - qp_1) - intF + int_uB_gamma - int_uS) << "\t\t,"
                //<< ((q_1 + q_1_surf) - (int_uB_init + int_uS_init) - intF - intF_surface) << "\t\t,"
                //<< (q_1 + q_1_surf) << '\n';

                qp_1_surface = q_1_surf;
                qp_1         = q_1;
            }

        #ifndef example2
            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));               // add first time dof
            sol += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof())); // add second time dof
            Fun_h funuh(Wh, sol);
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            Rn solsurf(WhGamma.get_nb_dof(), 0.);
            solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0));
            solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0 + WhGamma.get_nb_dof()));
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
                Fun_h sol(Wh, data_u0);

            #ifdef conservative
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_conservative" + std::to_string(iter + 1) + ".vtk");
            #else
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_nonconservative" + std::to_string(iter + 1) + ".vtk");
            #endif
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                // Expression2 uuh(sol, 0, op_id);
                // Expression2 uuex(uBex, 0, op_id);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);

            #ifdef conservative
                Paraview<Mesh> writer_surface(ThGamma,
                                              path_output_figures + "surface_conservative" + std::to_string(iter + 1) + ".vtk");
            #else
                Paraview<Mesh> writer_surface(ThGamma,
                                              path_output_figures + "surface_nonconservative" + std::to_string(iter + 1) + ".vtk");
            #endif
                writer_surface.add(s0h, "surface", 0, 1);
                Fun_h uSex(WhGamma, fun_uSurf, current_time);
                Fun_h fS(WhGamma, fun_rhsSurf, current_time);
                // Expression2 uuhsurf(sol, 0, op_id);
                // Expression2 uuex(uBex, 0, op_id);
                writer_surface.add(fabs(s0h.expr() - uSex.expr()), "surface_error");
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
