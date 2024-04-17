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
 * @note We consider a time-dependent bulk problem in a fictitious domain called Omega(t).

 *  Problem:

    * Omega_2(t):
    Find u in Omega_2(t) such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    in Omega_2(t).
                                        BC,  on Gamma(t).

 *  Numerical method:
    A space-time Cutfem, using the level-set method.

 *  Classical scheme: Standard FEM scheme.
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
#include "../num/gnuplot.hpp"
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "projection.hpp"
#include "../num/matlab.hpp"
#include "../num/redirectOutput.hpp"
#include "paraview.hpp"
#include "../problem/AlgoimIntegration.hpp"

using namespace globalVariable; // to access some globally defined constants

// Numerical examples

namespace Example1 {

const double D      = 0.01;
const double DGamma = 1.;

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17);
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
        return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * (x(0) - 0.5 - 0.28 * sin(M_PI * t)),
                                     2.0 * (x(1) - 0.5 + 0.28 * cos(M_PI * t)));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(pow(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))), 2) +
                      pow(2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))), 2));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))) / norm,
                  2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))) / norm);
    }
};

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // In Omega_2(t), there is no Neumann boundary condition, since the only boundary is Gamma(t)
    // to which only the coupling conditions apply
    return 0;
}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
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

R fun_one(double *P, const int i) { return 1.; }

} // namespace Example1

//* Set numerical example (options: "ex1", "ex2", or "ex3")
#define ex1
//* Set scheme for the method (options: "classical", "conservative")
#define conservative
//* Set stabilization method (options: "fullstab", "macro")
#define macro

#define use_h    // to set mesh size using the h parameter. Write use_n to decide
                 // using nx, ny.
#define use_tnot // write use_t to control dT manually. Otherwise it is set
                 // proportional to h.
#define conservation

// Setup two-dimensional class types
const int d = 2;

typedef MeshQuad2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

// Do not touch
#if defined(ex1)
using namespace Example1;
#elif defined(ex2)
// using namespace Example2;
using namespace Kite;
#elif defined(ex3)
using namespace Example3;
#elif defined(preuss)
using namespace Preuss;
#endif

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 6; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;             // starting mesh size
    double dT = 0.5;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;
    const  double tfinal = .1;

    // Time integration quadrature
    const size_t quadrature_order_time = 9;
    const QuadratureFormular1d &qTime(*Lobatto(quadrature_order_time));  // specify order of quadrature in time
    const Uint nbTime       = qTime.n;
    const Uint lastQuadTime = nbTime - 1;

    // Space integration quadrature 
    Levelset<2> phi;
    ProblemOption option;
    const int quadrature_order_space       = 9;
    option.order_space_element_quadrature_ = quadrature_order_space;
    AlgoimCutFEM<Mesh, Levelset<2>> convdiff(qTime, phi, option);

    // Global parameters
    const double tau_F_bulk = 1.;    // face stabilization
    const double tau_F_surf = 5.;    // face stabilization
    const double tau_G = 1.;         // interface stabilization
    
    const double delta_bulk = 0.3;   // macro parameter
    const double delta_surf = 0.5;   // macro parameter
    
    std::string ex, method, stab;

#ifdef ex1
    // Paths to store data
    ex = "example1";
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/henry/example1/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/henry/example1/paraview/";
#elif defined(ex2)
    const std::string path_output_data    = "../output_files/henry/example2/data/";
    const std::string path_output_figures = "../output_files/henry/example2/paraview/";
#endif

    // Create directory if not already existent
    // if (MPIcf::IamMaster()) {
    //     std::filesystem::create_directories(path_output_data);
    //     std::filesystem::create_directories(path_output_figures);
    // }

    // Data file to hold problem data
    #ifdef conservative
    method = "conservative";
    #else
    method = "non_conservative";
    #endif

    #ifdef fullstab
    stab = "full";
    #else
    stab = "macro_dsurf_" + std::to_string(delta_surf);
    #endif

    std::ofstream output_data(path_output_data + "data_" + method + "_" + stab + ".dat", std::ofstream::out);
    
    output_data << method << ",\t";
    output_data << stab << ",\t";
    output_data << "tau_F_bulk = " << tau_F_bulk << ",\t tau_F_surf = " << tau_F_surf << ",\t tau_G = " << tau_G << ",\t N = " << quadrature_order_time << ",\t T = " << tfinal << ",\t Example: " << ex;
    output_data << "\n---------------------\n";
    output_data << "h, \t dt,   L2(Omega(T)), L2(Gamma(T)), L2(Omega(t),0,T), L2(Gamma(t),0,T), e_c(T)\n";
    output_data.flush();

    // Data file to hold DOF indices
    std::ofstream indices(path_output_data + "indices.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors, errors_surf, errors_T, errors_T_surf, hs, nxs, nys, dts, omega, gamma,
        global_conservation_errors, reynold_error;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
#if defined(ex1)
        const double x0 = 0., y0 = 0., lx = 1., ly = 1.;
#elif defined(ex2)
        const double x0 = -1., y0 = -1., lx = 2., ly = 2.;
#endif

#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h                          = lx / (nx - 1);
#endif

        const Mesh Th(nx, ny, x0 - Epsilon, y0 - Epsilon, lx, ly);


#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 4;

        dT = h / divisionMeshSize;

        // double dT = tfinal / 32;
        // double dT = 0.5*std::pow(2, -j-1);
        //  double dT = tfinal / (j + 1);

        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;
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

        FESpace Vh(Th, DataFE<Mesh>::P2); // Background FE Space

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly); // FE Space in time
        const Uint ndfTime      = Ih[0].NbDoF();

        // Velocity field
        LagrangeQuad2 FEvelocity(2);
        FESpace VelVh(Th, FEvelocity);
        std::vector<Fun_h> vel(nbTime);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Visualization
        FESpace Lh(Th, DataFE<Mesh>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous, mass_last_previous_surf, mass_initial, mass_initial_surf;
        double intF = 0, int_outflow = 0, intF_surf = 0, intF_total = 0,
               intF_surf_total                = 0; // hold integrals of rhs and Neumann bcs
        double global_conservation_error = 0, local_conservation_error = 0, error_bulk = 0., error_surf = 0.,
               error_I = 0., error_I_surf = 0.;
        std::vector<double> local_conservation_errors;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * time_step;
            const TimeSlab &In(Ih[iter]);
            // const TimeSlab &In_interpolation(Ih_interpolation[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';
            std::cout << "dT = " << dT << '\n';

            // initialization of the interface in each quadrature point

            for (int i = 0; i < nbTime; ++i) {

                R tt  = In.Pt(R1(qTime(i).x));
                phi.t = tt;

                vel[i].init(VelVh, fun_velocity, tt);

                interface.init(i, Th, phi);

                ls[i].init(Lh, fun_levelSet, tt);
            }

            // Create active meshes
            ActiveMesh<Mesh> Thi(Th);

            Thi.truncate(interface, 1);

            //  Cut FE space
            CutSpace Wh(Thi, Vh);

            // Surface active mesh
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            // Surface FE space
            CutFESpace WhGamma(ThGamma, Vh);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Add time slab to surface space on top of the data
            convdiff.add(WhGamma, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);
            FunTest uS(WhGamma, 1), vS(WhGamma, 1);

            // Data for initial solution
            int idx_s0  = Wh.NbDoF() * In.NbDoF();      // index where uS^0 start in the solution array
            int idx_s1  = WhGamma.NbDoF() * In.NbDoF(); // degrees of freedom in surface x time space
            int tot_dof = (int)convdiff.get_nb_dof();   // total degrees of freedom

            Rn data_u0(convdiff.get_nb_dof(), 0.); // initial data total
            convdiff.initialSolution(data_u0);

            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));           // initial data bulk
            KN_<R> data_S0(data_u0(SubArray(WhGamma.NbDoF(), idx_s0))); // initial data surface


            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);
            Fun_h s0h(WhGamma, data_S0);

            // Variational formulation
#if defined(classical)
            
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);
            convdiff.addBilinear(+innerProduct(uS, vS), *interface(0), In, 0);

            // Impose initial condition
            if (iter == 0) {
                convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
                convdiff.addLinearExact(fun_uSurf, +innerProduct(1, vS), *interface(0), In, 0);
            } else {
                convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
                convdiff.addLinear(+innerProduct(s0h.expr(), vS), *interface(0), In, 0);
            }

            convdiff.addBilinear(innerProduct(dt(u), v) + innerProduct(D * grad(u), grad(v)), Thi, In);
            convdiff.addBilinear(innerProduct(dt(uS), vS) + innerProduct(DGamma * gradS(uS), gradS(vS)), interface, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(innerProduct((vel[i].exprList() * grad(u)), v), Thi, In, i);
                convdiff.addBilinear(innerProduct(uS * divS(vel[i]), vS) +
                                         innerProduct((vel[i].exprList() * grad(uS)), vS),
                                     interface, In, i);
            }

#elif defined(conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);
            convdiff.addBilinear(+innerProduct(uS, vS), *interface(lastQuadTime), In, (int)lastQuadTime);

            // Impose initial condition
            if (iter == 0) {
                convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
                convdiff.addLinearExact(fun_uSurf, +innerProduct(1, vS), *interface(0), In, 0);
            } else {
                convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
                convdiff.addLinear(+innerProduct(s0h.expr(), vS), *interface(0), In, 0);
            }

            convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)), Thi, In);
            convdiff.addBilinear(-innerProduct(uS, dt(vS)) + innerProduct(DGamma * gradS(uS), gradS(vS)), interface,
                                 In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Thi, In, i);
                convdiff.addBilinear(-innerProduct(uS, (vel[i].exprList() * grad(vS))), interface, In, i);
            }
#endif

            // Source function
            // convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);
            convdiff.addLinearExact(fun_rhsBulk, +innerProduct(1, v), Thi, In);
            convdiff.addLinearExact(fun_rhsSurf, +innerProduct(1, vS), interface, In);

            // Coupling boundary condition
            convdiff.addBilinear(innerProduct(u - uS, v - vS), interface, In);

#if defined(fullstab)

            // convdiff.addFaceStabilization(
            //     + innerProduct(h * tau_F_bulk * jump(grad(u) * n), jump(grad(v) * n)) 
            //     + innerProduct(h * h * h * tau_F_bulk * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In);

            convdiff.addPatchStabilization(
                + innerProduct(tau_F_bulk / h / h * jump(u), jump(v))
                , Thi
                , In
            );

            // convdiff.addFaceStabilization(
            //     + innerProduct(tau_F_surf * jump(grad(uS) * n), jump(grad(vS) * n)) 
            //     + innerProduct(h * h * tau_F_surf * jump(grad(grad(uS) * n) * n), jump(grad(grad(vS) * n) * n)),
            //     ThGamma, In);

            convdiff.addPatchStabilization(
                + innerProduct(tau_F_surf / h / h / h * jump(uS), jump(vS))
                , ThGamma
                , In
            );

            // convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(u), jump(v)), Thi, In);
            // convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(uS), jump(vS)), ThGamma, In);

#elif defined(macro)
            // MacroElementPartition<Mesh> TimeMacro(Thi, 0.3);
            // std::cout << TimeMacro.number_of_stabilized_edges << "\n";
            // std::cout << "number of stabilized edges: " << convdiff.get_number_of_stabilized_edges() << "\n";

            AlgoimMacro<Mesh, Levelset<2>> TimeMacro(Thi, 0.4, phi, In, qTime);
            TimeMacro.findSmallElement();
            TimeMacro.createMacroElement();
            TimeMacro.setInnerEdges();

            AlgoimMacroSurface<Mesh, Levelset<2>> TimeMacroSurf(ThGamma, 0.5, phi, In, qTime);
            TimeMacroSurf.findSmallElement();
            TimeMacroSurf.createMacroElement();
            TimeMacroSurf.setInnerEdges();

            // convdiff.addFaceStabilization(
            //     + innerProduct(h * tau_F * jump(grad(u) * n), jump(grad(v) * n)) 
            //     + innerProduct(h * h * h * tau_F * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n))
            //     , Thi
            //     , In
            //     , TimeMacro);

            convdiff.addPatchStabilization(
                + innerProduct(tau_F_bulk / h / h * jump(u), jump(v))
                , Thi
                , In
                , TimeMacro
            );

            // convdiff.addFaceStabilization(
            //     + innerProduct(tau_F * jump(grad(uS) * n), jump(grad(vS) * n)) 
            //     + innerProduct(h * h * tau_F * jump(grad(grad(uS) * n) * n), jump(grad(grad(vS) * n) * n))
            //     , ThGamma
            //     , In
            //     , TimeMacroSurf);

            convdiff.addPatchStabilization(
                + innerProduct(tau_F_surf / h / h / h * jump(uS), jump(vS))
                , ThGamma
                , In
                , TimeMacroSurf
            );

            if (iterations == 1 && h > 0.1) {
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
            }
#endif

            convdiff.addBilinear(+ innerProduct(tau_G * grad(uS) * n, grad(vS) * n)
                        + innerProduct(tau_G * h * h * grad(grad(uS) * n) * n, grad(grad(vS) * n) * n), interface, In);


            if (iter == total_number_iteration - 1) {

                matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");

                indices << idx_s0 << ", " << idx_s1
                        << "\n"; // export indices needed for computing the scaled condition number
                indices.flush(); // export simultaneously
            }

            // Solve linear system
            convdiff.solve("mumps");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Fun_h uh_t(Wh, In, data_u0); // FEM function in Pk(In) x Lagrange_m(Omega)
            // error_I += L2_norm_T(uh_t, fun_uBulkD, Thi, In, qTime, phi,
            //                      quadrature_order_space); // int_In ||u(t) - u_h(t)||_{Omega(t)} dt
            Rn data_uh_t(data_u0(SubArray(idx_s0, 0)));
            Fun_h uh_t(Wh, In, data_uh_t); // FEM function in Pk(In) x Lagrange_m(Gamma)
            error_I += L2_norm_T(uh_t, fun_uBulkD, Thi, In, qTime, phi,
                                 quadrature_order_space); // int_In ||u(t) - u_h(t)||_{Gamma(t)} dt

            //Rn data_uh_surf_t(WhGamma.get_nb_dof(), 0.);
            Rn data_uh_surf_t(data_u0(SubArray(idx_s1, idx_s0)));
            Fun_h uh_surf_t(WhGamma, In, data_uh_surf_t); // FEM function in Pk(In) x Lagrange_m(Gamma)
            error_I_surf += L2_norm_surf_T(uh_surf_t, fun_uSurf, interface, In, qTime, phi,
                                 quadrature_order_space); // int_In ||u(t) - u_h(t)||_{Gamma(t)} dt


            // Compute error in Reynold relation
            // {

            //     AlgoimCutFEM<Mesh, Levelset<2>> reynold(qTime, phi);

            //     reynold.initSpace(Wh, In);

            //     reynold.addBilinear(innerProduct(u, v), Thi, (int)lastQuadTime, In);
            //     if (iter == 0)
            //         reynold.addLinearExact(fun_uBulk, innerProduct(1, v), Thi, 0, In);
            //     else
            //         reynold.addLinear(innerProduct(b0h.expr(), v), Thi, 0, In);

            //     reynold.addBilinear(-innerProduct(dt(u), v) - innerProduct(u, dt(v)), Thi, In);

            //     for (int i = 0; i < nbTime; ++i) {
            //         reynold.addBilinear(-innerProduct((vel[i].exprList() * grad(u)), v) -
            //                                 innerProduct(u, (vel[i].exprList() * grad(v))),
            //                             Thi, In, i);
            //     }

            //     int N = Wh.NbDoF();
            //     Rn lhs(ndfTime * N);
            //     multiply(ndfTime * N, ndfTime * N, reynold.mat_, data_u0, lhs);

            //     lhs -= reynold.rhs_;

            //     reynold_error.at(j) = lhs.linfty();

            //     std::cout << " e_r^n = " << reynold_error.at(j) << '\n';
            // }

            // Compute area of domain in time quadrature point 0
            Fun_h funone(Wh, fun_one);

            double intGamma = integral_algoim(funone, In, interface, phi, 0, quadrature_order_space) / dT;
            double intOmega = integral_algoim(funone, Thi, phi, In, qTime, lastQuadTime, quadrature_order_space);

            gamma[j] = intGamma;
            omega[j] = intOmega;

            // Compute error of numerical solution
            Rn sol(Wh.get_nb_dof(), 0.);
            Fun_h funuh(Wh, sol);

            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));

            for (int n = 1; n < ndfTime; n++) {
                sol += data_u0(SubArray(Wh.get_nb_dof(), n * Wh.get_nb_dof()));
            }

            error_bulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1, quadrature_order_space);
            std::cout << " t_n -> || u-uex||_Omega(T) = " << error_bulk << '\n';

            Rn solsurf(WhGamma.get_nb_dof(), 0.);
            Fun_h funuh_surf(WhGamma, solsurf);

            solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0));

            for (int n = 1; n < ndfTime; n++) {
                solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0 + n * WhGamma.get_nb_dof()));
            }

            error_surf = L2_norm_surface(funuh_surf, fun_uSurf, *interface(lastQuadTime), current_time + dT, phi, 0, 1);
            std::cout << " t_n -> || u-uex||_Gamma(T) = " << error_surf << '\n';

            errors[j]      = error_bulk;
            errors_surf[j] = error_surf;

// Compute conservation error
// if (iterations == 1 && h > 0.01) {
#if defined(conservation)
            intF = integral_algoim(fun_rhsBulk, 0, Thi, phi, In, qTime,
                                   quadrature_order_space); // integrate source over In
            intF_surf = integral_algoim(fun_rhsSurf, In, interface, phi, 0,
                                   quadrature_order_space); // integrate flux boundary over In

            intF_total += intF;
            intF_surf_total += intF_surf;

            // intF = integral_algoim(f, 0, Thi, phi, In, qTime);        // integrate source over In
            // intG = integral_algoim(g_Neumann, In, interface, phi, 0); // integrate flux across boundary over In

            double mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime,
                                               quadrature_order_space); // mass in last quad point

            double mass_last_surf = integral_algoim(funuh_surf, *interface(lastQuadTime), 0, phi, In, qTime,
                                                    lastQuadTime, quadrature_order_space);

            if (iter == 0) {
                mass_initial       = integral_algoim(fun_uBulk, Thi, phi, In, qTime, 0, quadrature_order_space);
                mass_initial_surf = integral_algoim(fun_uSurf, *interface(0), 0, phi, In, qTime,
                                                    0, quadrature_order_space);
                mass_last_previous = mass_initial;
                mass_last_previous_surf = mass_initial_surf;
                // mass_last_previous = integral_algoim(b0h, Thi, phi, In, qTime, 0);
            }

            local_conservation_error  = (mass_last + mass_last_surf - mass_last_previous - mass_last_previous_surf - intF - intF_surf);
            // global_conservation_error += local_conservation_error;
            global_conservation_error = (mass_last + mass_last_surf - mass_initial - mass_initial_surf - intF_total - intF_surf_total);

            std::cout << "global_conservation_error: " << global_conservation_error << "\n";
            std::cout << "local_conservation_error: " << local_conservation_error << "\n";

            // output_data << std::setprecision(10);
            // output_data << current_time << "," << (mass_last - mass_last_previous) << "," << intF << "," << intG <<
            // ","
            //            << local_conservation_error << '\n';

            mass_last_previous = mass_last; // set current last to previous last for next time slab
                                            //}
            mass_last_previous_surf = mass_last_surf;

            global_conservation_errors[j] = std::fabs(global_conservation_error);
            local_conservation_errors.push_back(std::fabs(local_conservation_error));

#endif

            if ((iterations == 1) && (h > 0.01)) {
                Fun_h sol_h(Wh, sol);
                Paraview<Mesh> writerTh(Th, path_output_figures + "Th.vtk");
                writerTh.add(ls[0], "levelSet0", 0, 1);
                writerTh.add(ls[1], "levelSet1", 0, 1);
                writerTh.add(ls[2], "levelSet2", 0, 1);
                
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) + ".vtk");
                Paraview<Mesh> writer_surf(ThGamma, path_output_figures + "surf_" + std::to_string(iter + 1) + ".vtk");
                // writer.add(b0h, "bulk", 0, 1);
                // writer.add(sol_h, "bulk_end", 0, 1);

                writer.add(funuh, "bulk", 0, 1);
                writer_surf.add(funuh_surf, "surf", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                
                Fun_h uSex(WhGamma, fun_uSurf, current_time);
                Fun_h fS(WhGamma, fun_rhsSurf, current_time);

                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);

                writer_surf.add(uSex, "surf_exact", 0, 1);
                writer_surf.add(fS, "surf_rhs", 0, 1);

                //writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                writer.add(fabs(funuh.expr() - uBex.expr()), "bulk_error");
                writer_surf.add(fabs(funuh_surf.expr() - uSex.expr()), "surf_error");
                
                
                writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writer.writeFaceStab(Thi, 0, path_output_figures + "Edges" + std::to_string(iter + 1) + ".vtk");

                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 0, 2,
                                             path_output_figures + "AlgoimQuadrature_0_" + std::to_string(iter + 1) +
                                                 ".vtk");
                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 1, 2,
                                             path_output_figures + "AlgoimQuadrature_1_" + std::to_string(iter + 1) +
                                                 ".vtk");
                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 2, 2,
                                             path_output_figures + "AlgoimQuadrature_2_" + std::to_string(iter + 1) +
                                                 ".vtk");
            }

            // if (iterations > 1 && iter == total_number_iteration - 1)
            //     output_data << h << "," << dT << "," << errBulk << '\n';

            iter++;
        }

        errors_T[j] = std::sqrt(error_I);
        errors_T_surf[j] = std::sqrt(error_I_surf);

        output_data << h << "\t" << dT << "\t" << error_bulk << "\t" << error_surf << "\t" << errors_T[j] << "\t" << errors_T_surf[j] << "\t" << global_conservation_errors[j] << "\n";
        output_data.flush();

        std::cout << "error_T = " << errors_T[j] << "\n";
        std::cout << "error_T_surf = " << errors_T_surf[j] << "\n";

        // std::cout << "\n";
        // std::cout << "Local conservation error = [";
        // for (auto &err : local_conservation_errors) {

        //     std::cout << err;

        //     std::cout << ", ";
        // }
        // std::cout << "]\n";

// Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_t)
        dT *= 0.5;
#elif defined(use_h)
        h *= 0.5;
        // h *= sqrt(0.5);
#endif
    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "Errors Bulk = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_surf.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors In= [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_T.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors In Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_T_surf.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Global Conservation Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << global_conservation_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Reynold error = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << reynold_error.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "|Gamma_h| = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << gamma.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "|Omega_h| = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << omega.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "dT = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << dts.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << '\n';

    std::cout << "nx = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nxs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << '\n';

    std::cout << "ny = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nys.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}
