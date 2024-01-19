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
 * @note We consider a time-dependent bulk problem on Omega1 or Omega2.

 *  Problems:

    * Omega_1(t):
    Find u in Omega_1(t) such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    in Omega_1(t).
                                        BC1,  on Gamma(t),
                                        BC2,  on \partial\Omega (outer boundary).


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
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

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

    // Correct sign of normal vector at interface
    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (x - 1. / 2 - 0.28 * sin(pi * t))) /
               (250 * sqrt((pow(x - 1. / 2 - 0.28 * sin(pi * t), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2)))) -
           (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (y - 0.5 + 0.28 * cos(pi * t))) /
               (250 * sqrt((pow((x - 1. / 2 - 0.28 * sin(pi * t)), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2))));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

R fun_one(double *P, const int i) { return 1.; }

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

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
           (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 0.5)) / 5 +
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 0.5)) / 5;
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

namespace Example1_pure_diffusion {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (x - 1. / 2)) /
               (250 * sqrt((pow((x - 1. / 2), 2) + pow(y - 11. / 50, 2)))) -
           (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (y - 11. / 50)) /
               (250 * sqrt((pow((x - 1. / 2), 2) + pow(y - 11. / 50, 2))));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return 0.;
    else
        return 0.;
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

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * cos(pi * x) * cos(pi * y) * (100 * sin(2 * pi * t) - pi * cos(2 * pi * t))) / 125;
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

} // namespace Example1_pure_diffusion

namespace Example1_Omega1 {
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
    // return -sqrt((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.22)*(P[1]-0.22)) - 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17);
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const {
        return -((x(0) - 0.5 - 0.28 * sin(pi * t)) * (x(0) - 0.5 - 0.28 * sin(pi * t)) +
                 (x(1) - 0.5 + 0.28 * cos(pi * t)) * (x(1) - 0.5 + 0.28 * cos(pi * t)) - 0.17 * 0.17);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(-2.0 * (x(0) - 0.5 - 0.28 * sin(M_PI * t)),
                                     -2.0 * (x(1) - 0.5 + 0.28 * cos(M_PI * t)));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(double *P) const {
        R norm = sqrt(pow(2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))), 2) +
                      pow(2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))), 2));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(-2.0 * (P[0] - (0.5 + 0.28 * sin(M_PI * t))) / norm,
                  -2.0 * (P[1] - (0.5 - 0.28 * cos(M_PI * t))) / norm);
    }
};

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // Wrong sign of normal vector at interface
    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (x - 1. / 2 - 0.28 * sin(pi * t))) /
               (250 * sqrt((pow(x - 1. / 2 - 0.28 * sin(pi * t), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2)))) -
           (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (y - 0.5 + 0.28 * cos(pi * t))) /
               (250 * sqrt((pow((x - 1. / 2 - 0.28 * sin(pi * t)), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2))));
}

R fun_one(double *P, const int i) { return 1.; }

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

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * cos(M_PI * P[0]) * cos(M_PI * P[1]) * cos(2 * M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
           (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 0.5)) / 5 +
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 0.5)) / 5;
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

} // namespace Example1_Omega1

namespace Lehrenfeld {

double fun_one(double *P, const int i) { return 1.; }

double r0 = 1. + Epsilon;

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    return (x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y - r0 * r0;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - r0 * r0; }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] - (1 - P[1] * P[1]) * t) * (P[0] - (1 - P[1] * P[1]) * t) + P[1] * P[1] - r0 * r0;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2. * (x(0) - (1 - x(1) * x(1)) * t),
                                     4. * t * x(1) * (x(0) - t * (1 - x(1) * x(1))) + 2. * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(pow(2. * (P[0] - (1 - P[1] * P[1]) * t), 2) +
                      pow(4. * t * P[1] * (P[0] - t * (1 - P[1] * P[1])) + 2. * P[1], 2));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2. * (P[0] - (1 - P[1] * P[1]) * t) / norm,
                  (4. * t * P[1] * (P[0] - t * (1 - P[1] * P[1])) + 2. * P[1]) / norm);
    }
};

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return 0.;
    return -(pi * sin(pi * t) * sin(pi * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
             (2 * x + 2 * t * (y * y - 1)) * (2 * x + 2 * t * (y * y - 1))) /
               (200 * sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) -
           (pi * sin(pi * t) * sin(pi * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (2 * y + 4 * t * y * (x + t * (y * y - 1))) * (2 * y + 4 * t * y * (x + t * (y * y - 1)))) /
               (200 * sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return 1 - P[1] * P[1];
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double r0 = 1., x = P[0], y = P[1];
    // return cos(M_PI*sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) +
    // y*y)/r0)*sin(M_PI*t);
    return cos(M_PI * ((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) / (r0 * r0)) * sin(M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double r0 = 1., x = P[0], y = P[1];
    // return cos(M_PI*sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) +
    // y*y)/r0)*sin(M_PI*t);
    return cos(M_PI * ((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) / (r0 * r0)) * sin(M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return M_PI * cos(M_PI * t) * cos(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) +
           (M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y))) / 50 +
           (M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (8 * t * t * y * y + 4 * t * (x + t * (y * y - 1)) + 2)) /
               100 +
           (M_PI * M_PI * cos(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * sin(M_PI * t) *
            (2 * x + 2 * t * (y * y - 1)) * (2 * x + 2 * t * (y * y - 1))) /
               100 +
           (M_PI * M_PI * cos(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * sin(M_PI * t) *
            (2 * y + 4 * t * y * (x + t * (y * y - 1))) * (2 * y + 4 * t * y * (x + t * (y * y - 1)))) /
               100 +
           M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * (y * y - 1) *
               (2 * x + 2 * t * (y * y - 1)) -
           2 * M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
               (y * y - 1) * (x + t * (y * y - 1));
}

} // namespace Lehrenfeld

namespace Kite {
double R0 = .5;
double C  = 5. / 3;
double W  = 1. / 6;

double fun_one(double *P, const int i) { return 1.; }

double rho(double *P, const double t) { return (W - C * P[1] * P[1]) * t; }

double fun_levelSet(double *P, const int i, const double t) {
    return (P[0] - rho(P, t)) * (P[0] - rho(P, t)) + P[1] * P[1] - R0 * R0;
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] - (W - C * P[1] * P[1]) * t) * (P[0] - (W - C * P[1] * P[1]) * t) + P[1] * P[1] - R0 * R0;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2. * (x(0) - (W - C * x(1) * x(1)) * t),
                                     4. * C * t * x(1) * (x(0) - t * (W - C * x(1) * x(1))) + 2. * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(pow(2. * (P[0] - (W - C * P[1] * P[1]) * t), 2) +
                      pow(4. * C * t * P[1] * (P[0] - t * (W - C * P[1] * P[1])) + 2. * P[1], 2));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2. * (P[0] - (W - C * P[1] * P[1]) * t) / norm,
                  (4. * C * t * P[1] * (P[0] - t * (W - C * P[1] * P[1])) + 2. * P[1]) / norm);
    }
};

R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * cos(pi * t) +
           (pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t)) /
               (50 * (R0 * R0)) +
           (pi * pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
            pow((2 * x - 2 * t * (-C * y * y + W)), 2)) /
               (100 * (R0 * R0 * R0 * R0)) +
           (pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
            (8 * C * C * t * t * y * y + 4 * C * t * (x - t * (-C * y * y + W)) + 2)) /
               (100 * (R0 * R0)) +
           (pi * pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
            pow((2 * y + 4 * C * t * y * (x - t * (-C * y * y + W))), 2)) /
               (100 * (R0 * R0 * R0 * R0)) +
           (2 * pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
            (-C * y * y + W) * (x - t * (-C * y * y + W))) /
               (R0 * R0) -
           (pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) * (-C * y * y + W) *
            (2 * x - 2 * t * (-C * y * y + W))) /
               (R0 * R0);
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;

    return -(pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
             pow((2 * y + 4 * C * t * y * (x - t * (-C * y * y + W))), 2)) /
               (200 * (R0 * R0) * sqrt((pow((x - t * (-C * y * y + W)), 2) + y * y))) -
           (pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
            pow((2 * x - 2 * t * (-C * y * y + W)), 2)) /
               (200 * (R0 * R0) * sqrt((pow((x - t * (-C * y * y + W)), 2) + y * y)));
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return W - C * P[1] * P[1];
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    double r = (x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y;
    return cos(M_PI * (r) / (R0 * R0)) * sin(M_PI * t);

    // double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    // return cos(M_PI * (r + Epsilon) / (R0));
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double r = (x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y;
    return cos(M_PI * (r) / (R0 * R0)) * sin(M_PI * t);

    // double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    // return cos(M_PI * (r + Epsilon) / (R0));
}
} // namespace Kite

namespace Example4 {
// A circle moving in a circular trajectory
const double D        = 1.;
const double R0       = 0.17;
const double beta_max = M_PI * 0.5;

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - R0 * R0);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - R0 * R0);
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

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

R fun_one(double *P, const int i) { return 1.; }

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
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    double xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);

    double r   = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
    return cos(M_PI * r / R0) * sin(M_PI * t);

}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);

    double r   = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
    return cos(M_PI * r / R0) * sin(M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // automatic
    return D*((3.141592653589793*sin((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8))*sin(t*3.141592653589793)*1.0/sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0)*2.0E+8)/R0+(1.0/(R0*R0)*(3.141592653589793*3.141592653589793)*sin(t*3.141592653589793)*cos((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8))*pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12)/(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0)-(3.141592653589793*sin((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8))*sin(t*3.141592653589793)*pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*1.0/pow(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0,3.0/2.0)*4.0E+20)/R0-(3.141592653589793*sin((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8))*sin(t*3.141592653589793)*pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*1.0/pow(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0,3.0/2.0)*4.0E+20)/R0+(1.0/(R0*R0)*(3.141592653589793*3.141592653589793)*sin(t*3.141592653589793)*cos((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8))*pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12)/(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))+3.141592653589793*cos(t*3.141592653589793)*cos((3.141592653589793*sqrt(pow(y*5.0E+1+cos(t*3.141592653589793)*1.4E+1-2.5E+1,2.0)*4.0E+12+pow(x*-5.0E+1+sin(t*3.141592653589793)*1.4E+1+2.5E+1,2.0)*4.0E+12+1.0))/(R0*1.0E+8));

}

} // namespace Example4


// The below parameters can be varied according to the options to use different methods,
// different numerical examples, different subdomains and boundary conditions etc.

//* Quadrature method (options: "algoim", "quad", "triangle")
#define algoim

//* Set numerical example (options: "example1", "lehrenfeld")
#define example4
//* Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominatednot

//* Choose domain to solve on (options: "omega1", "omega2")
#define omega2
// If "omega1":
// Set type of BCs on outer boundary (options: "dirichlet1" or "neumann1")
#define dirichlet1
// Set type of BCs on interface (options: "dirichlet2" or "neumann2")
#define dirichlet2

// If "omega2":
// Set type of BCs on interface (options: "dirichlet", "neumann")
#define neumann
//* Set scheme for the method (options: "classical", "conservative")
#define conservative
//* Set stabilization method (options: "fullstab", "macro")
#define fullstab
//* Decide whether to solve for level set function, or to use exact (options:
// "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h    // to set mesh size using the h parameter. Write use_n to decide
                 // using nx, ny.
#define use_tnot // write use_t to control dT manually. Otherwise it is set
                 // proportional to h.

// Setup two-dimensional class types
const int d = 2;

typedef MeshQuad2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

// Do not touch
#ifdef example1
#ifdef omega1
using namespace Example1_Omega1;
#elif defined(omega2)
using namespace Example1;
#endif
#elif defined(lehrenfeld)
using namespace Lehrenfeld;
#elif defined(lehrenfeld2)
using namespace Lehrenfeld2;
#elif defined(kite)
using namespace Kite;
#elif defined(example4)
using namespace Example4;
#endif

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;            // starting mesh size
    double dT = 0.1;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#ifdef example1 
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/bulk/example1/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/bulk/example1/";
#elif defined(lehrenfeld)
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/bulk/lehrenfeld/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/bulk/lehrenfeld/";
#elif defined(lehrenfeld2)
    // Paths to store data
    const std::string path_output_data = "/NOBACKUP/smyrback/output_files/bulk/lehrenfeld2/";
    // const std::string path_output_figures = "../output_files/bulk/algoim/lehrenfeld2/paraview/";
#elif defined(kite)
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/bulk/kite/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/bulk/kite/";
#elif defined(example4)
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/bulk/example1/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/bulk/example1/";
#endif

    // Create directory if not already existent
    // if (MPIcf::IamMaster()) {
    //     std::filesystem::create_directories(path_output_data);
    //     std::filesystem::create_directories(path_output_figures);
    // }

    // Data file to hold problem data
    std::ofstream outputData(path_output_data + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors; // array to hold bulk errors
    std::array<double, iterations> hs;     // array to hold mesh sizes
    std::array<double, iterations> nxs;    // array to hold mesh sizes
    std::array<double, iterations> nys;    // array to hold mesh sizes
    std::array<double, iterations> dts;
    std::array<double, iterations> omega;
    std::array<double, iterations> gamma;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
#if defined(example1)
        const double lx = 1., ly = 1.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0. - Epsilon, 0. - Epsilon, lx, ly);
        // Mesh Th(12, 11, 0.249, 0. - Epsilon, 0.55, 0.45);
#elif defined(lehrenfeld)
        const double lx = 7., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3.5, -1.5, lx, ly);
        // Mesh Th(31, 31, -1.5, -1.5, 3., ly);

#elif defined(lehrenfeld2)
        const double lx = 2. + Epsilon, ly = 2. + Epsilon;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -1.12, -1.13, lx, ly);

#elif defined(kite)
        const double lx = 7., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3.523, -1.5108, lx, ly);

#elif defined(example4)
        const double lx = 1., ly = 1.;
        const double x0 = 0. - Epsilon, y0 = 0. - Epsilon;
        //const double lx = 4. * Example4::R0, ly = 4. * Example4::R0;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, x0, y0, lx, ly);
        //Mesh Th(nx, ny, 0. - 2.01 * Example4::R0, -2.01 * Example4::R0, lx, ly);


#endif

        // Parameters
        const double tfinal = .1; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 3;

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

#if defined(convection_dominated)
        const double D = .01;
#else
        const double D         = 1;
#endif

        const double lambda = 1.; // Nitsche's method penalty parameter

        // CG stabilization parameter
        const double tau1 = 0.1, tau2 = 0.1;

        FESpace Vh(Th, DataFE<Mesh>::P2);               // Background FE Space
        FESpace Vh_interpolation(Th, DataFE<Mesh>::P2); // for interpolating data

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly);               // FE Space in time
        FESpace1 Ih_interpolation(Qh, DataFE<Mesh1>::P2Poly); // for interpolating data

        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(7)); // specify order of quadrature in time
        // const QuadratureFormular1d &qTime(QF_Euler); // specify order of quadrature in time
        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Velocity field
        LagrangeQuad2 FEvelocity(2);

        FESpace VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Visualization
        FESpace Lh(Th, DataFE<Mesh>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

        Levelset<2> phi;
        AlgoimCutFEM<Mesh, Levelset<2>> convdiff(qTime, phi);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous;
        double intF = 0, int_outflow = 0, intG = 0; // hold integrals of rhs and Neumann bcs
        double errBulk = 0.;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * time_step;
            const TimeSlab &In(Ih[iter]);
            const TimeSlab &In_interpolation(Ih_interpolation[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';
            std::cout << "dT = " << dT << '\n';

            // initialization of the interface in each quadrature point

            for (int i = 0; i < nbTime; ++i) {

                R tt  = In.Pt(R1(qTime(i).x));
                phi.t = tt;
                interface.init(i, Th, phi);

                ls[i].init(Lh, fun_levelSet, tt);
            }

            // Create active meshes
            ActiveMesh<Mesh> Thi(Th);

            // #ifdef omega1
            Thi.truncate(interface, 1); // remove part with negative sign of level
                                        // #elif defined(omega2)
            //             Thi.truncate(interface, 1); // remove part with positive sign of level
            //  set to get inner domain
            // #endif
            //  Cut FE space
            CutSpace Wh(Thi, Vh);

            // AlgoimCutFEM<Mesh, Levelset<2>> initial_condition(Wh, phi);

            // initial_condition.initSpace(Wh, In);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // gnuplot::save(Th);
            // gnuplot::save<Mesh, Levelset<2>>(Thi, *interface(0), phi, "interface.dat", current_time);
            // gnuplot::save<Mesh>(*interface(0));
            // getchar();

            // Right hand side functions
            Fun_h f(Vh_interpolation, In_interpolation, fun_rhsBulk);
            Fun_h g(Vh_interpolation, In_interpolation, fun_uBulk); // create an FE-function of the exact bulk
                                                                    // solution Omega2
            // Fun_h g_Neumann(Vh, In, fun_neumann_Gamma); // computer Neumann BC
            Fun_h g_Neumann(Vh_interpolation, In_interpolation, fun_neumann_Gamma); // computer Neumann BC
            // Fun_h g_Neumann_left(Wh, In, fun_neumann_left);     // label 1
            // Fun_h g_Neumann_bottom(Wh, In, fun_neumann_bottom); // label 4
            // Fun_h g_Neumann_right(Wh, In, fun_neumann_right);   // label 2
            // Fun_h g_Neumann_top(Wh, In, fun_neumann_top);       // label 3

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);

            // Data for initial solution
            Rn data_u0(convdiff.get_nb_dof(), 0.); // initial data total
            convdiff.initialSolution(data_u0);
            // initial_condition.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0))); // initial data bulk

            if (iter == 0)
                interpolate(Wh, data_B0, fun_uBulkInit);

            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);

            // // Solve for initial condition
            // FunTest s(Wh, 1), r(Wh, 1);
            // initial_condition.addBilinearAlgoim(+innerProduct(s, r), Thi, 0, In);
            // initial_condition.addFaceStabilization(+innerProduct(0.01 * h * jump(grad(s) * n), jump(grad(r) * n)),
            //                                        Thi, In);
            // initial_condition.addLinearAlgoim(+innerProduct(b0h.expr(), r), Thi, 0, In);
            // initial_condition.solve("mumps");
            // data_u0 = initial_condition.rhs_;

            // Plot initial solution in paraview
            // #ifdef USE_MPI
            // if (iter == 0 && MPIcf::IamMaster()) {
            // #else

            // if (iter == 0) {
            //     // #endif
            //     Paraview<Mesh> writerInitial(Thi, path_output_figures + "BulkInitial.vtk");
            //     writerInitial.add(b0h, "bulk", 0, 1);

            //     // Add exact solutions
            //     Fun_h uBex(Wh, fun_uBulkD, 0.);
            //     Fun_h uRhs(Wh, fun_rhsBulk, 0.);

            //     writerInitial.add(uBex, "bulk_exact", 0, 1);
            //     writerInitial.add(uRhs, "bulk_rhs", 0, 1);
            //     // writerInitial.add(ls[0], "levelSet", 0, 1);
            // }

            //** Assembling linear and bilinear forms

            //* Time terms
#if defined(classical)

            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            convdiff.addBilinear(+innerProduct(dt(u), v) + innerProduct(D * grad(u), grad(v)) +
                                     innerProduct((vel.exprList() * grad(u)), v),
                                 Thi, In);

#elif defined(conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)) -
                                             innerProduct(u, (vel.exprList() * grad(v))),
                                         Thi, In);

            
            const double C = 0.1;
            double tau_supg = C*h*h/D;  

            // convdiff.addBilinear(innerProduct(dt(u) + vel.exprList()*grad(u) - D*div(grad(u)), tau_supg*(dt(v) + vel.exprList() * grad(v))), Thi, In);
            
            // convdiff.addLinear(innerProduct(f.expr(), tau_supg*(dt(v) + vel.exprList() * grad(v))), Thi, In);
            
#endif

            convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);

            //* Stabilization
            double stab_bulk_faces = tau1 * h;
            double stab_mass       = 0.; // tau1 * h;
            double stab_dt         = 0.; // tau1 * h;

#if defined(fullstab)
            // convdiff.addFaceStabilization(+innerProduct(stab_bulk_faces * jump(grad(u) * n), jump(grad(v) * n))
            //                               + innerProduct(tau1 * h*h*h*jump(grad(grad(u) * n)*n), jump(grad(grad(v) *
            //                               n)*n))
            // , Thi, In);

            // OLD
            // convdiff.addFaceStabilization(
            //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
            //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In);

            // NEW
            convdiff.addPatchStabilization(+innerProduct(tau1 / h / h * jump(u), jump(v)), Thi, In);

            // convdiff.addBilinear(+innerProduct(h * h * tau1 * grad(u) * n, grad(v) * n), interface, In);

            // double ccend = 1. / In.T.mesure() * 1. / qTime[lastQuadTime].a;
            // convdiff.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), ccend * jump(grad(v) * n)),
            // Thi,
            //                               In, lastQuadTime);

            // convdiff.addFaceStabilization(-innerProduct(stab_dt * jump(grad(u) * n), jump(grad(dt(v)) * n)), Thi,
            // In);

#elif defined(macro)
            MacroElementPartition<Mesh> TimeMacro(Thi, 0.45);
            convdiff.addFaceStabilization(
                +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
                Thi, In, TimeMacro);

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
            }
#endif

            //* Boundary conditions
// If Omega1
#ifdef omega1

//* Dirichlet on both outer and inner boundary
#if defined(dirichlet1) && defined(dirichlet2)
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

            // Dirichlet inner
            convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v)  // added penalty
                                 ,
                                 interface, In);

            convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                               interface, In);

//* Dirichlet on outer and Neumann on inner
#elif defined(dirichlet1) && defined(neumann2)
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

            // Neumann inner
            convdiff.addLinear(-innerProduct(g_Neumann.expr(), v), interface,
                               In); // note the negative sign because of the changed interface normal

            //* Neumann on outer and Dirichlet on inner
#elif defined(neumann1) && defined(dirichlet2)
            // Neumann outer

            // In Example 1, the Neumann function is zero on the outer boundary
            // convdiff.addLinear(innerProduct(g_Neumann.expr(), v), Thi, INTEGRAL_BOUNDARY, In);
            // convdiff.addLinear(innerProduct(g_Neumann_left.expr(), v), Thi, INTEGRAL_BOUNDARY, In,
            // (std::list<int>){4}); convdiff.addLinear(innerProduct(g_Neumann_bottom.expr(), v), Thi,
            // INTEGRAL_BOUNDARY, In, (std::list<int>){1}); convdiff.addLinear(innerProduct(g_Neumann_right.expr(), v),
            // Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){3}); convdiff.addLinear(innerProduct(g_Neumann_top.expr(),
            // v), Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){2});

#if defined(conservative) // extra term arises when using Reynold's transport theorem
            convdiff.addBilinear(innerProduct(vel.exprList() * u * n, v), Thi, INTEGRAL_BOUNDARY, In);
#endif

            // Dirichlet inner
            convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v)  // added penalty
                                 ,
                                 interface, In);

            convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                               interface, In);

//* Neumann on outer and Neumann on inner
#elif defined(neumann1) && defined(neumann2)
            // Neumann outer
            // convdiff.addLinear(innerProduct(g_Neumann_left.expr(), v), Thi, INTEGRAL_BOUNDARY, In,
            // (std::list<int>){1}); convdiff.addLinear(innerProduct(g_Neumann_bottom.expr(), v), Thi,
            // INTEGRAL_BOUNDARY, In, (std::list<int>){4}); convdiff.addLinear(innerProduct(g_Neumann_right.expr(), v),
            // Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){2}); convdiff.addLinear(innerProduct(g_Neumann_top.expr(),
            // v), Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){3});

#if defined(conservative)
            convdiff.addBilinear(innerProduct(vel.exprList() * u * n, v), Thi, INTEGRAL_BOUNDARY, In);
#endif
            // Neumann inner
            convdiff.addLinear(-innerProduct(g_Neumann.expr(), v), interface, In);

#endif

// If Omega2
#elif defined(omega2)
#ifdef neumann

            convdiff.addLinear(+innerProduct(g_Neumann.expr(), v), interface, In);

#elif defined(dirichlet)

            //* Nitsche's method:

            // LHS terms
            convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v)  // added penalty
                                 ,
                                 interface, In);

            // RHS terms
            convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                               interface, In);

#endif
#endif

            if (iter == total_number_iteration - 1) {
                matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
            }

            // Solve linear system
            convdiff.solve("mumps");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Compute area of domain in time quadrature point 0
            Fun_h funone(Wh, fun_one);

            double intGamma = integral_algoim(funone, In, interface, phi, 0) / dT;
            double intOmega = integral_algoim(funone, Thi, phi, In, qTime, lastQuadTime);
            gamma.at(j)     = intGamma;
            omega.at(j)     = intOmega;

            // Compute error of numerical solution
            Rn sol(Wh.get_nb_dof(), 0.);
            Fun_h funuh(Wh, sol);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));

            // errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, 0, phi, 0, 1);
            // std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << '\n';

            for (int n = 1; n < ndfTime; n++) {
                sol += data_u0(SubArray(Wh.get_nb_dof(), n * Wh.get_nb_dof()));
            }

            errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';
            errors.at(j) = errBulk;

            // Compute conservation error
            if (iterations == 1) {
                intF = integral_algoim(f, 0, Thi, phi, In, qTime);        // integrate source over In
                intG = integral_algoim(g_Neumann, In, interface, phi, 0); // integrate flux across boundary over In

                double mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime); // mass in last quad point

                if (iter == 0) {
                    mass_last_previous = integral_algoim(b0h, Thi, phi, In, qTime, 0);
                }

                outputData << std::setprecision(10);
                outputData << current_time << "," << (mass_last - mass_last_previous) << "," << intF << "," << intG
                           << "," << ((mass_last - mass_last_previous) - intF - intG) << '\n';

                mass_last_previous = mass_last; // set current last to previous last for next time slab
            }

            // #ifdef USE_MPI
            //          if ((iterations == 1) && MPIcf::IamMaster()) {
            // #else

            if ((iterations == 1)) {
                Fun_h sol_h(Wh, sol);
                Paraview<Mesh> writerTh(Th, path_output_figures + "Th.vtk");
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) + ".vtk");
                writer.add(b0h, "bulk", 0, 1);
                writer.add(sol_h, "bulk_end", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);

                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(g_Neumann, "neumann", 0, 1);
                writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);
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

            if (iterations > 1 && iter == total_number_iteration - 1)
                outputData << h << "," << dT << "," << errBulk << '\n';

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
