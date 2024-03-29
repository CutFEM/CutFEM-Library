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
 * @brief Time-dependent convection diffusion equation in Omega1 or Omega2.
 * @note Numerical method is the discontinuous Galerkin (dG) space-time CutFEM.

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

 * Numerical method:
 * A space-time discontinuous Cutfem, using the level-set method.

 * Classical scheme: Integration by parts on both diffusion and convection term.
 * Conservative scheme: Use Reynold's transport theorem instead. Slightly unclear
 * how to add a consistent numerical flux for the advection term.
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
#include "../algoim/quadrature_general.hpp"

using namespace globalVariable; // to access some globally defined constants

// Numerical examples
namespace Example1 {
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return -(sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17) - Epsilon;
    // return -sqrt((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.22)*(P[1]-0.22)) - 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return -(sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17) - Epsilon;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // Wrong sign of normal vector at interface
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
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17 - Epsilon;
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

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17 - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17);
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

// Normal x-direction
R n1(double *P, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return (P[0] - xc) / (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
}

// Normal y-direction
R n2(double *P, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return (P[1] - yc) / (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
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

    return (M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / 125 -
           (4 * M_PI * cos(M_PI * x) * cos(M_PI * y) * sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * x) * sin(M_PI * y) * (x - 0.5)) / 5 +
           (2 * M_PI * M_PI * cos(2 * M_PI * t) * cos(M_PI * y) * sin(M_PI * x) * (y - 0.5)) / 5;
}
} // namespace Example1_Omega1

namespace Lehrenfeld {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    double r0 = 1. + Epsilon;
    double x = P[0], y = P[1];

    return sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) - r0;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    double r0 = 1. + Epsilon;
    return sqrt(P[0] * P[0] + P[1] * P[1]) - r0;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
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
           2 * M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) +
           M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
               (8 * t * t * y * y + 4 * t * (x + t * (y * y - 1)) + 2) +
           M_PI * M_PI * cos(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * sin(M_PI * t) *
               (2 * x + 2 * t * (y * y - 1)) * (2 * x + 2 * t * (y * y - 1)) +
           M_PI * M_PI * cos(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * sin(M_PI * t) *
               (2 * y + 4 * t * y * (x + t * (y * y - 1))) * (2 * y + 4 * t * y * (x + t * (y * y - 1))) +
           M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * (y * y - 1) *
               (2 * x + 2 * t * (y * y - 1)) -
           2 * M_PI * sin(M_PI * t) * sin(M_PI * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
               (y * y - 1) * (x + t * (y * y - 1));
}

} // namespace Lehrenfeld

namespace Lehrenfeld_Convection_Dominated {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    double r0 = 1. + Epsilon;
    double x = P[0], y = P[1];

    return sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) - r0;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    double r0 = 1. + Epsilon;
    return sqrt(P[0] * P[0] + P[1] * P[1]) - r0;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

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

} // namespace Lehrenfeld_Convection_Dominated

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// The below parameters can be varied according to the options to use different methods,
// different numerical examples, different subdomains and boundary conditions etc.

//* Set numerical example (options: "example1", "lehrenfeld")
#define example1
//* Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominated

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
#if defined(omega1)
using namespace Example1_Omega1;
#elif defined(omega2)
using namespace Example1;
#endif
#endif
#if defined(lehrenfeld)
using namespace Lehrenfeld_Convection_Dominated;
#endif

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);
    // Mesh settings and data objects
    const size_t iterations = 9; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;             // starting mesh size
    double dT = 0.25;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#ifdef example1
    // Paths to store data
    const std::string path_output_data    = "../output_files/bulk/example1/data/";
    const std::string path_output_figures = "../output_files/bulk/example1/paraview/";
#elif defined(example2)
    // Paths to store data
    const std::string path_output_data    = "../output_files/bulk/example2/data/";
    const std::string path_output_figures = "../output_files/bulk/example2/paraview/";
#elif defined(lehrenfeld)
    // Paths to store data
    const std::string path_output_data    = "../output_files/bulk/lehrenfeld/data/";
    const std::string path_output_figures = "../output_files/bulk/lehrenfeld/paraview/";
#endif

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    // Data file to hold problem data
    std::ofstream outputData(path_output_data + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors;                  // array to hold bulk errors
    std::array<int, iterations> number_of_stabilized_edges; // array to count stabilized edges
    std::array<double, iterations> hs;                      // array to hold mesh sizes
    std::array<double, iterations> nxs;                     // array to hold mesh sizes
    std::array<double, iterations> nys;                     // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
#if defined(example1) || defined(example2)
        const double lx = 1., ly = 1.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0., 0., lx, ly);
#elif defined(lehrenfeld)
        const double lx = 7., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3.5, -1.5, lx, ly);
#endif

        // Parameters
        const double tfinal = .5; // Final time

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

        const double D = 0.01;

        double lambda = 1.; // Nitsche's method penalty parameter

        // DG interior penalty parameters
        const double lambda_A = 500, lambda_B = 0.5;

        // DG stabilization parameter
        const double tau0 = 2.5e-3, tau1 = 2.5e-3;

        FESpace2 Vh(Th, DataFE<Mesh>::P1dc); // discontinuous basis functions

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        FESpace2 Vh2(Th,
                     DataFE<Mesh>::P2);     // higher order space for interpolation
        FESpace2 Vh3(Th, DataFE<Mesh>::P3); // continuous basis functions

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        FESpace1 Ih2(Qh, DataFE<Mesh1>::P2Poly);
        // FESpace1 Ih3(Qh, DataFE<Mesh1>::P3Poly);

        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(3)); // specify order of quadrature in time
        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

// Velocity field
#if defined(example1) || defined(example2)
        Lagrange2 FEvelocity(1);
#elif defined(lehrenfeld)
        Lagrange2 FEvelocity(2);
#endif
        FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

#if defined(levelsetexact)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);
#elif defined(levelsetsolve)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet);
#endif

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double q0_0, q0_1, qp_0, qp_1;              // integral values to be computed
        double intF = 0, int_outflow = 0, intG = 0; // hold integrals of rhs and Neumann bcs

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

#if defined(levelsetexact)
                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);
#endif
                interface.init(i, Th, ls[i]);

#if defined(levelsetsolve)
                // We solve for the level-set using Crank-Nicholson in time
                if (i < lastQuadTime) {
                    LevelSet::move(ls[i], vel, vel, dt_levelSet, ls[i + 1]);
                }
#endif
            }

            // Create active meshes
            ActiveMesh<Mesh> Thi(Th);

            Thi.truncate(interface, -1); // remove part with negative sign of level

            // Cut FE space
            CutSpace Wh(Thi, Vh);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Vh, In, fun_rhsBulk);
            Fun_h g(Vh, In, fun_uBulk);                 // create an FE-function of the exact bulk
                                                        // solution Omega2
            Fun_h g_Neumann(Vh, In, fun_neumann_Gamma); // computer Neumann BC
            // Fun_h g_Neumann_left(Wh, In, fun_neumann_left);     // label 1
            // Fun_h g_Neumann_bottom(Wh, In, fun_neumann_bottom); // label 4
            // Fun_h g_Neumann_right(Wh, In, fun_neumann_right);   // label 2
            // Fun_h g_Neumann_top(Wh, In, fun_neumann_top);       // label 3

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);

            // Data for initial solution
            Rn data_u0; // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0))); // initial data bulk

            if (iter == 0)
                interpolate(Wh, data_B0, fun_uBulkInit);

            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);

            // Plot initial solution in paraview
            // #ifdef USE_MPI
            // if (iter == 0 && MPIcf::IamMaster()) {
            // #else
            if (iter == 0) {
                // #endif
                Paraview<Mesh> writerInitial(Thi, path_output_figures + "BulkInitialDG.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);

                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);
                Fun_h uRhs(Wh, fun_rhsBulk, 0.);

                writerInitial.add(uBex, "bulk_exact", 0, 1);
                writerInitial.add(uRhs, "bulk_rhs", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
            }

            //** Assembling linear and bilinear forms

            //* Time terms

#ifdef conservative

            convdiff.addBilinear(-innerProduct(u, dt(v)), Thi, In);

            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);

#else // classical scheme

            convdiff.addBilinear(+innerProduct(dt(u), v), Thi, In);

            // Time penalty term bulk LHS
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);

#endif

            //* Diffusion term
            convdiff.addBilinear(+innerProduct(D * grad(u), grad(v)), Thi, In);

            // Diffusion on inner edges (E_h)
            convdiff.addBilinear(-innerProduct(D * average(grad(u) * n), jump(v)) -
                                     innerProduct(D * jump(u), average(grad(v) * n)) +
                                     innerProduct(lambda_A / h * jump(u), jump(v)),
                                 Thi, INTEGRAL_INNER_EDGE_2D, In);

            //* Convection term
            convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Thi, In);

            // Convection on inner edges (E_h)
            convdiff.addBilinear(+innerProduct(average(vel * n * u), jump(v)) +
                                     innerProduct(0.5 * fabs(vel * n) * jump(u), jump(v)),
                                 Thi, INTEGRAL_INNER_EDGE_2D, In);

            //* Stabilization

#if defined(fullstab)

            convdiff.addFaceStabilization(+innerProduct(1. / h * tau0 * jump(u), jump(v)) +
                                              innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)),
                                          Thi, In);

#elif defined(macro)

            // TimeMacroElement<Mesh> TimeMacro(Thi, qTime, 0.16);
            MacroElementPartition<Mesh> TimeMacro(Thi, 0.45);

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
            }

            // Stabilization of the bulk
            // convdiff.mat_.clear();
            convdiff.addFaceStabilization(+innerProduct(1. / h * tau0 * jump(u), jump(v)) +
                                              innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)),
                                          Thi, In, TimeMacro);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // getchar();

#endif

            number_of_stabilized_edges.at(j) = convdiff.get_number_of_stabilized_edges();

            //* Boundary conditions

// If Omega1
#ifdef omega1

//* Dirichlet on both outer and inner boundary
#if defined(dirichlet1) && defined(dirichlet2)
            // Dirichlet outer

            convdiff.addBilinear(+innerProduct(u, 0.5 * fabs(vel * n) * v) + innerProduct(u, 0.5 * (vel * n) * v) -
                                     innerProduct(D * grad(u) * n, v)   // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v), // added penalty
                                 Thi, INTEGRAL_BOUNDARY, In);

            convdiff.addLinear(+innerProduct(g.expr(), 0.5 * fabs(vel * n) * v) -
                                   innerProduct(g.expr(), 0.5 * (vel * n) * v) -
                                   innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                               Thi, INTEGRAL_BOUNDARY, In);

            convdiff.addBilinear(+innerProduct(u, 0.5 * fabs(vel * n) * v) + innerProduct(u, 0.5 * (vel * n) * v) -
                                     innerProduct(D * grad(u) * n, v)   // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v)  // added penalty
                                 ,
                                 interface, In);

            convdiff.addLinear(+innerProduct(g.expr(), 0.5 * fabs(vel * n) * v) -
                                   innerProduct(g.expr(), 0.5 * (vel * n) * v) -
                                   innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
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

            convdiff.addBilinear(innerProduct(vel.exprList() * u * n, v), Thi, INTEGRAL_BOUNDARY, In);

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

#if defined(classical) && defined(omega2) && defined(dirichlet)
            convdiff.addBilinear(+innerProduct((vel * n) * u, 0.5 * v) + innerProduct(0.5 * fabs(vel * n) * u, v),
                                 interface, In);

            convdiff.addLinear(-innerProduct(g.expr(), (vel * n) * (0.5 * v)) +
                                   innerProduct(g.expr(), fabs(vel * n) * 0.5 * v),
                               interface, In);
#elif defined(classical) && defined(omega2) && defined(neumann)
            convdiff.addBilinear(+innerProduct((vel * n) * u, v), interface, In);
#endif

            // Add RHS in bulk
            convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);
#ifndef USE_MPI
            if (iter == total_number_iteration - 1)
                matlab::Export(convdiff.mat_[0], path_output_data + "mat_h_dg" + std::to_string(h) + "_" +
                                                     std::to_string(j + 1) + ".dat");
#elif defined(USE_MPI)

            if ((iter == total_number_iteration - 1) && MPIcf::IamMaster()) {
                matlab::Export(convdiff.mat_[0], path_output_data + "mat_dg_rank_" + std::to_string(MPIcf::my_rank()) +
                                                     "_h" + std::to_string(h) + "_" + std::to_string(j + 1) + ".dat");
            }

#endif

            // Solve linear system
            convdiff.solve("mumps");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Compute conservation error
            if (iterations == 1) {
                Fun_h funuh(Wh, data_u0);

                intF = integral(Thi, In, f, 0, qTime);
#if defined(conservative) && defined(omega1) && defined(neumann1) && defined(neumann2)
                auto outflow = (vel * n) * funuh.expr();
                int_outflow  = integral(Thi, In, (vel * n) * b0h.expr(), INTEGRAL_BOUNDARY,
                                        qTime); // integral(outflow, In, interface, 0);
                intG         = -integral(g_Neumann, In, interface, 0);
#endif

#if defined(omega2) && defined(neumann)
                intG = integral(g_Neumann, In, interface, 0);

#endif

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_u0(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Thi, funsol, 0, 0);
                sol2 += data_u0(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Thi, funsol, 0, lastQuadTime);

                if (iter == 0) {
                    q0_0 = q_0;
                    q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Thi, b0h, 0, 0);
                }

                outputData << std::setprecision(10);
                outputData << current_time << "," << (q_1 - qp_1) << "," << intF << "," << intG << ","
#if (defined(omega2) && dirichlet) || (defined(omega1) && defined(dirichlet1))
                    ;
#elif defined(omega2) && defined(neumann)
                           << ((q_1 - qp_1) - intF - intG) << '\n';
#elif defined(omega1) && defined(neumann1)
                           << ((q_1 - qp_1) - intF - intG + int_outflow) << '\n';
#endif
                qp_1 = q_1;
            }

            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));
            Fun_h funuh(Wh, sol);
            // double errBulk = L2normCut(funuh, fun_uBulkD, current_time, 0, 1);
            // std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << '\n';

            sol += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            errors.at(j) = errBulk;

            // #ifdef USE_MPI
            //          if ((iterations == 1) && MPIcf::IamMaster()) {
            // #else
            if ((iterations == 1)) {
                Fun_h sol(Wh, data_u0);
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_dg" + std::to_string(iter + 1) + ".vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                Expression2 uuh(sol, 0, op_id);
                Expression2 uuex(uBex, 0, op_id);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                // writer.add(fabs(uuh - uuex), "bulk_error");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);
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
        // h *= 0.5;
        h *= sqrt(0.5);
#endif
    }

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

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    std::cout << '\n';

    std::cout << "number of stabilized edges = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << number_of_stabilized_edges.at(i);
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
