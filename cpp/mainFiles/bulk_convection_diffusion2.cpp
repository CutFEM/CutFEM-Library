/*
   We consider a time-dependent bulk problem on Omega2.
   We consider problems with both Neumann and Dirichlet boundary conditions,
   and we consider scheme II, III and the Reynold scheme.

   // Problem:
   Find u in Omega2 such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    on Omega2.

    // Numerical method:
    A space-time Cutfem, using the level-set method,
    which allows for both dg and cg.

    Classical : Integration by parts on full convection term, no term added to make anti-symmetric
    Conservative: Reynold's transport theorem is used to make the bilinear form fulfill a conservation law
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
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "projection.hpp"
// #include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"

using namespace globalVariable;

namespace Preuss {
const double R0 = .5;
const double D  = 1.;
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

double rho(double *P, const double t) { return 1. / pi * std::sin(2 * pi * t); }

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    return -(std::sqrt(P[0] * P[0] + (P[1] - rho(P, t)) * (P[1] - rho(P, t))) - R0) - Epsilon;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {

    return -(std::sqrt(P[0] * P[0] + (P[1] - rho(P, 0)) * (P[1] - rho(P, 0))) - R0) - Epsilon;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) { return 0.; }

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    if (i == 0)
        return 0.;
    else
        return 2 * std::cos(2 * pi * t);
}

// // Normal x-direction
// R n1(double *P, const R t) {
//    R xc = 0.5 + 0.28 * std::sin(M_PI * t), yc = 0.5 - 0.28 * std::cos(M_PI * t);
//    return (P[0] - xc) /
//           (std::sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// // Normal y-direction
// R n2(double *P, const R t) {
//    R xc = 0.5 + 0.28 * std::sin(M_PI * t), yc = 0.5 - 0.28 * std::cos(M_PI * t);
//    return (P[1] - yc) /
//           (std::sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) {
    return std::pow(std::cos(pi * std::sqrt(P[0] * P[0] + P[1] * P[1]) / (2 * R0)), 2);
}

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * std::sin(2 * pi * t);
    double r   = std::sqrt(x * x + (y - rho) * (y - rho));
    return std::pow(std::cos(pi * r / (2 * R0)), 2);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * std::sin(2 * pi * t);
    double r   = std::sqrt(x * x + (y - rho) * (y - rho));
    return std::pow(std::cos(pi * r / (2 * R0)), 2);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // with eps
    return -D * ((1.0 / (R0 * R0) * (x * x) * (3.141592653589793 * 3.141592653589793) *
                  pow(cos((3.141592653589793 *
                           sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x +
                                1.0E-16)) /
                          (R0 * 2.0)),
                      2.0) *
                  (-1.0 / 2.0)) /
                     (pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16) +
                 (1.0 / (R0 * R0) * (x * x) * (3.141592653589793 * 3.141592653589793) *
                  pow(sin((3.141592653589793 *
                           sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x +
                                1.0E-16)) /
                          (R0 * 2.0)),
                      2.0)) /
                     (pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) * 2.0 + (x * x) * 2.0 +
                      2.0E-16) -
                 (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) *
                  pow(cos((3.141592653589793 *
                           sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x +
                                1.0E-16)) /
                          (R0 * 2.0)),
                      2.0) *
                  pow(y * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1, 2.0)) /
                     (pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) * 8.0 + (x * x) * 8.0 +
                      8.0E-16) +
                 (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) *
                  pow(sin((3.141592653589793 *
                           sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x +
                                1.0E-16)) /
                          (R0 * 2.0)),
                      2.0) *
                  pow(y * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1, 2.0)) /
                     (pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) * 8.0 + (x * x) * 8.0 +
                      8.0E-16) -
                 (3.141592653589793 *
                  cos((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  sin((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  1.0 / sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16) *
                  2.0) /
                     R0 +
                 ((x * x) * 3.141592653589793 *
                  cos((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  sin((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  1.0 /
                  pow(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16,
                      3.0 / 2.0)) /
                     R0 +
                 (3.141592653589793 *
                  cos((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  sin((3.141592653589793 *
                       sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                      (R0 * 2.0)) *
                  pow(y * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1, 2.0) * 1.0 /
                  pow(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16,
                      3.0 / 2.0)) /
                     (R0 * 4.0)) -
           (3.141592653589793 *
            cos((3.141592653589793 *
                 sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                (R0 * 2.0)) *
            sin((3.141592653589793 *
                 sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                (R0 * 2.0)) *
            cos(t * 3.141592653589793 * 2.0) * (y * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1) *
            1.0 / sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
               R0 +
           ((3.141592653589793 * 3.141592653589793) *
            cos((3.141592653589793 *
                 sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                (R0 * 2.0)) *
            sin((3.141592653589793 *
                 sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16)) /
                (R0 * 2.0)) *
            cos(t * 3.141592653589793 * 2.0) * (y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1) * 1.0 /
            sqrt(pow(y - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x + 1.0E-16) *
            6.366197723675814E-1) /
               R0;

    // without eps
    return (D * 1.0 / (R0 * R0) * 3.141592653589793 * 1.0 /
            std::pow(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x,
                     3.0 / 2.0) *
            (3.141592653589793 *
                 std::pow(std::cos((3.141592653589793 *
                                    std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1,
                                                       2.0) +
                                              x * x)) /
                                   (R0 * 2.0)),
                          2.0) *
                 std::pow(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x,
                          3.0 / 2.0) *
                 1.622592768292134E+32 -
             3.141592653589793 *
                 std::pow(std::sin((3.141592653589793 *
                                    std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1,
                                                       2.0) +
                                              x * x)) /
                                   (R0 * 2.0)),
                          2.0) *
                 std::pow(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x,
                          3.0 / 2.0) *
                 1.622592768292134E+32 +
             R0 * (x * x) *
                 std::sin((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 std::cos((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 3.245185536584267E+32 +
             R0 * (y * y) *
                 std::sin((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 std::cos((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 3.245185536584267E+32 +
             R0 *
                 std::sin((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 std::pow(std::sin(t * 3.141592653589793 * 2.0), 2.0) *
                 std::cos((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 3.28806039705713E+31 -
             R0 * y *
                 std::sin((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 std::sin(t * 3.141592653589793 * 2.0) *
                 std::cos((3.141592653589793 *
                           std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) +
                                     x * x)) /
                          (R0 * 2.0)) *
                 2.065949277590844E+32)) /
               3.245185536584267E+32 -
           (3.141592653589793 *
            std::sin(
                (3.141592653589793 *
                 std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x)) /
                (R0 * 2.0)) *
            std::cos(t * 3.141592653589793 * 2.0) *
            std::cos(
                (3.141592653589793 *
                 std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x)) /
                (R0 * 2.0)) *
            (y * 2.0 - std::sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1) * 1.0 /
            std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x)) /
               R0 +
           ((3.141592653589793 * 3.141592653589793) *
            std::sin(
                (3.141592653589793 *
                 std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x)) /
                (R0 * 2.0)) *
            std::cos(t * 3.141592653589793 * 2.0) *
            std::cos(
                (3.141592653589793 *
                 std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x)) /
                (R0 * 2.0)) *
            (y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1) * 1.0 /
            std::sqrt(std::pow(y - std::sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + x * x) *
            6.366197723675814E-1) /
               R0;
}

R fun_neumann_left(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

R fun_neumann_bottom(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

R fun_neumann_right(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

R fun_neumann_top(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

} // namespace Preuss

namespace Example1 {
/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * std::sin(M_PI * t), yc = 0.5 - 0.28 * std::cos(M_PI * t);
    // return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17);
    return -(std::sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17) - Epsilon; //! PUT BACK
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    // return ((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17);
    return -(std::sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17) - Epsilon; //! PUT BACK
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // Correct sign of normal vector at interface
    return -(pi * std::cos(2 * pi * t) * std::cos(pi * y) * std::sin(pi * x) * (x - 1. / 2 - 0.28 * std::sin(pi * t))) /
               (250 * std::sqrt((std::pow(x - 1. / 2 - 0.28 * std::sin(pi * t), 2) +
                                 std::pow(y - 0.5 + 0.28 * std::cos(pi * t), 2)))) -
           (pi * std::cos(2 * pi * t) * std::cos(pi * x) * std::sin(pi * y) * (y - 0.5 + 0.28 * std::cos(pi * t))) /
               (250 * std::sqrt((std::pow((x - 1. / 2 - 0.28 * std::sin(pi * t)), 2) +
                                 std::pow(y - 0.5 + 0.28 * std::cos(pi * t), 2))));
}

// Velocity field
R fun_velocity(double *P, const int i, const R t) {
    if (i == 0)
        return M_PI * (0.5 - P[1]);
    else
        return M_PI * (P[0] - 0.5);
}

// // Normal x-direction
// R n1(double *P, const R t) {
//    R xc = 0.5 + 0.28 * std::sin(M_PI * t), yc = 0.5 - 0.28 * std::cos(M_PI * t);
//    return (P[0] - xc) /
//           (std::sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// // Normal y-direction
// R n2(double *P, const R t) {
//    R xc = 0.5 + 0.28 * std::sin(M_PI * t), yc = 0.5 - 0.28 * std::cos(M_PI * t);
//    return (P[1] - yc) /
//           (std::sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.5 + 0.4 * std::cos(M_PI * P[0]) * std::cos(M_PI * P[1]); }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.5 + 0.4 * std::cos(M_PI * P[0]) * std::cos(M_PI * P[1]) * std::cos(2 * M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.5 + 0.4 * std::cos(M_PI * P[0]) * std::cos(M_PI * P[1]) * std::cos(2 * M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (M_PI * M_PI * std::cos(2 * M_PI * t) * std::cos(M_PI * x) * std::cos(M_PI * y)) / 125 -
           (4 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y) * std::sin(2 * M_PI * t)) / 5 -
           (2 * M_PI * M_PI * std::cos(2 * M_PI * t) * std::cos(M_PI * x) * std::sin(M_PI * y) * (x - 0.5)) / 5 +
           (2 * M_PI * M_PI * std::cos(2 * M_PI * t) * std::cos(M_PI * y) * std::sin(M_PI * x) * (y - 0.5)) / 5;
}

R fun_neumann_left(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * std::cos(2 * pi * t) * std::cos(pi * y) * std::sin(pi * x)) / 250;
}

R fun_neumann_bottom(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return (pi * std::cos(2 * pi * t) * std::cos(pi * x) * std::sin(pi * y)) / 250;
}

R fun_neumann_right(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * std::cos(2 * pi * t) * std::cos(pi * y) * std::sin(pi * x)) / 250;
}

R fun_neumann_top(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return -(pi * std::cos(2 * pi * t) * std::cos(pi * x) * std::sin(pi * y)) / 250;
}

} // namespace Example1

namespace Lehrenfeld {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    double r0 = 1.;
    double x = P[0], y = P[1];

    return -(std::sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) - r0);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    double r0 = 1.;
    return -(std::sqrt(P[0] * P[0] + P[1] * P[1]) - r0);
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const R t) {
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
    return std::cos(M_PI * std::sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) / r0) * std::sin(M_PI * t);
    // return std::cos(M_PI*((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/(r0*r0))*std::sin(M_PI*t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double r0 = 1., x = P[0], y = P[1];
    return std::cos(M_PI * std::sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) / r0) * std::sin(M_PI * t);
    // return std::cos(M_PI*((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/(r0*r0))*std::sin(M_PI*t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return M_PI * std::cos(M_PI * t) *
               std::cos(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) +
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y))) /
               std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y) +
           (M_PI * M_PI * std::sin(M_PI * t) *
            std::cos(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (2 * x + 2 * t * (y * y - 1)) * (2 * x + 2 * t * (y * y - 1))) /
               (4 * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) +
           (M_PI * M_PI * std::sin(M_PI * t) *
            std::cos(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (2 * y + 4 * t * y * (x + t * (y * y - 1))) * (2 * y + 4 * t * y * (x + t * (y * y - 1)))) /
               (4 * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) +
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (8 * t * t * y * y + 4 * t * (x + t * (y * y - 1)) + 2)) /
               (2 * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) -
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (2 * x + 2 * t * (y * y - 1)) * (2 * x + 2 * t * (y * y - 1))) /
               (4 * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y) *
                std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) -
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) *
            (2 * y + 4 * t * y * (x + t * (y * y - 1))) * (2 * y + 4 * t * y * (x + t * (y * y - 1)))) /
               (4 * ((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y) *
                std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) -
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * (y * y - 1) *
            (x + t * (y * y - 1))) /
               std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y) +
           (M_PI * std::sin(M_PI * t) *
            std::sin(M_PI * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y)) * (y * y - 1) *
            (2 * x + 2 * t * (y * y - 1))) /
               (2 * std::sqrt((x + t * (y * y - 1)) * (x + t * (y * y - 1)) + y * y));

    // return M_PI*std::cos(M_PI*t)*std::cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)) +
    // 2*M_PI*std::sin(M_PI*t)*std::sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)) +
    // M_PI*std::sin(M_PI*t)*std::sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(8*t*t*y*y + 4*t*(x + t*(y*y -
    // 1)) + 2) + M_PI*M_PI*std::cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*std::sin(M_PI*t)*(2*x + 2*t*(y*y
    // - 1))*(2*x + 2*t*(y*y - 1)) + M_PI*M_PI*std::cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y
    // - 1)) + y*y))*std::sin(M_PI*t)*(2*y + 4*t*y*(x + t*(y*y - 1)))*(2*y + 4*t*y*(x + t*(y*y - 1))) +
    // M_PI*std::sin(M_PI*t)*std::sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(2*x + 2*t*(y*y - 1))
    // - 2*M_PI*std::sin(M_PI*t)*std::sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(x + t*(y*y - 1));
}

} // namespace Lehrenfeld

namespace Lehrenfeld_6_2 {

// Level-set function
R fun_levelSet(double *P, const int i, const R t) {
    double r0 = .5 + Epsilon;
    double x = P[0], y = P[1];

    // double rho(R x, R y, R t) return 1./M_PI*std::sin(2*M_PI*t);
    // double r(R x, R y, R t) return std::sqrt((x-1./M_PI*std::sin(2*M_PI*t))*(x-1./M_PI*std::sin(2*M_PI*t)) + y*y);

    return -(std::sqrt((x - 1. / M_PI * std::sin(2 * M_PI * t)) * (x - 1. / M_PI * std::sin(2 * M_PI * t)) + y * y) -
             r0);
}

// Level-set function initial
R fun_levelSet(double *P, const int i) {
    double r0 = .5 + Epsilon;
    double x = P[0], y = P[1];

    return -(std::sqrt(x * x + y * y) - r0);
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const R t) {
    if (i == 0)
        return 2 * std::cos(2 * M_PI * t);
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }
} // namespace Lehrenfeld_6_2

namespace Lehrenfeld_6_2_Convection_Dominated {

// Level-set function
R fun_levelSet(double *P, const int i, const R t) {
    double r0 = .5 + Epsilon;
    double x = P[0], y = P[1];

    // double rho(R x, R y, R t) return 1./M_PI*std::sin(2*M_PI*t);
    // double r(R x, R y, R t) return std::sqrt((x-1./M_PI*std::sin(2*M_PI*t))*(x-1./M_PI*std::sin(2*M_PI*t)) + y*y);

    return -(std::sqrt((x - 1. / M_PI * std::sin(2 * M_PI * t)) * (x - 1. / M_PI * std::sin(2 * M_PI * t)) + y * y) -
             r0);
}

// Level-set function initial
R fun_levelSet(double *P, const int i) {
    double r0 = .5 + Epsilon;
    double x = P[0], y = P[1];

    return -(std::sqrt(x * x + y * y) - r0);
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const R t) {
    if (i == 0)
        return 2 * std::cos(2 * M_PI * t);
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double r0 = .5, x = P[0], y = P[1];

    return std::cos(M_PI *
                    ((x - 1. / M_PI * std::sin(2 * M_PI * t)) * (x - 1. / M_PI * std::sin(2 * M_PI * t)) + y * y) /
                    (r0 * r0)) *
           std::sin(M_PI * t);
    // return std::cos(M_PI*std::sqrt((x-1./M_PI*std::sin(2*M_PI*t))*(x-1./M_PI*std::sin(2*M_PI*t)) +
    // y*y)/r0)*std::sin(M_PI*t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double r0 = .5, x = P[0], y = P[1];

    return std::cos(M_PI *
                    ((x - 1. / M_PI * std::sin(2 * M_PI * t)) * (x - 1. / M_PI * std::sin(2 * M_PI * t)) + y * y) /
                    (r0 * r0)) *
           std::sin(M_PI * t);
    // return std::cos(M_PI*std::sqrt((x-1./M_PI*std::sin(2*M_PI*t))*(x-1./M_PI*std::sin(2*M_PI*t)) +
    // y*y)/r0)*std::sin(M_PI*t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];
    return 3.1415927 * std::cos(3.1415927 * t) *
               std::cos(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) +
           0.50265482 * std::sin(3.1415927 * t) *
               std::sin(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) +
           6.3165468 * y * y * std::sin(3.1415927 * t) *
               std::cos(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) +
           1.5791367 * std::sin(3.1415927 * t) *
               std::cos(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) *
               (2.0 * x - 0.63661977 * std::sin(6.2831853 * t)) * (2.0 * x - 0.63661977 * std::sin(6.2831853 * t)) -
           25.132741 * std::cos(6.2831853 * t) * std::sin(3.1415927 * t) *
               std::sin(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) *
               (2.0 * x - 0.63661977 * std::sin(6.2831853 * t)) +
           50.265482 * std::cos(6.2831853 * t) * std::sin(3.1415927 * t) *
               std::sin(12.566371 * (x - 0.31830989 * std::sin(6.2831853 * t)) *
                            (x - 0.31830989 * std::sin(6.2831853 * t)) +
                        12.566371 * y * y) *
               (x - 0.31830989 * std::sin(6.2831853 * t));
}

} // namespace Lehrenfeld_6_2_Convection_Dominated

namespace Lehrenfeld_6_3 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    return -(std::min(std::sqrt(x * x + (y - t + 3. / 4) * (y - t + 3. / 4)),
                      std::sqrt(x * x + (y - t + 3. / 4) * (y + t - 3. / 4))) -
             0.5);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {

    double x = P[0], y = P[1];

    return -(std::min(std::sqrt(x * x + (y + 3. / 4) * (y + 3. / 4)), std::sqrt(x * x + (y + 3. / 4) * (y - 3. / 4))) -
             0.5);
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const R t) {

    if ((P[1] > 0 && t <= 3. / 4) || (P[1] < 0 && t > 3. / 4)) {
        if (i == 0) {
            return 0.;
        } else
            return -1.;
    }

    else if ((P[1] <= 0 && t <= 3. / 4) || (P[1] > 0 && t > 3. / 4)) {
        if (i == 0) {
            return 0.;
        } else
            return 1.;
    }

    else {
        std::cout << "something wrong\n";
        return 0.;
    }
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) {

    // Sign(y)
    if (P[1] > 0)
        return 1.;
    else if (P[1] == 0)
        return 0.;
    else
        return -1.;
}

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) { return 0.; }

R fun_uBulkD(double *P, const int i, const int d, const R t) { return 0.; }

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

} // namespace Lehrenfeld_6_3

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

// Note: standard for this program is to solve the equations on the inner domain Omega 2.

// Choose Discontinuous or Continuous Galerkin method (options: "dg", "cg")
#define cg
// Set numerical example (options: "lehrenfeld_6_2", "lehrenfeld_6_3")
#define preuss

#define convection_dominatednot
// Set boundary condition type on Omega 2 (options: "dirichlet", "neumann" – note: neumann only works combined with
// example1)
#define neumann
// Set scheme for the dg method (options: "conservative", "classical" see thesis. Irrelevant if "cg" is defined instead
// of "dg")
#define classical
// Set stabilization method (options: "fullstab", "macro")
#define fullstab
// Decide whether to solve for level set function, or to use exact (options: "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h
#define use_tnot

#if defined(example1)
using namespace Example1;

#elif defined(lehrenfeld_6_2)
#ifdef convection_dominated
using namespace Lehrenfeld_6_2_Convection_Dominated;
#else
using namespace Lehrenfeld_6_2;
#endif
#elif defined(lehrenfeld_6_3)
using namespace Lehrenfeld_6_3;
#elif defined(preuss)
using namespace Preuss;
#endif

// dT = 0.125
// [0.182859, 0.138091, 0.114025] fullstab
// [0.155091, 0.112177, 0.024969] macro

// dT = 0.0625
// [0.138091, 0.114025, 0.108953] fullstab
// [0.112177, 0.024969, 0.0234635] macro

int main(int argc, char **argv) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Mesh settings and data objects
    const size_t iterations = 4; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 5, ny = 5;          // starting mesh size
    double h  = 0.2;             // starting mesh size
    double dT = 0.125;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#if defined(lehrenfeld_6_2) || defined(example1) || defined(preuss)
    // Paths to store data
    const std::string path_output  = "/NOBACKUP/smyrback/output_files/paper/preuss/data/";
    const std::string path_figures = "/NOBACKUP/smyrback/output_files/paper/preuss/paraview/";
#endif

    // Initialize MPI
    // MPIcf cfMPI(argc, argv);

    // Create directory if not already existent
    // if (MPIcf::IamMaster()) {
    //     std::experimental::filesystem::create_directories(path_output);
    //     std::experimental::filesystem::create_directories(path_figures);
    // }

    // Data file to hold problem data

    std::ofstream output_data(path_output + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors_bulk;        // array to hold bulk errors
    std::array<int, iterations> numb_stabilized_edges; // array to count stabilized edges
    std::array<double, iterations> reynold_error;
    std::array<double, iterations> hs; // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {
// Mesh size
// double h = std::pow(0.5, j+1);   //0.9*std::pow(0.5, j);    //lx/(nx);
// double h = std::sqrt(lx*lx/(nx*nx) + ly*ly/(ny*ny));

// Time
// double dT = std::pow(2, -j-2);    // Time step size

// Define background mesh
#if defined(example1) || defined(example2)
        const double lx = 1., ly = 1.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0., 0., lx, ly);
#elif defined(lehrenfeld_6_2)
        const double lx = 2., ly = 1.2;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -1., -0.6, lx, ly);
#elif defined(preuss)
        const double lx = 5., ly = 6.2;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -2.5, -3.1, lx, ly);
#endif

        // Parameters
        double tfinal = .1; // Final time

#ifdef use_t
        GTime::total_number_iteration = int(tfinal / dT);
#else
        int divisionMeshSize = 3;
        // int divisionMeshSize = 2*3*pi;
        // int divisionMeshSize = 18;

        double dT              = h / divisionMeshSize;
        // double dT = 3*h;
        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "Iteration " << j + 1 << "/" << iterations << std::endl;
        }

        std::cout << "h  = " << h << std::endl;
        std::cout << "nx = " << nx << std::endl;
        std::cout << "ny = " << ny << std::endl;
        std::cout << "dT = " << dT << std::endl;

#ifdef convection_dominated
        double A2 = 0.01;
#else
        double A2 = 1;
#endif
        double kappaTilde2 = 1;

        // Constants for penalty terms
        double tau_a2 = 500; // diffusion penalty scaling
        double tau_b2 = 10;  // convection penalty scaling

        // Bulk penalties
        double lambdaA = tau_a2 / h; // diffusion term
        double lambdaB = tau_b2;     // convection term

        // Penalty parameter for outer boundary
        double lambda = 500 / h; // only used when Dirichlet BCs apply

#ifdef dg
        // DG stabilization parameters
        double tau20 = 1e-1, tau21 = 1e-1; // bulk
        // DG Space
        FESpace2 Vh(Th, DataFE<Mesh>::P1dc); // discontinuous basis functions
#elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 1e-0;
        FESpace2 Vh(Th, DataFE<Mesh>::P1); // continuous basis functions
#endif

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        FESpace2 Vh2(Th,
                     DataFE<Mesh>::P2); // higher order space for interpolation

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
        Lagrange2 FEvelocity(2);
        FESpace2 VelVh(Th, FEvelocity);
        std::vector<Fun_h> vel(nbTime);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);

        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << std::endl;

        int iter = 0;
        double q0_0, q0_1, qp_0, qp_1;
        double intF = 0, intG = 0; // hold integrals of rhs and Neumann bcs

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * time_step;
            const TimeSlab &In(Ih[iter]);

            std::cout << " ------------------------------------------------------" << std::endl;
            std::cout << " ------------------------------------------------------" << std::endl;
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << std::endl;
            std::cout << " Time      \t : \t" << current_time << std::endl;
            std::cout << "dT = " << dT << std::endl;

            swap(ls[0], ls[lastQuadTime]);

            // computation of the interface in each of the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);

                vel[i].init(VelVh, fun_velocity, tt);

                interface.init(i, Th, ls[i]);
            }

            // Create active meshes
            ActiveMesh<Mesh> Kh2(Th);
            Kh2.truncate(interface, -1);

            // Cut FE space
            CutSpace Wh(Kh2, Vh);

            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Vh2, In, fun_rhsBulk);
            Fun_h g(Vh2, In, fun_uBulk); // FE-function of the exact bulk solution

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);

            // Data for initial solution
            Rn data_u0(convdiff.get_nb_dof(), 0.); // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0))); // initial data bulk

            if (iter == 0)
                interpolate(Wh, data_B0, fun_uBulkInit);

            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);

            // Plot initial solution in paraview
            if (iter == 0) {
                Paraview<Mesh> writerInitial(Kh2, path_figures + "BulkInitial.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);

                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);
                Fun_h uRhs(Wh, fun_rhsBulk, 0.);

                writerInitial.add(uBex, "bulk_exact", 0, 1);
                writerInitial.add(uRhs, "bulk_rhs", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
            }

            // * Assembling linear and bilinear forms

            // Partial derivative in time, bulk

// conservative scheme
#ifdef conservative

            convdiff.addBilinear(-innerProduct(u, dt(v)), Kh2, In);

            convdiff.addBilinear(+innerProduct(u, v), Kh2, (int)lastQuadTime, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Kh2, 0, In);

// classical scheme
#else

            convdiff.addBilinear(+innerProduct(dt(u), v), Kh2, In);
            convdiff.addBilinear(+innerProduct(u, v), Kh2, 0, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Kh2, 0, In);

#endif
            // Scheme for diffusion

            convdiff.addBilinear(+innerProduct(A2 * grad(u), grad(v)), Kh2, In);

#ifdef dg

            // Diffusion on inner edges (E_h)
            convdiff.addBilinear(-innerProduct(A2 * average(grad(u) * n), jump(v)) -
                                     innerProduct(A2 * jump(u), average(grad(v) * n)) +
                                     innerProduct(lambdaA * jump(u), jump(v)),
                                 Kh2, INTEGRAL_INNER_EDGE_2D, In);

            // Convection term
            convdiff.addBilinear(-innerProduct(u, (vel.expression() * grad(v))), Kh2, In);

            // Added terms
            convdiff.addBilinear(+innerProduct(average(vel * n * u), jump(v))
                                     //+ innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))
                                     + innerProduct(0.5 * fabs(vel * n) * jump(u), jump(v)),
                                 Kh2, INTEGRAL_INNER_EDGE_2D, In);

#elif defined(cg) && defined(classical)    // classic CG scheme

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(+innerProduct((vel[i].exprList() * grad(u)), v), Kh2, In, i);
            }
#elif defined(cg) && defined(conservative) // conservative CG scheme
            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Kh2, In, i);
                // convdiff.addBilinear(-innerProduct(u, (vel[i].expression() * grad(v))), Kh2, In, i);
            }
#endif

            //* Stabilization

#if defined(macro) || defined(macro2)

#ifdef macro
            double gamma = 0.125;
            TimeMacroElement<Mesh> TimeMacro(Kh2, qTime, gamma);
#endif

            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, path_figures + "Th" + std::to_string(iter + 1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                // writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

                // domain = 0,

                writerMacro.writeFaceStab(Kh2, 0,
                                          path_figures + "FullStabilization" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeActiveMesh(Kh2, path_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0, path_figures + "macro" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroInnerEdge(TimeMacro, 0,
                                                path_figures + "macro_inner_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroOutterEdge(TimeMacro, 0,
                                                 path_figures + "macro_outer_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeSmallElements(TimeMacro, 0,
                                               path_figures + "small_element" + std::to_string(iter + 1) + ".vtk");
            }

            // * Stabilization of the bulk

            // convdiff.mat_.clear();
            convdiff.addFaceStabilization(+innerProduct(1. / h * tau20 * jump(u), jump(v)) +
                                              innerProduct(h * tau21 * jump(grad(u) * n), jump(grad(v) * n)),
                                          Kh2, In, TimeMacro);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // getchar();

#elif defined(fullstab)

            // number_of_stabilized_edges.at(j) = convdiff.num_stabilized_edges(Th,
            // In);

            convdiff.addFaceStabilization(+innerProduct(1. / h * tau20 * jump(u), jump(v)) +
                                              innerProduct(h * tau21 * jump(grad(u) * n), jump(grad(v) * n))
                                          //+ innerProduct(tau20*(1+dT/h)/h/h*jump(u), jump(v))
                                          //+ innerProduct(tau21*h*h*jump(grad(u)), jump(grad(v)))
                                          ,
                                          Kh2, In);

#endif

            // numb_stabilized_edges.at(j) = convdiff.number_of_stabilized_edges;

            // Boundary conditions on interface

#ifdef neumann

            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
            convdiff.addLinear(+innerProduct(g_Neumann.expr(), v), interface, In);

#if defined(classical) && defined(dg)

            // convdiff.mat_.clear();

            convdiff.addBilinear(+innerProduct((vel * n) * u, v), interface, In);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // return 0;

#endif

#elif defined(dirichlet)
            convdiff.addBilinear(-innerProduct(A2 * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, A2 * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda * v)       // added penalty
                                 ,
                                 interface, In);

            // Diffusion inflow
            convdiff.addLinear(-innerProduct(g.expression(), A2 * grad(v) * n) +
                                   innerProduct(g.expression(), lambda * v),
                               interface, In);

            // Inflow and outflow terms on the boundary
            convdiff.addBilinear(+innerProduct(0.5 * (vel * n) * u, kappaTilde2 * v) +
                                     innerProduct(0.5 * fabs(vel * n) * u, kappaTilde2 * v),
                                 interface, In);

            convdiff.addLinear(-innerProduct(g.expression(), 0.5 * (vel * n) * v) +
                                   innerProduct(g.expression(), 0.5 * fabs(vel * n) * v),
                               interface, In);

#if defined(conservative) || defined(cg) // cg reynold, dg reynold and cg classic

            convdiff.addBilinear(-innerProduct((vel * n) * u, v) // from Reynold's transport theorem
                                 ,
                                 interface, In);
#endif

#endif

// Boundary conditions on outer boundary

// FIXME: Program crashes when using below
#if defined(omega1) && defined(example1)

#ifdef dirichlet1
            // solve with Dirichlet BCs

            convdiff.addBilinear(-innerProduct(A2 * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, A2 * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda * v)       // added penalty
                                 ,
                                 Kh2, INTEGRAL_BOUNDARY, In);

            // RHS on external boundary
#if defined(classical) && defined(dg)
            // NOTE: This gives errors
            convdiff.addLinear(-innerProduct(g.expression(), (vel * n) * v) -
                                   innerProduct(g.expression(), A2 * grad(v) * n) +
                                   innerProduct(g.expression(), lambda * v),
                               Kh2, INTEGRAL_BOUNDARY, In);

#elif defined(conservative) || defined(cg)
            convdiff.addLinear(-innerProduct(g.expression(), A2 * grad(v) * n) +
                                   innerProduct(g.expression(), lambda * v),
                               Kh2, INTEGRAL_BOUNDARY, In);

#endif

#elif defined(neumann1)
            // solve with Neumann BCs
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
#ifdef defined(classical)
            convdiff.addBilinear(+innerProduct((vel * n) * u, v), Kh2, INTEGRAL_BOUNDARY, In);
#endif

            convdiff.addLinear(+innerProduct(g_Neumann.expression(), v), Kh2, INTEGRAL_BOUNDARY, In);
#endif

#endif

            // Add RHS on bulk
            convdiff.addLinear(+innerProduct(f.expr(), v), Kh2, In);

            // Compute integrals
            Expression2 bhexp(b0h, 0, op_id);
            Expression2 ghexp(g, 0, op_id);
            Expression2 gdx(g, 0, op_dx);
            Expression2 gdy(g, 0, op_dy);
            // Expression2 normalx(n, 0, op_id);
            // Expression2 normaly(n, 1, op_id);

            intF = integral(Kh2, In, f, 0, qTime);
#ifdef neumann
            intG = integral(g_Neumann, In, interface, 0);
#elif defined(dirichlet)
            // double intu = integralSurf(b0h, In, qTime);
            // double intg = integralSurf(g, In, qTime);
            // double intDiff = lambda*(intu-intg);
            // double intGrad = integral<Mesh>((A2, grad(u)*n), b0h, In, qTime);
            double intu       = integral(bhexp, In, interface, 0);
            double intg       = integral(ghexp, In, interface, 0);
            double intVelb0h1 = 0.5 * integral(fabs(vel * n) * bhexp, In, interface, 0);
            double intVelb0h2 = 0.5 * integral((vel * n) * bhexp, In, interface, 0);
            double intVelg1   = 0.5 * integral(fabs(vel * n) * ghexp, In, interface, 0);
            double intVelg2   = 0.5 * integral((vel * n) * ghexp, In, interface, 0);
            double intVel     = intVelb0h1 - intVelb0h2 - (intVelg1 - intVelg2);
            double intGrad    = integral(A2 * gdx * n.x, In, interface, 0) + integral(A2 * gdy * n.y, In, interface, 0);
#endif

            // #ifdef conservative
            // #ifdef fullstab
            // #ifdef use_t
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_conservative_fullstab_h" +
            //                 std::to_string(h)
            //                 + "_"
            //                                                   + std::to_string(j + 1) + ".dat");
            // #else
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_conservative_fullstab_j" +
            //                 std::to_string(j + 1) + ".dat");
            // #endif
            // #elif defined(macro)
            // #ifdef use_t
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_conservative_macro_h" +
            //                 std::to_string(h) +
            //                 "_"
            //                                                   + std::to_string(j + 1) + ".dat");
            // #else
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_conservative_macro_j" + std::to_string(j
            //                 + 1)
            //                 + ".dat");
            // #endif
            // #endif
            // #elif defined(classical)
            // #ifdef fullstab
            // #ifdef use_t
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_classical_fullstab_h" +
            //                 std::to_string(h) +
            //                 "_"
            //                                                   + std::to_string(j + 1) + ".dat");
            // #else
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_classical_fullstab_j" + std::to_string(j
            //                 + 1)
            //                 + ".dat");
            // #endif
            // #elif defined(macro)
            // #ifdef use_t
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_,
            //                                path_output + "mat_classical_macro_h" + std::to_string(h) + "_" +
            //                                std::to_string(j + 1) + ".dat");
            // #else
            //             if (iter == GTime::total_number_iteration - 1)
            //                 matlab::Export(convdiff.mat_, path_output + "mat_classical_macro_j" + std::to_string(j +
            //                 1) +
            //                 ".dat");
            // #endif
            // #endif
            // #endif

            // Solve linear system
            convdiff.solve("umfpack");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Compute error in Reynold relation
            // Reynold error = [0.0199981, 0.00500015, 0.00125001]  // conservative
            // Reynold error = [0.0199981, 0.00500015, 0.00125001]  // classical
            {
                CutFEM<Mesh> reynold(qTime);
                reynold.initSpace(Wh, In);

                reynold.addBilinear(innerProduct(u, v), Kh2, (int)lastQuadTime, In);
                for (int i = 0; i < nbTime; ++i) {
                    reynold.addBilinear(innerProduct(dt(u), v) + innerProduct(u, dt(v)) +
                                            innerProduct((vel[i].exprList() * grad(u)), v) +
                                            innerProduct(u, (vel[i].exprList() * grad(v))),
                                        Kh2, In, i);
                }
                reynold.addLinear(innerProduct(b0h.expr(), v), Kh2, 0, In);
                int N = Wh.NbDoF();
                Rn lhs(2 * N);
                multiply(2 * N, 2 * N, reynold.mat_, data_u0, lhs);

                lhs -= reynold.rhs_;

                reynold_error.at(j) = lhs.linfty();
                std::cout << " e_r^n = " << reynold_error.at(j) << std::endl;
            }

            // Compute conservation error
            if (iterations == 1) {
                Fun_h funuh(Wh, data_u0);

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_u0(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Kh2, funsol, 0, 0);
                sol2 += data_u0(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Kh2, funsol, 0, lastQuadTime);

                if (iter == 0) {
                    q0_0 = q_0;
                    q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Kh2, b0h, 0, 0);
                }

                output_data << std::setprecision(10);
                output_data << current_time << "," << (q_1 - qp_1) << "," << intF << ","
#ifdef neumann
                            << intG << "," << ((q_1 - qp_1) - intF - intG) << std::endl;
#elif defined(dirichlet)
                            << intGrad << "," << intVel << "," << lambda * (intu - intg) << ","
                            << h * ((q_1 - qp_1) - intF - intGrad + intVel) + lambda * (intu - intg) << std::endl;
#endif
                qp_1 = q_1;
            }

            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));
            Fun_h funuh(Wh, sol);
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time, 0, 1);
            std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << std::endl;

            sol += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
            errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << std::endl;

            // R errBulk =  L2normCut(b0h, fun_uBulkD, GTime::current_time(), 0,
            // 1); std::cout << std::endl; std::cout << " L2 Error \t : \t" <<
            // errBulk << std::endl;

            // errors_bulk.at(j) = errBulk;
            errors_bulk.at(j) = errBulk;

            if ((iterations == 1)) {
                Paraview<Mesh> writerTh(Th, path_figures + "BackgroundMesh.vtk");
                Paraview<Mesh> writer(Kh2, path_figures + "Bulk" + std::to_string(iter + 1) + "DG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                // writer.add(ls[2], "levelSet2", 0, 1);
            }

            if (iterations > 1 && iter == total_number_iteration - 1)
                output_data << h << "," << dT << "," << errBulk << std::endl;

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

    std::cout << std::endl;
    std::cout << "Errors Bulk = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << errors_bulk.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << std::endl;
    std::cout << "Number of stabilized edges = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << numb_stabilized_edges.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << std::endl;
    std::cout << "Reynold error = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << reynold_error.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << std::endl;

    std::cout << "dT = [";
    for (int i = 0; i < iterations; i++) {
        std::cout << dts.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]"
              << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "[ms]" << std::endl;

    return 0;
}
