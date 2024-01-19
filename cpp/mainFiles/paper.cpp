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
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../num/matlab.hpp"
#include "../num/redirectOutput.hpp"
#include "paraview.hpp"
#include "../problem/AlgoimIntegration.hpp"

using namespace globalVariable; // to access some globally defined constants

// Numerical examples

// Numerical examples
namespace Example1 {
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

    // Automatic

    return D * 3.141592653589793 * cos(x * 3.141592653589793) * sin(y * 3.141592653589793) * 1.0 /
               sqrt(pow(y + cos(t * 3.141592653589793) * (7.0 / 2.5E+1) - 1.0 / 2.0, 2.0) +
                    pow(-x + sin(t * 3.141592653589793) * (7.0 / 2.5E+1) + 1.0 / 2.0, 2.0)) *
               (pow(cos(t * 3.141592653589793), 2.0) * 2.0 - 1.0) *
               (y * 2.0 + cos(t * 3.141592653589793) * (1.4E+1 / 2.5E+1) - 1.0) * (-1.0 / 5.0) +
           (D * 3.141592653589793 * cos(y * 3.141592653589793) * sin(x * 3.141592653589793) * 1.0 /
            sqrt(pow(y + cos(t * 3.141592653589793) * (7.0 / 2.5E+1) - 1.0 / 2.0, 2.0) +
                 pow(-x + sin(t * 3.141592653589793) * (7.0 / 2.5E+1) + 1.0 / 2.0, 2.0)) *
            (pow(cos(t * 3.141592653589793), 2.0) * 2.0 - 1.0) *
            (x * -2.0 + sin(t * 3.141592653589793) * (1.4E+1 / 2.5E+1) + 1.0)) /
               5.0;

    // Correct sign of normal vector at interface
    return -(pi * cos(2 * pi * t) * cos(pi * y) * sin(pi * x) * (x - 1. / 2 - 0.28 * sin(pi * t))) /
               (250 * sqrt((pow(x - 1. / 2 - 0.28 * sin(pi * t), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2)))) -
           (pi * cos(2 * pi * t) * cos(pi * x) * sin(pi * y) * (y - 0.5 + 0.28 * cos(pi * t))) /
               (250 * sqrt((pow((x - 1. / 2 - 0.28 * sin(pi * t)), 2) + pow(y - 0.5 + 0.28 * cos(pi * t), 2))));
}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
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

    // automatic
    return 3.141592653589793 * cos(x * 3.141592653589793) * cos(y * 3.141592653589793) *
               sin(t * 3.141592653589793 * 2.0) * (-4.0 / 5.0) -
           (3.141592653589793 * 3.141592653589793) * cos(t * 3.141592653589793 * 2.0) * cos(x * 3.141592653589793) *
               sin(y * 3.141592653589793) * (x - 1.0 / 2.0) * (2.0 / 5.0) +
           (3.141592653589793 * 3.141592653589793) * cos(t * 3.141592653589793 * 2.0) * cos(y * 3.141592653589793) *
               sin(x * 3.141592653589793) * (y - 1.0 / 2.0) * (2.0 / 5.0) +
           D * (3.141592653589793 * 3.141592653589793) * cos(t * 3.141592653589793 * 2.0) * cos(x * 3.141592653589793) *
               cos(y * 3.141592653589793) * (4.0 / 5.0);

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

namespace Example2 {
// Kite shape
const double R0       = .5;
const double C        = 5. / 3;
const double W        = 1. / 6;
const double D        = 0.01;
const double beta_max = W;

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

    return D * ((2 * pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t)) /
                    (R0 * R0) +
                (pi * pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
                 pow((2 * x - 2 * t * (-C * y * y + W)), 2)) /
                    (R0 * R0 * R0 * R0) +
                (pi * sin((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
                 (8 * C * C * t * t * y * y + 4 * C * t * (x - t * (-C * y * y + W)) + 2)) /
                    (R0 * R0) +
                (pi * pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * sin(pi * t) *
                 pow((2 * y + 4 * C * t * y * (x - t * (-C * y * y + W))), 2)) /
                    (R0 * R0 * R0 * R0)) +
           pi * cos((pi * (pow((x - t * (-C * y * y + W)), 2) + y * y)) / (R0 * R0)) * cos(pi * t) +
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
R fun_velocity(double *P, const int i, const double t) {
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
} // namespace Example2

namespace Example3 {
// N-sphere
const double R0       = .5;
const double D        = .01;
const double beta_max = 2.;

double fun_one(double *P, const int i) { return 1.; }

double rho(double *P, const double t) { return 1. / pi * sin(2 * pi * t); }

double fun_levelSet(double *P, const int i, const double t) {
    return (P[0] - rho(P, t)) * (P[0] - rho(P, t)) + P[1] * P[1] - R0 * R0;
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] - 1. / pi * sin(2 * pi * t)) * (P[0] - 1. / pi * sin(2 * pi * t)) + P[1] * P[1] - R0 * R0;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2. * (x(0) - 1. / pi * sin(2 * pi * t)), 2. * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(4. * (P[0] - 1. / pi * sin(2 * pi * t)) * (P[0] - 1. / pi * sin(2 * pi * t)) + 4. * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2. * (P[0] - 1. / pi * sin(2 * pi * t)) / norm, 2. * P[1] / norm);
    }
};

R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return D * ((4 * pi * sin(pi * t) *
                 sin((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) /
                     (R0 * R0))) /
                    (R0 * R0) +
                (4 * y * y * pi * pi * sin(pi * t) *
                 cos((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) /
                     (R0 * R0))) /
                    (R0 * R0 * R0 * R0) +
                (pi * pi * sin(pi * t) *
                 cos((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) /
                     (R0 * R0)) *
                 pow((2 * x - (5734161139222659 * sin(2 * pi * t)) / 9007199254740992), 2)) /
                    (R0 * R0 * R0 * R0)) +
           pi * cos(pi * t) *
               cos((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) /
                   (R0 * R0)) -
           (2 * pi * cos(2 * pi * t) * sin(pi * t) *
            sin((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) / (R0 * R0)) *
            (2 * x - (5734161139222659 * sin(2 * pi * t)) / 9007199254740992)) /
               (R0 * R0) +
           (5734161139222659 * pi * pi * cos(2 * pi * t) * sin(pi * t) *
            sin((pi * (pow((x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984), 2) + y * y)) / (R0 * R0)) *
            (x - (5734161139222659 * sin(2 * pi * t)) / 18014398509481984)) /
               (4503599627370496 * (R0 * R0));
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    if (i == 0)
        return 2 * cos(2 * pi * t);
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * sin(2 * pi * t);
    double r   = (x - rho) * (x - rho) + y * y;
    return cos(M_PI * (r) / (R0 * R0)) * sin(M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * sin(2 * pi * t);
    double r   = (x - rho) * (x - rho) + y * y;
    return cos(M_PI * (r) / (R0 * R0)) * sin(M_PI * t);
}
} // namespace Example3

namespace Preuss {
const double R0       = .5;
const double D        = 1.;
const double beta_max = 2.;

/* This works for running Test – i.e. a pure bulk problem on Omega_2. */

double fun_one(double *P, const int i) { return 1.; }

double fun_oneD(double *P, const int i, const int d, const R t) { return 1.; }

double rho(double *P, const double t) { return 1. / pi * std::sin(2 * pi * t); }

double fun_levelSet(double *P, const int i, const double t) {
    return P[0] * P[0] + (P[1] - rho(P, t)) * (P[1] - rho(P, t)) - R0 * R0;
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {

        return P[0] * P[0] + (P[1] - 1. / pi * sin(2 * pi * t)) * (P[1] - 1. / pi * sin(2 * pi * t)) - R0 * R0;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2. * x(0), 2. * (x(1) - 1. / pi * sin(2 * pi * t)));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        // R norm = sqrt(4. * (P[0] - 1. / pi * sin(2 * pi * t)) * (P[0] - 1. / pi * sin(2 * pi * t)) + 4. * P[1] *
        // P[1]);
        R norm = sqrt(4. * P[0] * P[0] + 4. * (P[1] - 1. / pi * sin(2 * pi * t)) * (P[1] - 1. / pi * sin(2 * pi * t)));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2. * P[0] / norm, 2. * (P[1] - 1. / pi * sin(2 * pi * t)) / norm);
    }
};

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) { return 0.; }

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    if (i == 0)
        return 0.;
    else
        return 2 * std::cos(2 * pi * t);
}

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

} // namespace Preuss

namespace Kite {

const double R0       = 1.;
const double C        = 1.;
const double W        = 1.;
const double D        = 1.;
const double beta_max = W;

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

    return D * ((3.141592653589793 *
                 sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * 1.0 / sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) /
                    R0 +
                (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) *
                 cos((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * pow(x * 2.0 - t * (W - C * (y * y)) * 2.0, 2.0)) /
                    (pow(x - t * (W - C * (y * y)), 2.0) * 4.0 + (y * y) * 4.0 + 4.0E-16) +
                (3.141592653589793 *
                 sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * 1.0 / sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16) *
                 ((C * C) * (t * t) * (y * y) * 8.0 + C * t * (x - t * (W - C * (y * y))) * 4.0 + 2.0)) /
                    (R0 * 2.0) -
                (3.141592653589793 *
                 sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * pow(y * 2.0 + C * t * y * (x - t * (W - C * (y * y))) * 4.0, 2.0) * 1.0 /
                 pow(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16, 3.0 / 2.0)) /
                    (R0 * 4.0) +
                (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) *
                 cos((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * pow(y * 2.0 + C * t * y * (x - t * (W - C * (y * y))) * 4.0, 2.0)) /
                    (pow(x - t * (W - C * (y * y)), 2.0) * 4.0 + (y * y) * 4.0 + 4.0E-16) -
                (3.141592653589793 *
                 sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * 3.141592653589793) * pow(x * 2.0 - t * (W - C * (y * y)) * 2.0, 2.0) * 1.0 /
                 pow(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16, 3.0 / 2.0)) /
                    (R0 * 4.0)) +
           3.141592653589793 *
               cos((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
               cos(t * 3.141592653589793) +
           (3.141592653589793 *
            sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
            sin(t * 3.141592653589793) * (W - C * (y * y)) * (x - t * (W - C * (y * y))) * 1.0 /
            sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) /
               R0 -
           (3.141592653589793 *
            sin((3.141592653589793 * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
            sin(t * 3.141592653589793) * (W - C * (y * y)) * (x * 2.0 - t * (W - C * (y * y)) * 2.0) * 1.0 /
            sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) /
               (R0 * 2.0);
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
R fun_velocity(double *P, const int i, const double t) {
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

    double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    return cos(M_PI * r / R0) * sin(M_PI * t);

    // double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    // return cos(M_PI * (r + Epsilon) / (R0));
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    return cos(M_PI * r / R0) * sin(M_PI * t);

    // double r = std::sqrt((x - (W - C * y * y) * t) * (x - (W - C * y * y) * t) + y * y);
    // return cos(M_PI * (r + Epsilon) / (R0));
}
} // namespace Kite

namespace N_sphere {
// N-sphere
const double R0       = .5;
const double D        = 1.;
const double beta_max = 2.;

double fun_one(double *P, const int i) { return 1.; }

double rho(double *P, const double t) { return 1. / pi * sin(2 * pi * t); }

double fun_levelSet(double *P, const int i, const double t) {
    return (P[0] - rho(P, t)) * (P[0] - rho(P, t)) + P[1] * P[1] - R0 * R0;
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] - 1. / pi * sin(2 * pi * t)) * (P[0] - 1. / pi * sin(2 * pi * t)) + P[1] * P[1] - R0 * R0;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2. * (x(0) - 1. / pi * sin(2 * pi * t)), 2. * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(4. * (P[0] - 1. / pi * sin(2 * pi * t)) * (P[0] - 1. / pi * sin(2 * pi * t)) + 4. * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2. * (P[0] - 1. / pi * sin(2 * pi * t)) / norm, 2. * P[1] / norm);
    }
};

R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return D * ((3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                     R0) *
                 sin(t * 3.141592653589793) * 1.0 /
                 sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16) * 2.0) /
                    R0 -
                (3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                     R0) *
                 sin(t * 3.141592653589793) *
                 pow(x * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1, 2.0) * 1.0 /
                 pow(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16,
                     3.0 / 2.0)) /
                    (R0 * 4.0) +
                (1.0 / (R0 * R0) * (y * y) * (3.141592653589793 * 3.141592653589793) *
                 cos((3.141592653589793 *
                      sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                     R0) *
                 sin(t * 3.141592653589793)) /
                    (pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16) +
                (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) *
                 cos((3.141592653589793 *
                      sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                     R0) *
                 sin(t * 3.141592653589793) *
                 pow(x * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1, 2.0)) /
                    (pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) * 4.0 + (y * y) * 4.0 +
                     4.0E-16) -
                ((y * y) * 3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                     R0) *
                 sin(t * 3.141592653589793) * 1.0 /
                 pow(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16,
                     3.0 / 2.0)) /
                    R0) +
           3.141592653589793 *
               cos((3.141592653589793 *
                    sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                   R0) *
               cos(t * 3.141592653589793) -
           (3.141592653589793 *
            sin((3.141592653589793 *
                 sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                R0) *
            cos(t * 3.141592653589793 * 2.0) * sin(t * 3.141592653589793) *
            (x * 2.0 - sin(t * 3.141592653589793 * 2.0) * 6.366197723675814E-1) * 1.0 /
            sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
               R0 +
           ((3.141592653589793 * 3.141592653589793) *
            sin((3.141592653589793 *
                 sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16)) /
                R0) *
            cos(t * 3.141592653589793 * 2.0) * sin(t * 3.141592653589793) *
            (x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1) * 1.0 /
            sqrt(pow(x - sin(t * 3.141592653589793 * 2.0) * 3.183098861837907E-1, 2.0) + y * y + 1.0E-16) *
            6.366197723675814E-1) /
               R0;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    if (i == 0)
        return 2 * cos(2 * pi * t);
    else
        return 0.;
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i) { return 0.; }

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * sin(2 * pi * t);
    double r   = std::sqrt((x - rho) * (x - rho) + y * y);
    return cos(M_PI * r / R0) * sin(M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double rho = 1. / pi * sin(2 * pi * t);
    double r   = std::sqrt((x - rho) * (x - rho) + y * y);
    return cos(M_PI * r / R0) * sin(M_PI * t);
}
} // namespace N_sphere

namespace Example4 {
// A circle moving in a circular trajectory
const double D        = 1.;
const double R0       = 0.17;
const double beta_max = M_PI * 0.5;

// double fun_tau(double *P, const int i, const R t) {
//     const double beta_2 = std::sqrt(std::pow(M_PI * (0.5 - P[1]), 2) + std::pow(M_PI * (P[0] - 0.5), 2));
//     const double Pe       = beta_2 * h / (2 * D);
//     const double xi       = 1. / std::tanh(Pe) - 1. / Pe;
//     const double tau_supg = h / (2 * beta_2) * xi;
// }

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
R fun_velocity(double *P, const int i, const double t) {
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

    double r = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
    return cos(M_PI * r / R0) * sin(M_PI * t);
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    double x = P[0], y = P[1];

    double xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);

    double r = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
    return cos(M_PI * r / R0) * sin(M_PI * t);
}

// RHS fB bulk
R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    // automatic
    return D * ((3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * 3.141592653589793) * 1.0 /
                 sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                      pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0) *
                 2.0E+8) /
                    R0 +
                (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) * sin(t * 3.141592653589793) *
                 cos((3.141592653589793 *
                      sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12) /
                    (pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                     pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0) -
                (3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * 3.141592653589793) * pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) *
                 1.0 /
                 pow(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0,
                     3.0 / 2.0) *
                 4.0E+20) /
                    R0 -
                (3.141592653589793 *
                 sin((3.141592653589793 *
                      sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * 3.141592653589793) * pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) *
                 1.0 /
                 pow(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0,
                     3.0 / 2.0) *
                 4.0E+20) /
                    R0 +
                (1.0 / (R0 * R0) * (3.141592653589793 * 3.141592653589793) * sin(t * 3.141592653589793) *
                 cos((3.141592653589793 *
                      sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12) /
                    (pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                     pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) +
           3.141592653589793 * cos(t * 3.141592653589793) *
               cos((3.141592653589793 *
                    sqrt(pow(y * 5.0E+1 + cos(t * 3.141592653589793) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * 3.141592653589793) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                   (R0 * 1.0E+8));
}

} // namespace Example4

//* Set numerical example (options: "ex1", "ex2", or "ex3")
#define ex1
//* Set scheme for the method (options: "classical", "conservative")
#define classical
//* Set stabilization method (options: "fullstab", "macro")
#define fullstab

#define use_h    // to set mesh size using the h parameter. Write use_n to decide
                 // using nx, ny.
#define use_tnot // write use_t to control dT manually. Otherwise it is set
                 // proportional to h.
#define conservationnot

#define reynoldnot

#define local   // compute and print local errors (meaning errors in all time instances during [0,T])

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
//using namespace Example2;
using namespace Kite;
#elif defined(ex3)
// using namespace Example3;
using namespace N_sphere;
#elif defined(preuss)
using namespace Preuss;
#elif defined(ex4)
using namespace Example4;
#endif

int main(int argc, char **argv) {
    // MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 1;         // number of mesh refinements   (set to 1 to run
                                         // only once and plot to paraview)
    int nx = 15, ny = 15;                // starting mesh size
    double h  = 0.1/2; // starting mesh size
    double dT = 0.25;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#ifdef ex1
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/example1/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/example1/paraview/";
#elif defined(ex2)
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/example2/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/example2/paraview/";
#elif defined(ex3)
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/example3/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/example3/paraview/";
#elif defined(ex4)
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/example4/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/example4/paraview/";
#elif defined(preuss)
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/preuss/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/preuss/paraview/";

#endif

    // Create directory if not already existent
    // if (MPIcf::IamMaster()) {
    // std::filesystem::create_directories(path_output_data);
    // std::filesystem::create_directories(path_output_figures);
    // }

    // Data file to hold problem data
    std::ofstream output_data(path_output_data + "data_temp.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors, errors_T, hs, nxs, ns_active, dts, omega, gamma, global_conservation_errors,
        reynold_error;
    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {       

        // Define background mesh
#if defined(ex1)
        const double x0 = 0.-Epsilon, y0 = 0. - Epsilon, lx = 1., ly = 1.;
        //const double x0 = 0. - 2.01 * Example1::R0, y0 = -2.01 * Example1::R0, lx = 4. * Example1::R0, ly = 4. * Example1::R0;
        // const double lx = 2., ly = 1.1;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        // Mesh Th(nx, ny, 0. - Epsilon, 0. - Epsilon, lx, ly);
        //Mesh Th(nx, ny, x0, y0, lx, ly);

#elif defined(ex2)
        //const double lx = 6. * Example2::R0, ly = 6. * Example2::R0, x0 = -3.01 * Example2::R0, y0 = -3.01 * Example2::R0;
        const double lx = 7., ly = 3., x0 = -3.5-Epsilon, y0 = -1.5 - Epsilon; // like paper
        //const double lx = 6., ly = 4., x0 = -2.-Epsilon, y0 = -2. - Epsilon;     // works like the above
        //const double lx = 6., ly = 6., x0 = -3.-Epsilon, y0 = -3. - Epsilon;    // works better? for some reason
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        // Mesh Th(nx, ny, -3.01 * Example2::R0, -3.01 * Example2::R0, lx, ly);
        // Mesh Th(nx, ny, -3.500123, -1.5, 7., 3.);
        //Mesh Th(nx, ny, -2.500123, -1.5, 7., 3.);
        // Mesh Th(12, 11, -.85, -0.8, 2.05, 1.6);
#elif defined(ex3)
        // const double lx = 4. * Example3::R0, ly = 3. * Example3::R0;
        //  const double lx = 2., ly = 1.1;
        const double lx = 2., ly = 1.2;
        double x0 = -1. - Epsilon, y0 = -0.6 - Epsilon;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        // Mesh Th(nx, ny, -2.01 * Example3::R0, -1.5 * Example3::R0, lx, ly);
        // Mesh Th(nx, ny, -1.-Epsilon, -0.6-Epsilon, lx, ly);

#elif defined(ex4)
        const double lx = 1., ly = 1.;
        const double x0 = 0. - Epsilon, y0 = 0. - Epsilon;
        // const double lx = 4. * Example4::R0, ly = 4. * Example4::R0;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        // Mesh Th(nx, ny, 0. - Epsilon, 0. - Epsilon, lx, ly);
        // Mesh Th(nx, ny, 0. - 2.01 * Example4::R0, -2.01 * Example4::R0, lx, ly);

#elif defined(preuss)
        const double lx = 3. * Preuss::R0, ly = 4. * Preuss::R0;
        const double x0 = -1.51 * Preuss::R0, y0 = -1.5 * Preuss::R0;
        // const double lx = 2., ly = 1.1;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        // Mesh Th(nx, ny, -1.51 * Preuss::R0, -1.5 * Preuss::R0, lx, ly);
#endif

        Mesh Th(nx, ny, x0, y0, lx, ly);

        // Parameters
        const double tfinal = .5; // Final time
                                  // const double tfinal = .5; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 3;

        double dT = h / divisionMeshSize;
        //dT = 0.0078125 / std::sqrt(2);

        //dT = 0.5 * std::pow(2, -j-1 - 1);
        //  double dT = tfinal / 32;

        // double dT = tfinal / (j + 1);

        total_number_iteration = int(tfinal / dT);
#endif
        dT        = tfinal / total_number_iteration;
        time_step = dT;

        hs.at(j)  = h;
        nxs.at(j) = nx;
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

        // const double D = N_sphere::D;

        // CG stabilization parameter
        // const double tau1 = 0.1 * (D + beta_max), tau2 = 0.1 * (D + beta_max);
        const double tau1 = .1, tau2 = 0.1;

        FESpace Vh(Th, DataFE<Mesh>::P3); // Background FE Space
        // FESpace Vh_interpolation(Th, DataFE<Mesh>::P3); // for interpolating data

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P3Poly); // FE Space in time
        // FESpace1 Ih_interpolation(Qh, DataFE<Mesh1>::P3Poly); // for interpolating data

        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(7)); // specify order of quadrature in time
        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

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

        Levelset<2> phi;
        ProblemOption option;
        const int quadrature_order_space       = 9;
        option.order_space_element_quadrature_ = quadrature_order_space;
        AlgoimCutFEM<Mesh, Levelset<2>> convdiff(qTime, phi, option);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous, mass_initial;
        double intF = 0, int_outflow = 0, intG = 0, intF_total = 0,
               intG_total                = 0; // hold integrals of rhs and Neumann bcs
        double global_conservation_error = 0, local_conservation_error = 0, errBulk = 0., error_I = 0.;
        std::vector<double> local_errors, local_conservation_errors, global_conservation_errors_t;

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

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            ns_active.at(j) = Thi.get_nb_element();

            // Right hand side functions
            // Fun_h f(Vh_interpolation, In_interpolation, fun_rhsBulk);
            // Fun_h g_Neumann(Vh_interpolation, In_interpolation, fun_neumann_Gamma); // computer Neumann BC

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

            // Variational formulation
#if defined(classical)
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);
            // Impose initial condition
            if (iter == 0) {
                convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
            } else {
                convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            }
            convdiff.addBilinear(+innerProduct(dt(u), v) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(+innerProduct((vel[i].exprList() * grad(u)), v), Thi, In, i);
            }

#elif defined(conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);

            // Impose initial condition
            if (iter == 0) {
                convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
            } else {
                convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            }

            convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Thi, In, i);
            }

#endif

            // Source function
            // convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);
            convdiff.addLinearExact(fun_rhsBulk, +innerProduct(1, v), Thi, In);

            // Neumann boundary condition
            // convdiff.addLinear(+innerProduct(g_Neumann.expr(), v), interface, In);
            convdiff.addLinearExact(fun_neumann_Gamma, +innerProduct(1, v), interface, In);

            // Stabilization
            double stab_bulk_faces = 0.;
            double stab_mass       = 0.; // tau1 * h;
            double stab_dt         = 0.; // tau1 * h;

#if defined(fullstab)

            // convdiff.addFaceStabilization(
            //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
            //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In);

            convdiff.addPatchStabilization(+innerProduct(tau1 / h / h * jump(u), jump(v)), Thi, In);
            // convdiff.addPatchStabilization(+innerProduct(tau1 / 10000 / h * jump(u), jump(v)), Thi, In, 0);

#elif defined(macro)
            // MacroElementPartition<Mesh> TimeMacro(Thi, 0.3);
            // std::cout << TimeMacro.number_of_stabilized_edges << "\n";
            // std::cout << "number of stabilized edges: " << convdiff.get_number_of_stabilized_edges() << "\n";

            AlgoimMacro<Mesh, Levelset<2>> TimeMacro(Thi, 0.5, phi, In, qTime);
            // convdiff.addFaceStabilization(
            //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
            //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In, TimeMacro);

            convdiff.addPatchStabilization(+innerProduct(tau1 / h / h * jump(u), jump(v)), Thi, In, TimeMacro);

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

            if (iter == total_number_iteration - 1) {
                matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
            }

            // Solve linear system
            convdiff.solve("umfpack");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            Fun_h uh_t(Wh, In, data_u0); // FEM function in Pk(In) x Lagrange_m(Omega)
            // Rn zrs(convdiff.get_nb_dof(), 0.); // initial data total
            // Fun_h zrs_fun(Wh, In, zrs);
            // Rn halves(convdiff.get_nb_dof(), 0.5); // initial data total
            // Fun_h halves_fun(Wh, In, halves);
            // error_I = std::pow(L2_norm_T(zrs_fun, fun_oneD, Thi, In, qTime, phi), 2) / dT; // int_In ||u(t) -
            // error_I = std::pow(L2_norm_T(halves_fun, fun_oneD, Thi, In, qTime, phi), 2) / dT; // int_In ||u(t) -
            // u_h(t)||_{Omega(t)} dt
            error_I += L2_norm_T(uh_t, fun_uBulkD, Thi, In, qTime, phi,
                                 quadrature_order_space); // int_In ||u(t) - u_h(t)||_{Omega(t)} dt

            std::cout << " t_n -> || u-uex||_(In x Omega)^2 = " << error_I << '\n';

// Compute error in Reynold relation
#if defined(reynold)
            AlgoimCutFEM<Mesh, Levelset<2>> reynold(qTime, phi);

            reynold.initSpace(Wh, In);

            reynold.addBilinear(innerProduct(u, v), Thi, (int)lastQuadTime, In);
            if (iter == 0)
                reynold.addLinearExact(fun_uBulk, innerProduct(1, v), Thi, 0, In);
            else
                reynold.addLinear(innerProduct(b0h.expr(), v), Thi, 0, In);

            reynold.addBilinear(-innerProduct(dt(u), v) - innerProduct(u, dt(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                reynold.addBilinear(-innerProduct((vel[i].exprList() * grad(u)), v) -
                                        innerProduct(u, (vel[i].exprList() * grad(v))),
                                    Thi, In, i);
            }

            int N = Wh.NbDoF();
            Rn lhs(ndfTime * N);
            multiply(ndfTime * N, ndfTime * N, reynold.mat_, data_u0, lhs);

            lhs -= reynold.rhs_;

            reynold_error.at(j) = lhs.linfty();

            std::cout << " e_r^n = " << reynold_error.at(j) << '\n';
#endif

            // Compute area of domain in time quadrature point 0
            Fun_h funone(Wh, fun_one);

            // double intGamma = integral_algoim(funone, In, interface, phi, 0, quadrature_order_space) / dT;
            // double intOmega = integral_algoim(funone, Thi, phi, In, qTime, lastQuadTime, quadrature_order_space);

            // gamma[j] = intGamma;
            // omega[j] = intOmega;

            // Compute error of numerical solution
            Rn sol(Wh.get_nb_dof(), 0.);
            Fun_h funuh(Wh, sol);

            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));

            // errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, 0, phi, 0, 1, quadrature_order_space);
            // std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << '\n';

            for (int n = 1; n < ndfTime; n++) {
                sol += data_u0(SubArray(Wh.get_nb_dof(), n * Wh.get_nb_dof()));
            }

            errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1, quadrature_order_space);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            errors[j] = errBulk;
            local_errors.push_back(errBulk);

// Compute conservation error
// if (iterations == 1 && h > 0.01) {
#if defined(conservation)
            intF = integral_algoim(fun_rhsBulk, 0, Thi, phi, In, qTime,
                                   quadrature_order_space); // integrate source over In
            intG = integral_algoim(fun_neumann_Gamma, In, interface, phi, 0,
                                   quadrature_order_space); // integrate flux boundary over In

            intF_total += intF;
            intG_total += intG;

            // intF = integral_algoim(f, 0, Thi, phi, In, qTime);        // integrate source over In
            // intG = integral_algoim(g_Neumann, In, interface, phi, 0); // integrate flux across boundary over In

            double mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime,
                                               quadrature_order_space); // mass in last quad point

            if (iter == 0) {
                mass_initial       = integral_algoim(fun_uBulk, Thi, phi, In, qTime, 0, quadrature_order_space);
                mass_last_previous = mass_initial;
                // mass_last_previous = integral_algoim(b0h, Thi, phi, In, qTime, 0);
            }

            local_conservation_error  = ((mass_last - mass_last_previous) - intF - intG);
            // global_conservation_error += local_conservation_error;
            global_conservation_error = (mass_last - mass_initial - intF_total - intG_total);

            std::cout << "global_conservation_error: " << global_conservation_error << "\n";
            std::cout << "local_conservation_error: " << local_conservation_error << "\n";

            // outputData << std::setprecision(10);
            // outputData << current_time << "," << (mass_last - mass_last_previous) << "," << intF << "," << intG <<
            // ","
            //            << local_conservation_error << '\n';

            mass_last_previous = mass_last; // set current last to previous last for next time slab
                                            //}

            global_conservation_errors[j] = std::fabs(global_conservation_error);
            local_conservation_errors.push_back(std::fabs(local_conservation_error));
            global_conservation_errors_t.push_back(std::fabs(global_conservation_error));

#endif

            if ((iterations == 1) && (h > 0.01)) {
                Fun_h sol_h(Wh, sol);
                Paraview<Mesh> writerTh(Th, path_output_figures + "Th.vtk");
                Paraview<Mesh> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) + ".vtk");
                writer.add(b0h, "bulk", 0, 1);
                writer.add(sol_h, "bulk_end", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);

                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                // writer.add(g_Neumann, "neumann", 0, 1);
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

            // if (iterations > 1 && iter == total_number_iteration - 1)
            //     outputData << h << "," << dT << "," << errBulk << '\n';

            iter++;
        }
        errors_T[j] = std::sqrt(error_I);

        if (iterations > 1) {
            output_data << h << "," << dT << "," << errBulk << "," << errors_T[j] << "," << ns_active.at(j) << '\n';
            output_data.flush();
        }

        std::cout << "error_T = " << errors_T[j] << "\n";

    #ifdef local
        std::cout << "\n";
        std::cout << "Local errors = [";
        for (auto &err : local_errors) {

            std::cout << err;

            std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "\n";
        std::cout << "Local conservation error = [";
        for (auto &err : local_conservation_errors) {

            std::cout << err;

            std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "\n";
        std::cout << "Global conservation errors = [";
        for (auto &err : global_conservation_errors_t) {

            std::cout << err;

            std::cout << ", ";
        }
        std::cout << "]\n";
    #endif

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

    std::cout << "ns_active_mesh = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << ns_active.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}
