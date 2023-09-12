/*
  Solve the time-dependent convection-diffusion equation along a 2D surface using a space-
  time CutFEM with algoim for quadrature.
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
#include <vector>
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
#include "../problem/AlgoimIntegration.hpp"
#include <ranges>
#include <numeric>

using namespace globalVariable;

namespace Example1 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return (P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - 0.17 * 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return (P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22) - 0.17 * 0.17;
}

template <int N> struct Levelset {

    double t;

    // level set function
    // template <typename T> T operator()(const algoim::uvector<T, N> &x) const { return (x(0) - 0.5 - 0.28 * sin(pi *
    // t)) * (x(0) - 0.5 - 0.28 * sin(pi * t)) + (x(1) - 0.5 + 0.28 * cos(pi * t)) * (x(1) - 0.5 + 0.28 * cos(pi * t)) -
    // 0.17*0.17;}
    template <typename V> typename V::value_type operator()(const V &x) const {
        return (x[0] - 0.5 - 0.28 * sin(pi * t)) * (x[0] - 0.5 - 0.28 * sin(pi * t)) +
               (x[1] - 0.5 + 0.28 * cos(pi * t)) * (x[1] - 0.5 + 0.28 * cos(pi * t)) - 0.17 * 0.17;
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
R fun_levelSet(double *P, const int i, const R t) { return (P[0] - t) * (P[0] - t) + P[1] * P[1] - 2; }
R fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - 2; }

// R fun_levelSet(double *P, const int i, const R t) { return (P[0] - t) * (P[0] - t) + P[1] * P[1] - 2 - Epsilon; }
// R fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - 2 - Epsilon; }

R fun_rhs(double *P, const int cc, const R t) {
    R x = P[0], y = P[1];

    return -1. / 4 * exp(-t / 4) *
           (y * ((x - t) * (x - t) + y * y - 4.) / (pow(((x - t) * (x - t) + y * y), 3. / 2))); // u = exp(-t/4)...

    // return -(y * exp(-4 * t) * (4 * t * t - 8 * t * x + 4 * x * x + 4 * y * y - 1)) /
    //        pow((t * t - 2 * t * x + x * x + y * y), 3. / 2);        // u = exp(-4t)...
}

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] - t) * (P[0] - t) + P[1] * P[1] - 2;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * (x(0) - t), 2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(4. * (P[0] - t) * (P[0] - t) + 4. * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2.0 * (P[0] - t) / norm, 2.0 * P[1] / norm);
    }
};

} // namespace Shi1

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
// R fun_levelSet(double *P, const int i, const R t) {
//     return (P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1] - 1 - Epsilon;
// }
// R fun_levelSet(double *P, const int i) { return sqrt((P[0] + 0.5) * (P[0] + 0.5) + P[1] * P[1]) - 1 - Epsilon; }

R fun_levelSet(double *P, const int i, R t) { return (P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1] - 1; }
R fun_levelSet(double *P, const int i) { return (P[0] + 0.5) * (P[0] + 0.5) + P[1] * P[1] - 1; }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename V> typename V::value_type operator()(const V &P) const {
        return (P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + P[1] * P[1] - 1;
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * (x(0) + 0.5 - 2 * t), 2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(std::span<double> P) const {
        R norm = sqrt(4. * (P[0] + 0.5 - 2 * t) * (P[0] + 0.5 - 2 * t) + 4. * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2.0 * (P[0] + 0.5 - 2 * t) / norm, 2.0 * P[1] / norm);
    }
};

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

// Set numerical example (options: "example1", "shi1", "shi2", "deckelnick", "deckelnick2")
#define example1
// Set scheme for the dg method (options: "conservative", "classical" see
// thesis. Irrelevant if "cg" is defined instead of "dg")
#define conservative

#define levelsetexact

#define use_h

// Setup two-dimensional class types
const int d = 2;

typedef MeshQuad2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

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

    ProblemOption option;
    option.order_space_element_quadrature_ = 7;

    // Mesh settings and data objects
    const size_t iterations = 3; // number of mesh refinements   (set to 1 to
                                 // run only once and plot to paraview)
    int nx = 20, ny = 15;        // starting mesh size (only apply if use_n is defined)
    // double h  = 0.1 * pow(0.5, 5) * sqrt(0.5); // starting mesh size
    double h  = 0.1; // starting mesh size
    double dT = 0.0625;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#if defined(frachon1)
    // Paths to store data
    const std::string path_output_data = "../output_files/surface/algoim/frachon1/data/";
    const std::string path_figures     = "../output_files/surface/algoim/frachon1/paraview/";
#elif defined(zahedi1)
    const std::string path_output_data = "../output_files/surface/algoim/zahedi1/data/";
    const std::string path_figures     = "../output_files/surface/algoim/zahedi1/paraview/";
#elif defined(zahedi2)
    const std::string path_output_data = "../output_files/surface/algoim/zahedi2/data/";
    const std::string path_figures     = "../output_files/surface/algoim/zahedi2/paraview/";
#elif defined(shi1)
    const std::string path_output_data = "../output_files/surface/algoim/shi1/data/";
    const std::string path_figures     = "../output_files/surface/algoim/shi1/paraview/";
#elif defined(shi2)
    const std::string path_output_data = "../output_files/surface/algoim/shi2/data/";
    const std::string path_figures     = "../output_files/surface/algoim/shi2/paraview/";
#elif defined(shi3)
    const std::string path_output_data = "../output_files/surface/algoim/shi3/data/";
    const std::string path_figures     = "../output_files/surface/algoim/shi3/paraview/";
#elif defined(example1)
    const std::string path_output_data = "../output_files/surface/algoim/example1/data/";
    const std::string path_figures     = "../output_files/surface/algoim/example1/paraview/";
#elif defined(deckelnick)
    const std::string path_output_data = "../output_files/surface/algoim/deckelnick/data/";
    const std::string path_figures     = "../output_files/surface/algoim/deckelnick/paraview/";
#elif defined(deckelnick2)
    const std::string path_output_data = "../output_files/surface/algoim/deckelnick2/data/";
    const std::string path_figures     = "../output_files/surface/algoim/deckelnick2/paraview/";
#elif defined(deckelnick2toshi1)
    const std::string path_output_data = "../output_files/surface/algoim/deckelnick2_to_shi1/data/";
    const std::string path_figures     = "../output_files/surface/algoim/deckelnick2_to_shi1/paraview/";
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
    std::array<double, iterations> average_errors;     // array to hold bulk errors
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
        Mesh Th(nx, ny, -3+Epsilon, -3+Epsilon, lx, ly);

        // std::string f = "../mesh/square_seb_"+std::to_string(j+1)+".msh";
        // Mesh Th(f.c_str());
#elif defined(deckelnick) || defined(deckelnick2)
        // const double lx = 4.8123, ly = 4.8353;
        //const double lx = 8.8123, ly = 4.8353;
        //const double lx = 10.8123, ly = 4.8353;
        //const double lx = 8.01, ly = 6.01; 
        //const double lx = 8., ly = 6.; // as in paper
        const double lx = 4.8, ly = 4.8;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif

        //Mesh Th(nx, ny, -2.435, -2.597, lx, ly);
        //Mesh Th(nx, ny, -3-Epsilon, -3-Epsilon, lx, ly);    // as in Shi paper
        Mesh Th(nx, ny, -2.4+Epsilon, -2.4+Epsilon, lx, ly);    // as in Deckelnick paper
#endif

        // Parameters
        double tfinal = .1; // Final time

        int divisionMeshSize = 3;

        double dT = h / divisionMeshSize;
        //double dT = h * h;

        total_number_iteration = int(tfinal / dT);
        
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
        double tau0 = 0, tau1 = 0.1, tau2 = 0.1;

        // Background FE Space, Time FE Space & Space-Time Space
        FESpace Vh(Th, DataFE<Mesh>::P3);       // Background FE space
        FESpace Vh_interpolation(Th, DataFE<Mesh>::P2);      // For interpolating data

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly);     // Time FE space
        FESpace1 Ih_interpolation(Qh, DataFE<Mesh1>::P2Poly);    // for interpolating data
        

        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(5));

        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Velocity field
        LagrangeQuad2 FEvelocity(1);

        FESpace VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Convection-Diffusion Problem Object
        Levelset<2> phi;    
        AlgoimCutFEM<Mesh, Levelset<2>> surfactant(qTime, phi, option);
        //AlgoimCutFEM<Mesh, Levelset<2>> surfactant(qTime, phi);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << "\n";

        int iter = 0;
        // double q0_0, q0_1, qp_0, qp_1;
        double q_init0, q_init1, qp1;
        double intF = 0, intG = 0, intGamma = 0; // hold integrals of rhs and Neumann bcs
        double errL2 = 0.;

        std::vector<double> gamma_length, error_t;
        std::vector<int> numb_stabilized_edges, numb_of_active_elements;
        
        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double tid            = iter * time_step;

            const TimeSlab &In(Ih[iter]);
            const TimeSlab &In_interpolation(Ih_interpolation[iter]);
            

            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << "\n";
            std::cout << " Time      \t : \t" << current_iteration * time_step << "\n";

            // computation of the interface in the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

                R tt = In.Pt(R1(qTime(i).x));
                phi.t = tt;
                interface.init<Levelset<2>>(i, Th, phi);
            }

            // Create active meshes
            ActiveMesh<Mesh> ThGamma(Th);
            
            ThGamma.createSurfaceMesh(interface);

            CutSpace Wh(ThGamma, Vh);

            // Data for initial solution
            surfactant.initSpace(Wh, In);

            // Interpolate data
            //Fun_h funrhs(Vh2, In2, fun_rhs);
            Fun_h funrhs(Vh_interpolation, In_interpolation, fun_rhs);
            
            Rn datau0(surfactant.get_nb_dof(), 0.);
            surfactant.initialSolution(datau0); 
            KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(), 0)));
            if (iter == 0)
                interpolate(Wh, datas0, fun_init_surfactant);
            
            Fun_h u0(Wh, datau0);


            if (iter == 0) {
                Paraview<Mesh> writer(ThGamma, path_figures + "surfactant_initial" + ".vtk");
                Fun_h uS_ex(Wh, fun_sol_surfactant, 0);
                writer.add(u0, "surfactant", 0, 1);
                writer.add(uS_ex, "surfactant_exact", 0, 1);
            }

            Normal n;
            FunTest u(Wh, 1), v(Wh, 1);

            // phi.t = tid;
            //  gnuplot::save(Th);
            //  gnuplot::save<Mesh, Levelset<2>>(ThGamma, *interface(0), phi, "interface.dat",tid);
            //  gnuplot::save<Mesh>(*interface(0));
            //  getchar();

#ifdef conservative
            
            surfactant.addBilinear(+innerProduct(u, v), *interface(lastQuadTime), In, lastQuadTime);
            surfactant.addLinear(+innerProduct(u0.expr(), v), *interface(0), In, 0);
 

            // surfactant.addBilinear(+innerProduct(0.1*dT*u, v), *interface(0), In, 0);
            // surfactant.addLinear(+innerProduct(0.1*dT*u0.expr(), v), *interface(0), In, 0);

            surfactant.addBilinear(
                - innerProduct(u, dt(v)) 
                + innerProduct(D * gradS(u), gradS(v)) 
                - innerProduct(u, (vel.exprList() * grad(v))),
                interface, In);

// classical scheme
#elif defined(classical)

            surfactant.addBilinear(innerProduct(u, v), *interface(0), In, 0);
            surfactant.addLinear(innerProduct(u0.expr(), v), *interface(0), In, 0);

            surfactant.addBilinear(
                + innerProduct(dt(u), v)
                + innerProduct(D * gradS(u), gradS(v)) 
                + innerProduct((vel.exprList() * grad(u)), v)
                , interface, In);

#endif

            // Stabilization
            double stab_surf_face      = h * tau1;
            double stab_surf_interface = h * h * tau2;
            double stab_mass           = 0.;//tau1 * h;
            double stab_dt             = 0.;//tau1 * dT;

            surfactant.addFaceStabilization(+innerProduct(tau1 * jump(grad(u) * n), jump(grad(v) * n))+
                    innerProduct(h * h * tau1 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
                                            ThGamma, In);
            
            // stabilize in last quadrature point
            double ccend = 1. / In.T.measure() * 1. / qTime[lastQuadTime].a;
            double ccmid = 1. / In.T.measure() * 1. / qTime[lastQuadTime - 1].a;
            
            surfactant.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), ccend * jump(grad(v) * n)),
                                            ThGamma, In, lastQuadTime);
            
            // surfactant.addFaceStabilizationSpecial(+innerProduct(stab_mass * jump(grad(u) * n), jump(grad(v) * n)),
            //                                 ThGamma, lastQuadTime, In);
            
            // surfactant.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), jump(grad(v) * n)),
            //                                 ThGamma, lastQuadTime, In);

            // surfactant.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), jump(grad(v) * n)),
            //                                 ThGamma, 0, In);
            // surfactant.addFaceStabilization(+innerProduct(stab_mass * jump(grad(u) * n), jump(grad(v) * n)),
            //                                 ThGamma, 1, In);

            

            surfactant.addFaceStabilization(-innerProduct(stab_dt * jump(grad(u) * n), jump(grad(dt(v)) * n)),
            ThGamma, In);

            // surfactant.addFaceStabilization(-innerProduct(stab_mass * jump(grad(u) * n), jump((vel.exprList() *
            // grad(v)))), ThGamma,
            //                                 In);

            surfactant.addBilinear(+ innerProduct(tau2 * h * h * grad(u) * n, grad(v) * n)
                                   + innerProduct(tau2 * h * h * h * h * grad(grad(u) * n) * n, grad(grad(v) * n) * n), interface, In);

            // Add RHS on surface
            surfactant.addLinear(+innerProduct(funrhs.expr(), v), interface, In);

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

                //KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
                //datau0 = dw;
                datau0 = surfactant.rhs_;

                surfactant.saveSolution(datau0);

                // Compute error
                Rn sol(Wh.get_nb_dof(), 0.);
                sol += datau0(SubArray(Wh.get_nb_dof(), 0));
                // sol += datau0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
                for (int n = 1; n < ndfTime; n++) {
                    sol += datau0(SubArray(Wh.get_nb_dof(), n * Wh.get_nb_dof()));
                }

                Fun_h funuh_0(Wh, datau0);
                Fun_h funuh(Wh, sol);

                Fun_h funone(Wh, fun_one, 0.);

                intF     = integral_algoim(funrhs, In, interface, phi, 0);
                intGamma = integral_algoim(funone, *interface(0), 0, phi, tid);
                // std::cout << std::setprecision(16);
                std::cout << "intGamma = " << intGamma << "\n";
                std::cout << "length error = " << fabs(intGamma - 2 * 0.17 * pi) << "\n";

                //errL2 = L2_norm_surface(funuh_0, fun_sol_surfactant, *interface(0), tid, phi, 0, 1);
                errL2 = L2_norm_surface(funuh_0.expr(), fun_sol_surfactant, *interface(0), In, qTime, 0, phi);
                std::cout << " t_n -> || u-uex||_2 = " << errL2 << "\n";
                //errL2 = L2_norm_surface(funuh, fun_sol_surfactant, *interface(lastQuadTime), tid + dT, phi, 0, 1);
                
                errL2 = L2_norm_surface(funuh.expr(), fun_sol_surfactant, *interface(lastQuadTime), In, qTime,
                lastQuadTime, phi);
                std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << "\n";

                // Conservation error
                double q0 = integral_algoim<Levelset<2>, Fun_h>(funuh_0, *interface(0), 0, phi, In, qTime, 0);
                double q1 = integral_algoim<Levelset<2>, Fun_h>(funuh, *interface(lastQuadTime), 0, phi, In, qTime,
                                                                lastQuadTime);
                if (iter == 0) {
                    q_init0 = q0;
                    // q_init1 = q1;
                    qp1     = q1;
                    // q_init1 = integral_algoim<Levelset<2>, Fun_h>(u0, *interface(0), 0, phi);
                }

                output_data << std::setprecision(10);
                output_data << tid << "," << (q1 - qp1) << "," << intF << "," << ((q1 - qp1) - intF) << "\n";
                qp1 = q1;


                error_t.push_back(errL2);
                numb_stabilized_edges.push_back(surfactant.get_number_of_stabilized_edges());
                numb_of_active_elements.push_back(ThGamma.get_nb_element());
                errors.at(j)         = errL2;
                gamma_length_h.at(j) = intGamma;

                if (iterations == 1) {
                    Fun_h sol(Wh, datau0);

                    Paraview<Mesh> ThB(Th, path_figures + "Th.vtk");
                    Paraview<Mesh> writer(ThGamma, path_figures + "surfactant" + std::to_string(iter + 1) + ".vtk");

                    Fun_h uS_ex(Wh, fun_sol_surfactant, tid);
                    writer.add(u0, "surfactant", 0, 1);
                    writer.add(uS_ex, "surfactant_exact", 0, 1);
                    writer.add(fabs(u0.expr() - uS_ex.expr()), "surfactant_error");
                    writer.writeActiveMesh(ThGamma, path_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                    writer.writeFaceStab(ThGamma, 0, path_figures + "Edges" + std::to_string(iter + 1) + ".vtk");
                    writer.writeAlgoimQuadrature(ThGamma, phi, In, qTime, 0, 2,
                                                 path_figures + "algoim_quadrature_0_" + std::to_string(iter + 1) +
                                                     ".vtk");
                    writer.writeAlgoimQuadrature(ThGamma, phi, In, qTime, 1, 2,
                                                 path_figures + "algoim_quadrature_1_" + std::to_string(iter + 1) +
                                                     ".vtk");
                    writer.writeAlgoimQuadrature(ThGamma, phi, In, qTime, 2, 2,
                                                 path_figures + "algoim_quadrature_2_" + std::to_string(iter + 1) +
                                                     ".vtk");
                }

                if (iterations > 1 && iter == total_number_iteration - 1)
                    output_data << h << "," << dT << "," << errL2 << "\n";

                iter++;

            }


            std::cout << "\n";
            std::cout << "Error t = [";
            for (auto &err : error_t) {

                std::cout << err;

                std::cout << ", ";
            }
            std::cout << "]\n";
            std::cout << "Elements in the active mesh = [";
            for (auto &elems : numb_of_active_elements) {

                std::cout << elems;

                std::cout << ", ";
            }
            std::cout << "]\n";
            std::cout << "Number of stabilized edges = [";
            for (auto &numb : numb_stabilized_edges) {

                std::cout << numb;

                std::cout << ", ";
            }
            std::cout << "]\n\n";

            average_errors.at(j) = std::reduce(error_t.begin(), error_t.end()) / static_cast<float>(error_t.size());

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

        std::cout << std::setprecision(16);
        std::cout << "\n";
        std::cout << "Errors = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << errors.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "Average errors = [";
        for (int i = 0; i < iterations; i++) {
            std::cout << average_errors.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "h = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << hs.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "nx = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << nxs.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "ny = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << nys.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "dT = [";
        for (int i = 0; i < iterations; i++) {

            std::cout << dts.at(i);
            if (i < iterations - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "error in length of Gamma(t) vs h = [";
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
