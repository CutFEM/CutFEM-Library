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
const double D        = 0.01;
const double R0       = 0.17;
const double beta_max = M_PI * 0.5;

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
using namespace Example2;
#elif defined(ex3)
using namespace Example3;
#endif

int main(int argc, char **argv) {
    // MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.1;             // starting mesh size
    double dT = 0.1;

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

#endif

    // Create directory if not already existent
    // // if (MPIcf::IamMaster()) {
    // std::filesystem::create_directories(path_output_data);
    // std::filesystem::create_directories(path_output_figures);
    //}

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
    std::array<double, iterations> global_conservation_errors;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
#if defined(ex1)
        const double lx = 1., ly = 1.;
        // const double lx = 2., ly = 1.1;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0. - Epsilon, 0. - Epsilon, lx, ly);
#elif defined(ex2)
        const double lx = 6. * Example2::R0, ly = 6. * Example2::R0;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3.01 * Example2::R0, -3.01 * Example2::R0, lx, ly);
        // Mesh Th(12, 11, -.85, -0.8, 2.05, 1.6);
#elif defined(ex3)
        const double lx = 4. * Example3::R0, ly = 3. * Example3::R0;
        // const double lx = 2., ly = 1.1;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -2.01 * Example3::R0, -1.5 * Example3::R0, lx, ly);
        // Mesh Th(nx, ny, -1., -0.6, lx, ly);
#endif

        // Parameters
        // const double tfinal = 2.; // Final time
        const double tfinal = .1; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        const int divisionMeshSize = 3;

        double dT = h / divisionMeshSize;
        // double dT = tfinal / 6;

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

        // const double D = N_sphere::D;

        // CG stabilization parameter
        // const double tau1 = 0.1 * (D + beta_max), tau2 = 0.1 * (D + beta_max);
        const double tau1 = 0.1, tau2 = 0.1;

        FESpace Vh(Th, DataFE<Mesh>::P2);               // Background FE Space
        FESpace Vh_interpolation(Th, DataFE<Mesh>::P3); // for interpolating data

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly);               // FE Space in time
        FESpace1 Ih_interpolation(Qh, DataFE<Mesh1>::P3Poly); // for interpolating data

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
        option.order_space_element_quadrature_ = 5;
        AlgoimCutFEM<Mesh, Levelset<2>> convdiff(qTime, phi, option);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous;
        double intF = 0, int_outflow = 0, intG = 0; // hold integrals of rhs and Neumann bcs
        double global_conservation_error = 0;
        double errBulk                   = 0.;

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

            // Right hand side functions
            Fun_h f(Vh_interpolation, In_interpolation, fun_rhsBulk);
            Fun_h g_Neumann(Vh_interpolation, In_interpolation, fun_neumann_Gamma); // computer Neumann BC

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
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            convdiff.addBilinear(+innerProduct(dt(u), v) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(+innerProduct((vel[i].exprList() * grad(u)), v), Thi, In, i);
            }

#elif defined(conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);
            convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Thi, In, i);
            }
#endif

            // Source function
            convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);

            // Neumann boundary condition
            convdiff.addLinear(+innerProduct(g_Neumann.expr(), v), interface, In);

            // Stabilization
            double stab_bulk_faces = 0.;
            double stab_mass       = 0.; // tau1 * h;
            double stab_dt         = 0.; // tau1 * h;

#if defined(fullstab)

            // convdiff.addFaceStabilization(
            //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
            //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In);
            
            convdiff.addPatchStabilization(
                +innerProduct(tau1/h/h * jump(u), jump(v))
                , Thi, In);

#elif defined(macro)
            // MacroElementPartition<Mesh> TimeMacro(Thi, 0.3);
            // std::cout << TimeMacro.number_of_stabilized_edges << "\n";
            // std::cout << "number of stabilized edges: " << convdiff.get_number_of_stabilized_edges() << "\n";

            AlgoimMacro<Mesh, Levelset<2>> TimeMacro(Thi, 0.5, phi, In, qTime);

            // convdiff.addFaceStabilization(
            //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
            //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     Thi, In, TimeMacro);

            convdiff.addPatchStabilization(
                +innerProduct(tau1/h/h * jump(u), jump(v))
                , Thi, In, TimeMacro);

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

            // Compute area of domain in time quadrature point 0
            Fun_h funone(Wh, fun_one);

            double intGamma = integral_algoim(funone, In, interface, phi, 0) / dT;
            double intOmega = integral_algoim(funone, Thi, phi, In, qTime, lastQuadTime);
            gamma[j]        = intGamma;
            omega[j]        = intOmega;

            // Compute error of numerical solution
            Rn sol(Wh.get_nb_dof(), 0.);
            Fun_h funuh(Wh, sol);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));

            errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, 0, phi, 0, 1);
            std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << '\n';

            for (int n = 1; n < ndfTime; n++) {
                sol += data_u0(SubArray(Wh.get_nb_dof(), n * Wh.get_nb_dof()));
            }

            errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';
            errors[j] = errBulk;

            // Compute conservation error
            // if (iterations == 1 && h > 0.01) {
            intF = integral_algoim(f, 0, Thi, phi, In, qTime);        // integrate source over In
            intG = integral_algoim(g_Neumann, In, interface, phi, 0); // integrate flux across boundary over In

            double mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime); // mass in last quad point

            if (iter == 0) {
                mass_last_previous = integral_algoim(b0h, Thi, phi, In, qTime, 0);
            } else {
                global_conservation_error += ((mass_last - mass_last_previous) - intF - intG);
            }
            std::cout << "global_conservation_error: " << global_conservation_error << "\n";

            outputData << std::setprecision(10);
            outputData << current_time << "," << (mass_last - mass_last_previous) << "," << intF << "," << intG << ","
                       << ((mass_last - mass_last_previous) - intF - intG) << '\n';

            mass_last_previous = mass_last; // set current last to previous last for next time slab
            //}

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

            global_conservation_errors[j] = std::fabs(global_conservation_error);

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
