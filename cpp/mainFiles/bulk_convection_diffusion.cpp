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
 * @note We consider a time-dependent bulk problem on Omega2.
   
 *  Problem:
    Find u in Omega_2(t) such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    in Omega_2(t).
                                        + BCS,  on Gamma(t).

 *  Numerical method:
    A space-time Cutfem, using the level-set method,
    which allows for both dg and cg.

 *  Classical scheme: Integration by parts on convection term if dg, 
    otherwise just integration by parts on diffusion term.
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

using namespace globalVariable;     // to access some globally defined constants

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
//     R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
//     return (P[0] - xc) / (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
// }

// // Normal y-direction
// R n2(double *P, const R t) {
//     R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
//     return (P[1] - yc) / (sqrt((P[1] - yc) * (P[1] - yc) + (P[0] - xc) * (P[0] - xc)));
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
} // namespace Example1

namespace Example1_Omega1 {

// Level-set function
double fun_levelSet(double *P, const int i, const R t) {
    R xc = 0.5 + 0.28 * sin(M_PI * t), yc = 0.5 - 0.28 * cos(M_PI * t);
    return sqrt((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc)) - 0.17;
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    return sqrt((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.22) * (P[1] - 0.22)) - 0.17;
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

    return -(sqrt((x - (1 - y * y) * t) * (x - (1 - y * y) * t) + y * y) - r0);
}

// Level-set function initial
double fun_levelSet(double *P, const int i) {
    double r0 = 1. + Epsilon;
    return -(sqrt(P[0] * P[0] + P[1] * P[1]) - r0);
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

// Note: standard for this program is to solve the equations on the inner domain
// Omega 2. To solve on Omega 1, define the word "omega1" for the pre-processor
// (only works for example 1)

// Choose Discontinuous or Continuous Galerkin method (options: "dg", "cg")
#define cg
// Set numerical example (options: "example1", "lehrenfeld")
#define example1
// Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominated
// Set boundary condition type on Omega 2 (options: "dirichlet", "neumann" –
// note: neumann only works combined with example1)
#define neumann
// Set scheme for the dg method (options: "classical", "conservative".
// Irrelevant if "cg" is defined instead of "dg")
#define conservative
// Set stabilization method (options: "fullstab", "macro")
#define fullstab
// Decide whether to solve for level set function, or to use exact (options:
// "levelsetsolve", "levelsetexact")
#define levelsetexact
// Solve on Omega_1 (options: "omega1" or anything else to solve on Omega 2)
#define omega2
// If "omega1" is defined, set type of BCs on outer boundary (options:
// "dirichlet1" or "neumann1")
#define dirichlet1

#define use_h       // to set mesh size using the h parameter. Write use_n to decide using nx, ny.
#define use_tnot    // write use_t to control dT manually. Otherwise it is set proportional to h.

#ifdef example1
#ifdef omega1
using namespace Example1_Omega1;
#else
using namespace Example1; // on Omega 2
#endif
#endif
#if defined(lehrenfeld)
using namespace Lehrenfeld_Convection_Dominated;
#endif

int main(int argc, char **argv) {

    // Mesh settings and data objects
    const std::size_t iterations = 5;    // number of mesh refinements   (set to 1 to run
                                    // only once and plot to paraview)
    int nx = 15, ny = 15;           // starting mesh size
    double h  = 0.1; // starting mesh size
    double dT = 0.125;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

#ifdef example1
    // Paths to store data
    const std::string pathOutputFolder  = "../output_files/bulk/example1/data/";
    const std::string pathOutputFigures = "../output_files/bulk/example1/paraview/";
#elif defined(example2)
    // Paths to store data
    const std::string pathOutputFolder  = "../output_files/bulk/example2/data/";
    const std::string pathOutputFigures = "../output_files/bulk/example2/paraview/";
#elif defined(lehrenfeld)
    // Paths to store data
    const std::string pathOutputFolder  = "../output_files/bulk/lehrenfeld/data/";
    const std::string pathOutputFigures = "../output_files/bulk/lehrenfeld/paraview/";
#endif

// Initialize MPI
#ifdef USE_MPI
    MPIcf cfMPI(argc, argv);
#endif
    // Create directory if not already existent
    //    if (MPIcf::IamMaster()) {
    //       std::experimental::filesystem::create_directories(pathOutputFolder);
    //       std::experimental::filesystem::create_directories(pathOutputFigures);
    //    }

    // Data file to hold problem data

    std::ofstream outputData(pathOutputFolder + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors;                  // array to hold bulk errors
    std::array<int, iterations> number_of_stabilized_edges; // array to count stabilized edges
    std::array<double, iterations> reynold_error; // array to count stabilized edges
    std::array<double, iterations> hs;                      // array to hold mesh sizes
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
        const int divisionMeshSize = 1;        

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

#ifdef convection_dominated
        const double A2 = 0.01;
#else
        const double A2 = 1;
#endif

#ifdef dg    
        // Constants for penalty terms
        const double tau_a2 = 500; // diffusion penalty scaling
        const double tau_b2 = 10;  // convection penalty scaling

        // Bulk penalties
        const double lambdaA = tau_a2 / h; // diffusion term
        const double lambdaB = tau_b2;     // convection term

        // Penalty parameter for outer boundary
        const double lambda = 500 / h; // only used when Dirichlet BCs apply

        // DG stabilization parameters
        const double tau20 = 1e-1, tau21 = 1e-1; // bulk
        // DG Space
        FESpace2 Vh(Th, DataFE<Mesh>::P1dc); // discontinuous basis functions
#elif defined(cg)
        // CG stabilization parameters
        const double tau20 = 0, tau21 = 0.1;
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
        const QuadratureFormular1d &qTime(*Lobatto(3));     // specify order of quadrature in time
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
        double q0_0, q0_1, qp_0, qp_1;  // integral values to be computed
        double intF = 0, intG = 0;      // hold integrals of rhs and Neumann bcs

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
            ActiveMesh<Mesh> Kh2(Th);
            Kh2.truncate(interface, 1);     // remove part with positive sign of level set to get inner domain

            // Cut FE space
            CutSpace Wh(Kh2, Vh);

            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Vh, In, fun_rhsBulk);
            Fun_h g(Vh, In, fun_uBulk); // create an FE-function of the exact bulk
                                         // solution Omega2

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
                Paraview<Mesh> writerInitial(Kh2, pathOutputFigures + "BulkInitial.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);

                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);
                Fun_h uRhs(Wh, fun_rhsBulk, 0.);

                writerInitial.add(uBex, "bulk_exact", 0, 1);
                writerInitial.add(uRhs, "bulk_rhs", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
            }

//// Assembling linear and bilinear forms

// Partial derivative in time, bulk

// reynold scheme
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
            convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Kh2, In);

            // Added terms
            convdiff.addBilinear(+innerProduct((vel * n) * average(u), jump(v))
                                     //+ innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))
                                     + innerProduct(0.5 * fabs(vel * n) * jump(u), jump(v)),
                                 Kh2, INTEGRAL_INNER_EDGE_2D, In);

#elif defined(cg) && defined(classical) // classic CG scheme
            convdiff.addBilinear(+innerProduct((vel.exprList() * grad(u)), v), Kh2, In);

#elif defined(cg) && defined(conservative) // classic CG scheme
            convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Kh2, In);

#endif

            // Stabilization

#ifdef macro

            TimeMacroElement<Mesh> TimeMacro(Kh2, qTime, 0.125);

            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, pathOutputFigures + "Th" + std::to_string(iter + 1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

                // domain = 0,

                writerMacro.writeFaceStab(Kh2, 0,
                                          pathOutputFigures + "FullStabilization" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeActiveMesh(Kh2, pathOutputFigures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0,
                                              pathOutputFigures + "macro" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroInnerEdge(
                    TimeMacro, 0, pathOutputFigures + "macro_inner_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeMacroOutterEdge(
                    TimeMacro, 0, pathOutputFigures + "macro_outer_edge" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeSmallElements(TimeMacro, 0,
                                               pathOutputFigures + "small_element" + std::to_string(iter + 1) + ".vtk");
            }

            // Stabilization of the bulk
            // convdiff.mat_.clear();
            convdiff.addFaceStabilization(+innerProduct(1. / h * tau20 * jump(u), jump(v)) +
                                              innerProduct(h * tau21 * jump(grad(u)), jump(grad(v))),
                                          Kh2, In, TimeMacro);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // getchar();

#elif defined(fullstab)

            convdiff.addFaceStabilization(+innerProduct(1. / h * tau20 * jump(u), jump(v)) +
                                              innerProduct(h * tau21 * jump(grad(u)), jump(grad(v)))
                                          ,
                                          Kh2, In);

#endif
            number_of_stabilized_edges.at(j) = convdiff.get_number_of_stabilized_edges();
            
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

            intF = integral(Kh2, In, f, 0, qTime);
#ifdef neumann
            intG = integral(g_Neumann, In, interface, 0);
#endif

            if (iter == total_number_iteration - 1)
                matlab::Export(convdiff.mat_[0],
                               pathOutputFolder + "mat_h" + std::to_string(h) + "_" + std::to_string(j + 1) + ".dat");

            // Solve linear system
            convdiff.solve("umfpack");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Compute error in Reynold relation
            {
                CutFEM<Mesh> reynold(qTime);
                reynold.initSpace(Wh, In);

                reynold.addBilinear(innerProduct(u, v), Kh2, (int)lastQuadTime, In);
                reynold.addLinear(innerProduct(b0h.expr(), v), Kh2, 0, In);
                reynold.addBilinear(-innerProduct(dt(u), v) - innerProduct(u, dt(v))
                                        - innerProduct((vel.exprList() * grad(u)), v)
                                        - innerProduct(u, (vel.exprList() * grad(v))),
                                    Kh2, In);

                int N = Wh.NbDoF();
                Rn lhs(2 * N);
                multiply(2 * N, 2 * N, reynold.mat_, data_u0, lhs);

                lhs -= reynold.rhs_;

                reynold_error.at(j) = lhs.linfty();
                std::cout << " e_r^n = " << reynold_error.at(j) << '\n';
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

                outputData << std::setprecision(10);
                outputData << current_time << "," << (q_1 - qp_1) << "," << intF << ","
#ifdef neumann
                           << intG << "," << ((q_1 - qp_1) - intF - intG) << '\n';
#elif defined(dirichlet)
                           << intGrad << "," << intVel << "," << lambda * (intu - intg) << ","
                           << h * ((q_1 - qp_1) - intF - intGrad + intVel) + lambda * (intu - intg) << '\n';
#endif
                qp_1 = q_1;
            }

            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));
            Fun_h funuh(Wh, sol);
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time, 0, 1);
            std::cout << " t_{n-1} -> || u-uex||_2 = " << errBulk << '\n';

            sol += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
            errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            errors.at(j) = errBulk;

            // #ifdef USE_MPI
            //          if ((iterations == 1) && MPIcf::IamMaster()) {
            // #else
            if ((iterations == 1)) {
                Fun_h sol(Wh, data_u0);
                Paraview<Mesh> writer(Kh2, pathOutputFigures + "Bulk" + std::to_string(iter + 1) + "DG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                Expression2 uuh(sol, 0, op_id);
                Expression2 uuex(uBex, 0, op_id);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                //writer.add(fabs(uuh - uuex), "bulk_error");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                // writer.add(ls[2], "levelSet2", 0, 1);
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
    std::cout << "Number of stabilized edges = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << number_of_stabilized_edges.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

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
