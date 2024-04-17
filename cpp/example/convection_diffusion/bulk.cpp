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
    Find u in Omega(t) such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    in Omega(t).
                            n*D*grad(u) = 0,    on Gamma(t).

 *  Numerical method:
    A space-time Cutfem.

 *  Two schemes are implemented:
 *  Classical scheme,
 *  Conservative scheme.
 
 *  Two types of stabilization is possible: face-based ghost penalty and patch-based ghost penalty.
 *  Moreover, each of them can be used with either full stabilization or macroelement stabilization.
*/

// Dependencies
#include "../tool.hpp"
#include "../problem/AlgoimIntegration.hpp"
#include "../num/matlab.hpp"
#include <string>
using namespace globalVariable; // to access some globally defined constants

// Numerical examples
namespace Circle {
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

} // namespace Circle

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

using mesh_t        = MeshQuad2;
using fun_test_t    = TestFunction<mesh_t>;
using fct_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using fespace_t     = GFESpace<mesh_t>;
using cut_fespace_t = CutFESpace<mesh_t>;

std::string example, method, stabilization;

// Define example, method and stabilization

#define kite     // circle/kite
#define conservative   // classical/conservative
#define macro       // fullstab/macro


int main(int argc, char **argv) {

// Do not touch
#if defined(circle)
    example = "circle";
    using namespace Circle;
#elif defined(kite)
    example = "kite";
    using namespace Kite;
#else
    #error "No example defined"
#endif

#if defined(classical)
    method = "classical";
#elif defined(conservative)
    method = "conservative";
#else
    #error "No method defined"
#endif

#if defined(fullstab)
    stabilization = "fullstab";
#elif defined(macro)
    stabilization = "macro";
#else
    #error "No stabilization defined"
#endif

#ifdef USE_MPI
    MPIcf cfMPI(argc, argv);
#endif
                                            
    const  size_t iterations = 1;           // number of mesh refinements   
    double h                = 0.1;          // starting mesh size
    int    nx, ny;                          

    const double cfl_number = 5./18;
    int total_number_iteration;
    double dT = 0.1;
    const double t0 = 0., tfinal = .5;

    // Time integration quadrature
    const size_t quadrature_order_time = 5;
    const QuadratureFormular1d &qTime(*Lobatto(quadrature_order_time));  // specify order of quadrature in time
    const Uint nbTime       = qTime.n;
    const Uint lastQuadTime = nbTime - 1;

    // Space integration quadrature 
    Levelset<2> phi;
    ProblemOption option;
#ifdef USE_MPI
    option.solver_name_ = "mumps";
#else
    option.solver_name_ = "umfpack";
#endif
    const int quadrature_order_space       = 5;
    option.order_space_element_quadrature_ = quadrature_order_space;
    AlgoimCutFEM<mesh_t, Levelset<2>> convdiff(qTime, phi, option);

    // Method parameters
    const double tau = 1.;     // stabilization constant
    const double delta = 0.5;   // macro parameter
    
    // Paths to store data
    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/paper/" + example + "/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/paper/" + example + "/paraview/";

    // Data file to hold problem data
    std::ofstream output_data(path_output_data + "data_" + method + "_" + stabilization + ".dat", std::ofstream::out);
    output_data << "example, method, stabilization, delta, tau, N_time, N_space, T\n";
    output_data << example << ", " << method << ", " << stabilization << ", ";
    #if defined(macro)
        output_data << delta << ", ";
    #else
        output_data << "0, ";
    #endif
    output_data << tau << ", " << quadrature_order_time << ", " << quadrature_order_space << ", " << tfinal;
    output_data << "\n---------------------\n";
    output_data << "h, dt, L2(Omega(T)), L2(L2(Omega(t), 0, T)), e_c(T)\n";
    output_data.flush();

    // Arrays to hold data
    std::array<double, iterations> hs, dts, errors, errors_T, global_conservation_errors;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {       

        // Define background mesh
        #if defined(circle)
        const double x0 = 0. - Epsilon, y0 = 0. - Epsilon;
        const double lx = 1., ly = 1.;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;

        #elif defined(kite)
        const double x0 = -3.5-Epsilon, y0 = -1.5 - Epsilon;
        const double lx = 7., ly = 3.;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;

        #endif

        mesh_t Th(nx, ny, x0, y0, lx, ly);

        dT = cfl_number*h;
        total_number_iteration = int(tfinal / dT);
        dT        = tfinal / total_number_iteration;

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
        std::cout << "dT = " << dT << '\n';

        fespace_t Vh(Th, DataFE<mesh_t>::P2);       // Background FE Space

        // 1D Time mesh
        double final_time = total_number_iteration * dT;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P2Poly); // FE Space in time
        const Uint ndof_time_slab = Ih[0].NbDoF();

        // Velocity field
        LagrangeQuad2 FEvelocity(2);
        fespace_t VelVh(Th, FEvelocity);
        std::vector<fct_t> vel(nbTime);

        // Declare time dependent interface
        TimeInterface<mesh_t> interface(qTime);

        // Visualization of the level set
        fespace_t Lh(Th, DataFE<mesh_t>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<fct_t> ls(nbTime);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        double mass_last_previous = 0., mass_initial = 0., mass_last = 0.;
        double intF = 0, intG = 0, intF_total = 0, intG_total = 0;
        double global_conservation_error = 0, errBulk = 0., error_I = 0.;
        std::vector<double> global_conservation_errors_t;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * dT;
            const TimeSlab &In(Ih[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * dT << '\n';
            std::cout << "dT = " << dT << '\n';

            // Initialization of the interface in each quadrature point
            for (int i = 0; i < nbTime; ++i) {

                R tt  = In.Pt(R1(qTime(i).x));
                phi.t = tt;

                vel[i].init(VelVh, fun_velocity, tt);

                interface.init(i, Th, phi);

                ls[i].init(Lh, fun_levelSet, tt);

            }

            // Create active meshes
            activemesh_t Thi(Th);

            Thi.truncate(interface, 1);

            //  Cut FE space
            cut_fespace_t Wh(Thi, Vh);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;

            // Test and Trial functions
            fun_test_t u(Wh, 1), v(Wh, 1);

            std::vector<double> data_init(convdiff.get_nb_dof(), 0.); // initial data total
            std::span<double> data_init_span(data_init);
            std::span<double> data_uh0 = std::span<double>(data_init.data(), Wh.NbDoF());

            if (iter == 0) {
                interpolate(Wh, data_uh0, fun_uBulkInit);
            } else {
                convdiff.initialSolution(data_init_span);
            }

            std::vector<double> data_all(data_init);
            std::span<double> data_uh = std::span<double>(data_all.data(), Wh.NbDoF()*In.NbDoF());

            fct_t uh0(Wh, data_uh0);
            fct_t uh(Wh, In, data_uh);

            // Variational formulation
            #if defined(classical)
                convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);
                // Impose initial condition
                if (iter == 0) {
                    convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
                } else {
                    convdiff.addLinear(+innerProduct(uh0.expr(), v), Thi, 0, In);
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
                    convdiff.addLinear(+innerProduct(uh0.expr(), v), Thi, 0, In);
                }

                convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)), Thi, In);

                for (int i = 0; i < nbTime; ++i) {
                    convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Thi, In, i);
                }
            #endif

            // Source function
            convdiff.addLinearExact(fun_rhsBulk, +innerProduct(1, v), Thi, In);

            // Stabilization

            #if defined(fullstab)

                // convdiff.addFaceStabilization(
                //     +innerProduct(h * tau * jump(grad(u) * n), jump(grad(v) * n)) +
                //         innerProduct(h * h * h * tau * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
                //     Thi, In);

                convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(u), jump(v)), Thi, In);
            

            #elif defined(macro)

                AlgoimMacro<mesh_t, Levelset<2>> TimeMacro(Thi, delta, phi, In, qTime);
                TimeMacro.findSmallElement();
                TimeMacro.createMacroElement();
                TimeMacro.setInnerEdges();

                // convdiff.addFaceStabilization(
                //     +innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)) +
                //         innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
                //     Thi, In, TimeMacro);

                convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(u), jump(v)), Thi, In, TimeMacro);

                if (iterations == 1 && h > 0.1) {
                    Paraview<mesh_t> writerMacro(Th, path_output_figures + "Th" + std::to_string(iter + 1) + ".vtk");
                    writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                    writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                    writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

                    // domain = 0

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
                matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + "_" + method + "_" + stabilization  +  ".dat");
            }

            // Solve linear system
            convdiff.solve();


            std::span<double> rhs = std::span<double>(convdiff.rhs_.data(), convdiff.get_nb_dof());
            data_all.assign(rhs.begin(), rhs.end());
            std::span<double> data_all_span(data_all);
            convdiff.saveSolution(data_all_span);

            // Compute (int_In ||u-uex||^2 dt)^2
            fct_t uh_t(Wh, In, data_all); 
            error_I += L2_norm_T(uh_t, fun_uBulkD, Thi, In, qTime, phi, quadrature_order_space);
            std::cout << " t_n -> || u-uex||_(In x Omega)^2 = " << error_I << '\n';

            // Compute error of numerical solution
            std::vector<double> sol(Wh.get_nb_dof());

            for (int n = 0; n < ndof_time_slab; n++) {
                std::vector<double> u_dof_n(data_uh.begin() + n * Wh.get_nb_dof(), data_uh.begin() + (n + 1) * Wh.get_nb_dof());
                std::transform(sol.begin(), sol.end(), u_dof_n.begin(), sol.begin(), std::plus<double>()); // sum up all dofs
            }

            fct_t funuh(Wh, sol);

            errBulk = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1, quadrature_order_space);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            errors[j] = errBulk;

            // Compute conservation error
            intF = integral_algoim(fun_rhsBulk, 0, Thi, phi, In, qTime,
                                   quadrature_order_space); // integrate source over In
            intG = integral_algoim(fun_neumann_Gamma, In, interface, phi, 0,
                                   quadrature_order_space); // integrate flux boundary over In

            intF_total += intF;
            intG_total += intG;

            mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime, quadrature_order_space);

            if (iter == 0) {
                mass_initial       = integral_algoim(fun_uBulk, Thi, phi, In, qTime, 0, quadrature_order_space);
                mass_last_previous = mass_initial;
            }

            global_conservation_error = (mass_last - mass_initial - intF_total - intG_total);

            std::cout << "global_conservation_error: " << global_conservation_error << "\n";

            mass_last_previous = mass_last; // set current last to previous last for next time slab
                                            //}

            global_conservation_errors[j] = std::fabs(global_conservation_error);
            global_conservation_errors_t.push_back(std::fabs(global_conservation_error));

            // Write paraview files
            if ((iterations == 1) && (h > 0.01)) {
                Paraview<mesh_t> writerTh(Th, path_output_figures + "Th.vtk");
                Paraview<mesh_t> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) + ".vtk");
                writer.add(uh0, "bulk_0", 0, 1);
                writer.add(funuh, "bulk_N", 0, 1);

                fct_t uBex(Wh, fun_uBulk, current_time);
                fct_t uBex_N(Wh, fun_uBulk, current_time+dT);
                fct_t fB(Wh, fun_rhsBulk, current_time);

                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                writer.add(fabs(uh0.expr() - uBex.expr()), "bulk_error_0");
                writer.add(fabs(funuh.expr() - uBex_N.expr()), "bulk_error_N");
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);
                writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
                writer.writeFaceStab(Thi, 0, path_output_figures + "Edges" + std::to_string(iter + 1) + ".vtk");

                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, 0, 2,
                                             path_output_figures + "AlgoimQuadrature_0_" + std::to_string(iter + 1) +
                                                 ".vtk");

                writer.writeAlgoimQuadrature(Thi, phi, In, qTime, lastQuadTime, 2,
                                             path_output_figures + "AlgoimQuadrature_N_" + std::to_string(iter + 1) +
                                                 ".vtk");
            }

            iter++;

        }

        errors_T[j] = std::sqrt(error_I);

        if (iterations > 1) {
            output_data << h << "," << dT << "," << errBulk << "," << errors_T[j] << "," << global_conservation_errors[j] << '\n';
            output_data.flush();
        }

        std::cout << "error_T = " << errors_T[j] << "\n";

        if (iterations == 1) {
            std::cout << "\n";
            std::cout << "Global conservation errors = [";
            for (auto &err : global_conservation_errors_t) {

                std::cout << err;

                std::cout << ", ";
            }
            std::cout << "]\n";
        }

        h *= 0.5;
        //dT *= 0.5;

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

    return 0;
}
