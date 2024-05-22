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
 *  Non-conservative scheme,
 *  Conservative scheme.

 *  Two types of stabilization is possible: face-based ghost penalty and patch-based ghost penalty.
 *  Moreover, each of them can be used with either full stabilization or macroelement stabilization.
*/

// Dependencies
#include "../cutfem.hpp"
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
    return D * ((M_PI *
                 sin((M_PI *
                      sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * M_PI) * 1.0 /
                 sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                      pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0) *
                 2.0E+8) /
                    R0 +
                (1.0 / (R0 * R0) * (M_PI * M_PI) * sin(t * M_PI) *
                 cos((M_PI *
                      sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12) /
                    (pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                     pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0) -
                (M_PI *
                 sin((M_PI *
                      sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * M_PI) * pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) *
                 1.0 /
                 pow(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0,
                     3.0 / 2.0) *
                 4.0E+20) /
                    R0 -
                (M_PI *
                 sin((M_PI *
                      sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 sin(t * M_PI) * pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) *
                 1.0 /
                 pow(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0,
                     3.0 / 2.0) *
                 4.0E+20) /
                    R0 +
                (1.0 / (R0 * R0) * (M_PI * M_PI) * sin(t * M_PI) *
                 cos((M_PI *
                      sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                           pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
                     (R0 * 1.0E+8)) *
                 pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12) /
                    (pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                     pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) +
           M_PI * cos(t * M_PI) *
               cos((M_PI *
                    sqrt(pow(y * 5.0E+1 + cos(t * M_PI) * 1.4E+1 - 2.5E+1, 2.0) * 4.0E+12 +
                         pow(x * -5.0E+1 + sin(t * M_PI) * 1.4E+1 + 2.5E+1, 2.0) * 4.0E+12 + 1.0)) /
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

    return D * ((M_PI *
                 sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * 1.0 / sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) /
                    R0 +
                (1.0 / (R0 * R0) * (M_PI * M_PI) *
                 cos((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * pow(x * 2.0 - t * (W - C * (y * y)) * 2.0, 2.0)) /
                    (pow(x - t * (W - C * (y * y)), 2.0) * 4.0 + (y * y) * 4.0 + 4.0E-16) +
                (M_PI *
                 sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * 1.0 / sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16) *
                 ((C * C) * (t * t) * (y * y) * 8.0 + C * t * (x - t * (W - C * (y * y))) * 4.0 + 2.0)) /
                    (R0 * 2.0) -
                (M_PI *
                 sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * pow(y * 2.0 + C * t * y * (x - t * (W - C * (y * y))) * 4.0, 2.0) * 1.0 /
                 pow(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16, 3.0 / 2.0)) /
                    (R0 * 4.0) +
                (1.0 / (R0 * R0) * (M_PI * M_PI) *
                 cos((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * pow(y * 2.0 + C * t * y * (x - t * (W - C * (y * y))) * 4.0, 2.0)) /
                    (pow(x - t * (W - C * (y * y)), 2.0) * 4.0 + (y * y) * 4.0 + 4.0E-16) -
                (M_PI *
                 sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
                 sin(t * M_PI) * pow(x * 2.0 - t * (W - C * (y * y)) * 2.0, 2.0) * 1.0 /
                 pow(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16, 3.0 / 2.0)) /
                    (R0 * 4.0)) +
           M_PI *
               cos((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
               cos(t * M_PI) +
           (M_PI *
            sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
            sin(t * M_PI) * (W - C * (y * y)) * (x - t * (W - C * (y * y))) * 1.0 /
            sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) /
               R0 -
           (M_PI *
            sin((M_PI * sqrt(pow(x - t * (W - C * (y * y)), 2.0) + y * y + 1.0E-16)) / R0) *
            sin(t * M_PI) * (W - C * (y * y)) * (x * 2.0 - t * (W - C * (y * y)) * 2.0) * 1.0 /
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

std::vector<const GTypeOfFE<mesh_t> *> FE_space = {&DataFE<mesh_t>::P1, &DataFE<mesh_t>::P2, &DataFE<mesh_t>::P3};
std::vector<const GTypeOfFE<Mesh1> *> FE_time   = {&DataFE<Mesh1>::P0Poly, &DataFE<Mesh1>::P1Poly, &DataFE<Mesh1>::P2Poly,
                                                   &DataFE<Mesh1>::P3Poly};


// Define example, method and stabilization, and polynomial order in time and space

#define circle         // example (circle/kite)
#define conservative // method (conservative/non_conservative)
#define macro        // stabilization (fullstab/macro)
#define K 2          // polynomial order in time (0/1/2/3)
#define M 2          // polynomial order in space (1/2/3)


int main(int argc, char **argv) {

    std::string example, method, stabilization;

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

#if defined(non_conservative)
    method = "non_conservative";
#elif defined(conservative)
    method              = "conservative";
#else
#error "No method defined"
#endif

#if defined(fullstab)
    stabilization = "fullstab";
#elif defined(macro)
    stabilization       = "macro";
#else
#error "No stabilization defined"
#endif

    //MPIcf cfMPI(argc, argv);

    const int k = K;
    const int m = M;

    assert(k >= 0 && k <= 3);
    assert(m >= 1 && m <= 3);

    const size_t iterations = 3;   // number of mesh refinements
    double h                = 0.1; // starting mesh size
    int nx, ny;

    const double cfl_number = 1./3;
    int total_number_iteration;
    double dT       = 0.1;
    const double t0 = 0., tfinal = .1;

    // Time integration quadrature
    const size_t quadrature_order_time = 5;
    const QuadratureFormular1d &qTime(*Lobatto(quadrature_order_time)); // specify order of quadrature in time
    const Uint nbTime       = qTime.n;
    const Uint lastQuadTime = nbTime - 1;

    // Space integration quadrature
    Levelset<2> phi;
    ProblemOption option;
#ifdef USE_MPI
    option.solver_name_ = "mumps";
#else
    //option.solver_name_ = "umfpack";
    option.solver_name_ = "mumps";
#endif
    const int quadrature_order_space       = 5;
    option.order_space_element_quadrature_ = quadrature_order_space;
    AlgoimCutFEM<mesh_t, Levelset<2>> convdiff(qTime, phi, option);

    // Method parameters
    const double tau   = 1.;  // stabilization constant
    const double delta = 0.5; // macro parameter

    assert(tau > 0.);
    assert(delta > 0. && delta <= 1.);
    
    // Data file to hold problem data
    std::ofstream output_data("../cpp/example/convection_diffusion/output_bulk.dat", std::ofstream::out);
    output_data << "Example: " << example << "\n";
    output_data << "Method: " << method << "\n";
    output_data << "Polynomial order time: " << k << "\n";
    output_data << "Polynomial order space: " << m << "\n";
    output_data << "Stabilization: " << stabilization << "\n";
    if (stabilization == "macro")
        output_data << "Macroelement parameter: " << delta << "\n";
    output_data << "Stabilization constant: " << tau << "\n";
    output_data << "Quadrature order time: " << quadrature_order_time << "\n";
    output_data << "Quadrature order space: " << quadrature_order_space << "\n";
    output_data << "Tfinal: " << tfinal << "\n\n";
    output_data << "h" << ",\t\t\t\t\t\t" <<  "dt" << ",\t\t\t\t\t\t" << "L2(Omega(T))" << ",\t\t\t" << "H1(Omega(T))" << ",\t\t\t" << "L2(L2(Omega(t), 0, T))" << ",\t" << "L2(H1(Omega(t), 0, T))" << ",\t" << "e_c(T)\n";

    output_data.flush();

    // Arrays to hold data
    std::array<double, iterations> hs, dts, L2_errors, H1_errors, L2L2_errors, L2H1_errors, global_conservation_errors;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

// Define background mesh
#if defined(circle)
        const double x0 = 0. - Epsilon, y0 = 0. - Epsilon;
        const double lx = 1., ly = 1.;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;

#elif defined(kite)
        const double x0 = -3.5 - Epsilon, y0 = -1.5 - Epsilon;
        const double lx = 7., ly = 3.;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;

#endif

        mesh_t Th(nx, ny, x0, y0, lx, ly);

        dT                     = cfl_number * h;
        total_number_iteration = int(tfinal / dT);
        dT                     = tfinal / total_number_iteration;

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

        fespace_t Vh(Th, *FE_space[m-1]);
        // 1D Time mesh
        double final_time = total_number_iteration * dT;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);

        // 1D Time space
        FESpace1 Ih(Qh, *FE_time[k]);
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

        int iter                  = 0;
        double mass_last_previous = 0., mass_initial = 0., mass_last = 0.;
        double intF = 0, intG = 0, intF_total = 0, intG_total = 0;
        double global_conservation_error = 0, L2_error = 0., H1_error = 0., L2L2_error = 0., L2H1_error = 0.;
        std::vector<double> global_conservation_errors_t, L2_errors_t, H1_errors_t;

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
            std::span<double> data_uh0 = std::span<double>(data_init.data(), Wh.NbDoF());

            if (iter == 0) {
                interpolate(Wh, data_uh0, fun_uBulkInit);
            } else {
                convdiff.initialSolution(data_init);
            }

            std::vector<double> data_all(data_init);
            std::span<double> data_uh = std::span<double>(data_all.data(), Wh.NbDoF() * In.NbDoF());

            fct_t uh0(Wh, data_uh0);
            fct_t uh(Wh, In, data_uh);

// Variational formulation
#if defined(non_conservative)
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

#endif

            // if (iter == total_number_iteration - 1) {
            //     matlab::Export(convdiff.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + "_" + method +
            //                                          "_" + stabilization + ".dat");
            // }

            // Solve linear system
            convdiff.solve();

            std::span<double> rhs = std::span<double>(convdiff.rhs_.data(), convdiff.get_nb_dof());
            data_all.assign(rhs.begin(), rhs.end());
            convdiff.saveSolution(data_all);

            // Compute (int_In ||u-uex||^2 dt)^2
            fct_t uh_t(Wh, In, data_all);
            fct_t u_t(Wh, In, fun_uBulk);

            L2L2_error += L2L2_norm(uh_t, fun_uBulkD, Thi, In, qTime, phi, quadrature_order_space);
            L2H1_error += L2H1_norm(uh_t, u_t, Thi, In, qTime, phi, quadrature_order_space);
            std::cout << " t_n -> || u-uex||_(L2L2)^2 = " << L2L2_error << '\n';
            std::cout << " t_n -> || u-uex||_(L2H1)^2 = " << L2H1_error << '\n';

            // Compute error of numerical solution
            std::vector<double> sol(Wh.get_nb_dof());

            for (int n = 0; n < ndof_time_slab; n++) {
                std::vector<double> u_dof_n(data_uh.begin() + n * Wh.get_nb_dof(),
                                            data_uh.begin() + (n + 1) * Wh.get_nb_dof());
                std::transform(sol.begin(), sol.end(), u_dof_n.begin(), sol.begin(),
                               std::plus<double>()); // sum up all dofs
            }

            fct_t funuh(Wh, sol);

            auto uhdx = dx(funuh.expr());
            auto uhdy = dy(funuh.expr());

            fct_t u_exact(Wh, fun_uBulkD, current_time + dT);
            auto udx = dx(u_exact.expr());
            auto udy = dy(u_exact.expr());

            L2_error = L2_norm_cut(funuh, fun_uBulkD, In, qTime, lastQuadTime, phi, 0, 1, quadrature_order_space);
            H1_error = std::sqrt(integral_algoim((funuh.expr() - u_exact.expr())*(funuh.expr() - u_exact.expr()) + (uhdx - udx)*(uhdx - udx) + (uhdy-udy)*(uhdy-udy), Thi, phi, In, qTime, lastQuadTime, quadrature_order_space));

            L2_errors[j] = L2_error;
            H1_errors[j] = H1_error;
            L2_errors_t.push_back(L2_error);            
            H1_errors_t.push_back(H1_error);

            std::cout << " t_n -> || u-uex||_L2 = " << L2_error << '\n';
            std::cout << " t_n -> || u-uex||_H1 = " << H1_error << '\n';

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

            global_conservation_errors[j] = std::fabs(global_conservation_error);
            global_conservation_errors_t.push_back(std::fabs(global_conservation_error));

            iter++;
        }

        L2L2_errors[j] = std::sqrt(L2L2_error);
        L2H1_errors[j] = std::sqrt(L2H1_error);

        output_data << std::fixed << std::setprecision(16) << h << ",\t\t" << dT << ",\t\t" << L2_error << ",\t\t" << H1_error << ",\t\t" << L2L2_errors[j] << ",\t\t" << L2H1_errors[j] << ",\t\t" << std::scientific << global_conservation_errors[j] << '\n';
        output_data.flush();
        
        std::cout << "error_L2L2 = " << L2L2_errors[j] << "\n";
        std::cout << "error_L2H1 = " << L2H1_errors[j] << "\n";

        if (iterations == 1) {
            std::cout << "\n";
            std::cout << "L2 errors = [";
            for (auto &err : L2_errors_t) {

                std::cout << err;

                std::cout << ", ";
            }
            std::cout << "]\n";

            std::cout << "\n";
            std::cout << "H1 errors = [";
            for (auto &err : H1_errors_t) {

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
        }

        h *= 0.5;
        // dT *= 0.5;
    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "L2 Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << L2_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "H1 Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << H1_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "L2L2 Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << L2L2_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "L2H1 Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << L2H1_errors.at(i);
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
