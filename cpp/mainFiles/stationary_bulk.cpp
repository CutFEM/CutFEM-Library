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

namespace ConvectionDiffusion {

// RHS for variable in Omega 1
double fun_rhs1(double *P) {
    return -exp(-pow(P[0], 2) - pow(P[1], 2) + 1) *
           (12 * pow(P[0], 4) * P[1] + 3 * pow(P[0], 3) + 8 * pow(P[0], 2) * pow(P[1], 3) - 48 * pow(P[0], 2) * P[1] -
            9 * P[0] * pow(P[1], 2) - 4 * pow(P[1], 5) + 16 * pow(P[1], 3));
}

// RHS for variable in Omega 2
double fun_rhs2(double *P) {
    return -2 * exp(-pow(P[0], 2) - pow(P[1], 2) + 1) *
           (6 * pow(P[0], 4) * P[1] + 3 * pow(P[0], 3) + 4 * pow(P[0], 2) * pow(P[1], 3) - 24 * pow(P[0], 2) * P[1] -
            9 * P[0] * pow(P[1], 2) - 2 * pow(P[1], 5) + 8 * pow(P[1], 3));
}

// RHS for bulk variables
double fun_rhsBulk(double *P, int elementComp) {
    return fun_rhs1(P);
    // return fun_rhs2(P);
}

// Exact solution in the bulk
double fun_uBulk(double *P, int elementComp, int domain) {
    return exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - pow(P[1], 3));
    // return 2.0*exp(1.0-P[0]*P[0]-P[1]*P[1])*( 3.0*P[0]*P[0]*P[1] - pow(P[1],3));
}

double fun_one(double *P, const int i) { return 1.; }

double fun_g_Neumann(double *P, int elementComp) {
    R x = P[0], y = P[1];
    // Outward pointing w.r.t Omega1
    return (y * exp(-x * x - y * y + 1) *
            (6 * x * x * x * x + 4 * x * x * y * y - 9 * x * x - 2 * y * y * y * y + 3 * y * y)) /
           sqrt((x * x + y * y));
}

// Level-set function
double fun_levelSet(double *P, const int i) { return +((P[0]) * (P[0]) + (P[1]) * (P[1]) - 1.0); }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const {
        return -(x(0) * x(0) + x(1) * x(1) - 1.);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(-2.0 * x(0), -2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(double *P) {
        R norm = sqrt(4.0 * P[0] * P[0] + 4.0 * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(-2.0 * P[0] / norm, -2.0 * P[1] / norm);
    }
};

// Velocity Field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return P[1];
    else
        return -P[0];
}
} // namespace ConvectionDiffusion

namespace ConvectionDiffusionConstVel {

// RHS for variable in Omega 1
double fun_rhsBulk(double *P, int elementComp) {
    double x = P[0], y = P[1];
    return -(exp(-x * x - y * y + 1) * (60 * pow(x, 4) * y + 12 * pow(x, 3) * y + 40 * x * x * pow(y, 3) +
                                        18 * x * x * y * y - 240 * x * x * y - 9 * x * x - 4 * x * pow(y, 3) -
                                        12 * x * y - 20 * pow(y, 5) - 6 * pow(y, 4) + 80 * pow(y, 3) + 9 * y * y)) /
           5;
}

// Exact solution in the bulk
double fun_uBulk(double *P, int elementComp, int domain) {
    return exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - pow(P[1], 3));
    // return 2.0*exp(1.0-P[0]*P[0]-P[1]*P[1])*( 3.0*P[0]*P[0]*P[1] - pow(P[1],3));
}

double fun_one(double *P, const int i) { return 1.; }

double fun_g_Neumann(double *P, int elementComp) {
    R x = P[0], y = P[1];
    // Omega1
    return (y * exp(-x * x - y * y + 1) * (6 * pow(x, 4) + 4 * x * x * y * y - 9 * x * x - 2 * pow(y, 4) + 3 * y * y)) /
           sqrt(x * x + y * y);
}

// Level-set function
double fun_levelSet(double *P, const int i) { return +((P[0]) * (P[0]) + (P[1]) * (P[1]) - 1.0); }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const {
        return -(x(0) * x(0) + x(1) * x(1) - 1.);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(-2.0 * x(0), -2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(double *P) {
        R norm = sqrt(4.0 * P[0] * P[0] + 4.0 * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(-2.0 * P[0] / norm, -2.0 * P[1] / norm);
    }
};

// Velocity Field
R fun_velocity(double *P, const int i) {
    if (i == 0)
        return 0.4;
    else
        return 0.6;
}
} // namespace ConvectionDiffusionConstVel

namespace Diffusion {

// RHS for variable in Omega 1
double fun_rhsBulk(double *P, int elementComp) {
    double x = P[0], y = P[1];
    
    return -4 * y * std::exp(-x * x - y * y + 1) *
           (3 * x * x * x * x + 2 * x * x * y * y - 12 * x * x - y * y * y * y + 4 * y * y);
}

// Exact solution in the bulk
double fun_uBulk(double *P, int elementComp, int domain) {
    return std::exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - std::pow(P[1], 3));
    // return 2.0*exp(1.0-P[0]*P[0]-P[1]*P[1])*( 3.0*P[0]*P[0]*P[1] - pow(P[1],3));
}

double fun_one(double *P, const int i) { return 1.; }

double fun_g_Neumann(double *P, int elementComp) {
    R x = P[0], y = P[1];
    // on Omega2
    // return -(y * exp(-x * x - y * y + 1) *
    //          (6 * x * x * x * x + 4 * x * x * y * y - 9 * x * x - 2 * y * y * y * y + 3 * y * y)) /
    //        std::sqrt(x * x + y * y);

    // on Omega1
    return (y * exp(-x * x - y * y + 1) * (6 * pow(x, 4) + 4 * x * x * y * y - 9 * x * x - 2 * pow(y, 4) + 3 * y * y)) /
           std::sqrt(x * x + y * y);
}

// Level-set function
double fun_levelSet(double *P, const int i) { return +((P[0]) * (P[0]) + (P[1]) * (P[1]) - 1.0); }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const {
        return -(x(0) * x(0) + x(1) * x(1) - 1.);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(-2.0 * x(0), -2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(double *P) {
        R norm = sqrt(4.0 * P[0] * P[0] + 4.0 * P[1] * P[1]);
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(-2.0 * P[0] / norm, -2.0 * P[1] / norm);
    }
};

// Velocity Field
R fun_velocity(double *P, const int i) { return 0.; }
} // namespace Diffusion

namespace Poisson {

double fun_one(double *P, const int i) { return 1.; }
double fun_rhsBulk(double *P, int elementComp) { return 40 * pow(pi, 2) * sin(2 * pi * P[0]) * sin(4 * pi * P[1]); }
double fun_uBulk(double *P, int elementComp, int domain) { return 2 * sin(2 * pi * P[0]) * sin(4 * pi * P[1]); }
double fun_g_Neumann(double *P, int elementComp) { return 0.; }
double fun_levelSet(double *P, const int i) {
    return -((P[0] - 0.5) * (P[0] - 0.5) + (P[1] - 0.5) * (P[1] - 0.5) - 0.25 * 0.25);
}
// Velocity Field
R fun_velocity(double *P, const int i) { return 0.; }
template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const {
        return ((x(0) - 0.5) * (x(0) - 0.5) + (x(1) - 0.5) * (x(1) - 0.5) - .25 * .25);
    }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * (x(0) - 0.5), 2.0 * (x(1) - 0.5));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(double *P) {
        R norm = sqrt(4.0 * (P[0] - 0.5) * (P[0] - 0.5) + 4.0 * (P[1] - 0.5) * (P[1] - 0.5));
        // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
        return R2(2.0 * (P[0] - 0.5) / norm, 2.0 * (P[1] - 0.5) / norm);
    }
};

} // namespace Poisson


using namespace ConvectionDiffusion;

#define algoim         // options: "algoim", "quad", "triangle"
#define neumann
#define use_h

// Setup two-dimensional class types
const int d = 2;
#if defined(algoim) || defined(quad)
typedef MeshQuad2 Mesh;
#else
typedef Mesh2 Mesh;
#endif
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;
typedef ExpressionFunFEM<Mesh> Expression;

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // Mesh settings and data objects
    const size_t iterations = 1; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h = 0.1;              // starting mesh size

    // Paths to store data
    const std::string path_output_data    = "../output_files/bulk/algoim/statbulk/data/";
    const std::string path_output_figures = "../output_files/bulk/algoim/statbulk/paraview/";

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    // Data file to hold problem data
    std::ofstream outputData(path_output_data + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors; // array to hold bulk errors
    std::array<double, iterations> hs;     // array to hold mesh sizes
    std::array<double, iterations> nxs;    // array to hold mesh sizes
    std::array<double, iterations> nys;    // array to hold mesh sizes
    std::array<double, iterations> omega;
    std::array<double, iterations> gamma;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Define background mesh
        const double lx = 3., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -1.5-Epsilon, -1.5-Epsilon, lx, ly);

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;

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

        // CG stabilization parameter
        const double tau1 = .01;

        const double D = 1., lambda = 1e1;

        FESpace Vh(Th, DataFE<Mesh>::P1);       // continuous basis functions

        // Velocity field
#if defined(algoim) || defined(quad)
        LagrangeQuad2 FEvelocity(1);
#elif defined(triangle)
        Lagrange2 FEvelocity(1);
#endif
        FESpace VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace Lh(Th, DataFE<Mesh>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Th, levelSet);

        // Create active meshes
        ActiveMesh<Mesh> Thi(Th);

        Thi.truncate(interface, -1); // remove part with negative sign of level

        // Cut FE space
        CutSpace Wh(Thi, Vh);

#if defined(algoim)
        // Declare algoim interface
        Levelset<2> phi;
        AlgoimCutFEM<Mesh, Levelset<2>> convdiff(Wh, phi);
#elif defined(triangle) || defined(quad)
        CutFEM<Mesh> convdiff(Wh);
#endif

        // Objects needed for the weak form
        Normal n;
        Tangent t;

        // Right hand side functions
        Fun_h f(Vh, fun_rhsBulk);
        Fun_h g(Vh, fun_uBulk); // create an FE-function of the exact solution

        Fun_h g_Neumann(Vh, fun_g_Neumann); // computer Neumann BC

        // Test and Trial functions
        FunTest u(Wh, 1), v(Wh, 1);

        // Data for initial solution
        Rn data_u0(convdiff.get_nb_dof(), 0.);                  // initial data total
        KN_<R> data_B0(convdiff.rhs_(SubArray(Wh.NbDoF(), 0))); // data bulk

        // Make function objects to use in innerProducts
        Fun_h b0h(Wh, data_B0);

        // gnuplot::save(Th);
        // gnuplot::save<Mesh, Levelset<2>>(Thi, interface, phi, 0, "interface.dat");
        // gnuplot::save<Mesh>(interface);
        // getchar();
        

    #if defined(algoim)
        convdiff.addBilinearAlgoim(+innerProduct(D * grad(u), grad(v)), Thi);
        convdiff.addBilinearAlgoim(innerProduct((vel.exprList() * grad(u)), v), Thi);
        convdiff.addLinearAlgoim(+innerProduct(f.expr(), v), Thi);
    #else
        convdiff.addBilinear(+innerProduct(D * grad(u), grad(v)), Thi);
        convdiff.addBilinear(innerProduct((vel.exprList() * grad(u)), v), Thi);
        convdiff.addLinear(+innerProduct(f.expr(), v), Thi);
    #endif

        //* Stabilization
        convdiff.addFaceStabilization(innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)), Thi);

        //* Boundary conditions

        // Inner boundary

        // Neumann (won't work if we solve on Omega2, since then the solution is not unique)
        convdiff.addLinear(innerProduct(g_Neumann.expr(), v), interface);

        // // Dirichlet
        // convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
        //                          - innerProduct(u, D * grad(v) * n) // added to make symmetric
        //                          + innerProduct(u, lambda / h * v)  // added penalty
        //                      ,
        //                      interface);

        // // RHS terms
        // convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
        //                    interface);

        // Outer boundary

        //* Nitsche's method:

        // LHS terms
        convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                 - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                 + innerProduct(u, lambda / h * v)  // added penalty
                             ,
                             Thi, INTEGRAL_BOUNDARY);

        // RHS terms
        convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v), Thi,
                           INTEGRAL_BOUNDARY);




        matlab::Export(convdiff.mat_[0], path_output_data + "mat_rank_" + std::to_string(MPIcf::my_rank()) + "_" +
                                             std::to_string(j + 1) + ".dat");
        //matlab::Export(convdiff.mat_[0], "mat_P1.dat");

        // Solve linear system
        convdiff.solve("mumps");

#if defined(algoim)

        // Compute area of domain in time quadrature point 0
        Fun_h funone(Wh, fun_one);
        double intGamma = integral_algoim<Levelset<2>, Fun_h>(funone, interface, 0, phi);
        double intOmega = integral_algoim<Levelset<2>, Fun_h>(funone, Thi, phi, 0);

        double errBulk = L2_norm_cut(b0h, fun_uBulk, phi, 0, 1);

        std::cout << "|| u-uex||_2 = " << errBulk << '\n';

#elif defined(triangle) || defined(quad)
        Fun_h funone(Wh, fun_one);
        double intGamma = integral(funone, interface, 0);
        double intOmega = integral(Thi, funone, 0);
        
        double errBulk = L2normCut(b0h, fun_uBulk, 0, 1);
        std::cout << "|| u-uex||_2 = " << errBulk << '\n';
#endif
        errors.at(j) = errBulk;

        gamma.at(j)     = std::fabs(intGamma - 2*pi);
        omega.at(j)     = std::fabs(intOmega-(pi));

        if ((iterations == 1)) {
            Fun_h sol(Wh, data_u0);
            Paraview<Mesh> background_mesh(Th, path_output_figures + "Th.vtk");
            ExpressionFunFEM<Mesh> vx(vel, 0, op_id);
            ExpressionFunFEM<Mesh> vy(vel, 1, op_id);
            background_mesh.add(vx, "velx");
            background_mesh.add(vy, "vely");

            Paraview<Mesh> writer(Thi, path_output_figures + "bulk.vtk");
            writer.add(b0h, "bulk", 0, 1);
            Fun_h uBex(Wh, fun_uBulk);
            Fun_h fB(Wh, fun_rhsBulk);
            Expression uuh(sol, 0, op_id);
            Expression uuex(uBex, 0, op_id);
            writer.add(uBex, "bulk_exact", 0, 1);
            writer.add(fB, "bulk_rhs", 0, 1);
            writer.add(g_Neumann, "neumann", 0, 1);

            writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
            writer.add(levelSet, "levelSet", 0, 1);
            writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh.vtk");
            writer.writeFaceStab(Thi, 0, path_output_figures + "Faces.vtk");
#if defined(algoim)
            writer.writeAlgoimQuadrature(Thi, phi, -1, path_output_figures + "AlgoimQuadrature.vtk");
#endif
        }

// Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
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
