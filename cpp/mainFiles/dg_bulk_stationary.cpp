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
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../num/matlab.hpp"
#include "../num/redirectOutput.hpp"
#include "paraview.hpp"

using namespace globalVariable; // to access some globally defined constants

namespace Circle {

// RHS for bulk variables
double fun_rhsBulk(double *P, int elementComp, const R t) {
    R x = P[0], y = P[1];
    return -(exp(-x * x - y * y - t / 4 + 1) *
             (120 * x * x * x * x * y + 96 * x * x * x * y + 80 * x * x * y * y * y + 72 * x * x * y * y -
              465 * x * x * y - 36 * x * x - 32 * x * y * y * y - 96 * x * y - 40 * y * y * y * y * y -
              24 * y * y * y * y + 155 * y * y * y + 36 * y * y)) /
           10;
}

double fun_uBulk(double *P, int elementComp, const R t) {
    return 2.0 * exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - pow(P[1], 3)) * exp(-0.25 * t);
}

double fun_uBulkD(double *P, int elementComp, int domain, const R t) {
    return 2.0 * exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - pow(P[1], 3)) * exp(-0.25 * t);
}

double fun_uBulkInit(double *P, int elementComp) {
    return 2.0 * exp(1.0 - P[0] * P[0] - P[1] * P[1]) * (3.0 * P[0] * P[0] * P[1] - pow(P[1], 3));
}
double fun_neumann_Gamma(double *P, const int i, const R t) { return 0.; }

// Level-set function
double fun_levelSet(double *P, const int i) {
    R r0 = 0.5, r1 = 0.15;
    return -(sqrt((P[0]) * (P[0]) + (P[1]) * (P[1])) - 1);
}

double fun_velocity(double *P, const int i) {
    if (i == 0)
        return 0.8;
    else
        return 0.6;
}
} // namespace Circle

namespace Flower {

    // RHS for variable in Omega 1
    double fun_rhs1(double *P) {
        return -exp(- pow(P[0],2) - pow(P[1],2) + 1)*(12*pow(P[0],4)*P[1] + 3*pow(P[0],3) + 8*pow(P[0],2)*pow(P[1],3) - 48*pow(P[0],2)*P[1] - 9*P[0]*pow(P[1],2) - 4*pow(P[1],5) + 16*pow(P[1],3));
    }

    // RHS for variable in Omega 2
    double fun_rhs2(double *P) {
        return -2*exp(- pow(P[0],2) - pow(P[1],2) + 1)*(6*pow(P[0],4)*P[1] + 3*pow(P[0],3) + 4*pow(P[0],2)*pow(P[1],3) - 24*pow(P[0],2)*P[1] - 9*P[0]*pow(P[1],2) - 2*pow(P[1],5) + 8*pow(P[1],3));
    }

    // RHS for bulk variables
    double fun_rhsBulk(double *P, int elementComp) {
        return fun_rhs1(P);
        //return fun_rhs2(P);
    }

    // Exact solution in the bulk
    double fun_uBulk(double *P,  int elementComp, int domain) {
        return exp(1.0-P[0]*P[0]-P[1]*P[1])*( 3.0*P[0]*P[0]*P[1] - pow(P[1],3));
        //return 2.0*exp(1.0-P[0]*P[0]-P[1]*P[1])*( 3.0*P[0]*P[0]*P[1] - pow(P[1],3));
    }

    double g_Neumann(double *P, int elementComp) {
        R x = P[0], y = P[1];
        return -(y*exp(- x*x - y*y + 1)*(6*x*x*x*x + 4*x*x*y*y - 9*x*x - 2*y*y*y*y + 3*y*y))/sqrt((x*x + y*y));
    }

    // Level-set function
    double fun_levelSet(double *P, const int i) {
        R r0 = 0.5, r1 = 0.15;
        return (sqrt((P[0])*(P[0]) + (P[1])*(P[1])) - r0 - r1*cos(5*atan2(P[1], P[0])));
    }

    // Velocity Field
    R fun_velocity(double *P, const int i){
        if(i == 0) return P[1];
        else       return -P[0];
    }
}


// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// The below parameters can be varied according to the options to use different methods,
// different numerical examples, different subdomains and boundary conditions etc.

//* Set numerical example (options: "example1", "lehrenfeld")
#define example1
//* Choose domain to solve on (options: "omega1", "omega2")
#define omega1
// If "omega1":
// Set type of BCs on outer boundary (options: "dirichlet1" or "neumann1")
#define dirichlet1
// Set type of BCs on interface (options: "dirichlet2" or "neumann2")
#define dirichlet2

// If "omega2":
// Set type of BCs on interface (options: "dirichlet", "neumann")
#define neumann

//* Set stabilization method (options: "fullstab", "macro")
#define fullstab

#define use_h // to set mesh size using the h parameter. Write use_n to decide
              // using nx, ny.

using namespace Flower;

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);
    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h = 0.1;              // starting mesh size

    // Paths to store data
    const std::string path_output_data    = "../output_files/bulk/stationary/data/";
    const std::string path_output_figures = "../output_files/bulk/stationary/paraview/";

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

        hs.at(j) = h;

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

        const double D = 1.;

        double lambda = 1.; // Nitsche's method penalty parameter

        // DG interior penalty parameters
        const double lambda_A = 500, lambda_B = 0.5;

        // DG stabilization parameter
        const double tau0 = 1., tau1 = 1.;

        FESpace2 Vh(Th, DataFE<Mesh>::P1dc); // discontinuous basis functions

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        FESpace2 Vh2(Th,
                     DataFE<Mesh>::P2);     // higher order space for interpolation
        FESpace2 Vh3(Th, DataFE<Mesh>::P3); // continuous basis functions

// Velocity field
        Lagrange2 FEvelocity(1);
        FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Th, levelSet);

        // Create active meshes
        ActiveMesh<Mesh> Thi(Th);

        Thi.truncate(interface, -1); // remove part with negative sign of level

        // Cut FE space
        CutSpace Wh(Thi, Vh);

        CutFEM<Mesh> convdiff(Wh);

        // Objects needed for the weak form
        Normal n;
        Tangent t;

        // Right hand side functions
        Fun_h f(Vh, fun_rhsBulk);
        Fun_h g(Vh, fun_uBulk);                 // create an FE-function of the exact bulk
                                                    // solution Omega2
        Fun_h g_Neumann(Vh, g_Neumann); // computer Neumann BC

        // Test and Trial functions
        FunTest u(Wh, 1), v(Wh, 1);

        // Data for initial solution
        Rn data_u0; // initial data total
        KN_<R> data_B0(convdiff.rhs_(SubArray(Wh.NbDoF(),0))); // data bulk

        // Make function objects to use in innerProducts
        Fun_h b0h(Wh, data_B0);

        //** Assembling linear and bilinear forms

        //* Diffusion term
        convdiff.addBilinear(+innerProduct(D * grad(u), grad(v)), Thi);

        // Diffusion on inner edges (E_h)
        convdiff.addBilinear(-innerProduct(D * average(grad(u) * n), jump(v)) -
                                    innerProduct(D * jump(u), average(grad(v) * n)) +
                                    innerProduct(lambda_A / h * jump(u), jump(v)),
                                Thi, INTEGRAL_INNER_EDGE_2D);

        //* Convection term
        convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Thi);

        // Convection on inner edges (E_h)
        convdiff.addBilinear(+innerProduct(average(vel * n * u), jump(v)) +
                                    innerProduct(0.5 * fabs(vel * n) * jump(u), jump(v)),
                                Thi, INTEGRAL_INNER_EDGE_2D);

        //* Stabilization

#if defined(fullstab)

        convdiff.addFaceStabilization(+innerProduct(1. / h * tau0 * jump(u), jump(v)) +
                                            innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)),
                                        Thi);

#elif defined(macro)

        MacroElementPartition<Mesh> TimeMacro(Thi, 0.2);

        // Visualize macro elements
        if (iterations == 1 && h > 0.01) {
            Paraview<Mesh> writerMacro(Th, path_output_figures + "Th.vtk");
            writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
            writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
            writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

            // domain = 0,

            writerMacro.writeFaceStab(
                Thi, 0, path_output_figures + "FullStabilization.vtk");
            writerMacro.writeActiveMesh(Thi,
                                        path_output_figures + "ActiveMesh.vtk");
            writerMacro.writeMacroElement(TimeMacro, 0,
                                            path_output_figures + "macro.vtk");
            writerMacro.writeMacroInnerEdge(
                TimeMacro, 0, path_output_figures + "macro_inner_edge.vtk");
            writerMacro.writeMacroOutterEdge(
                TimeMacro, 0, path_output_figures + "macro_outer_edge.vtk");
            writerMacro.writeSmallElements(
                TimeMacro, 0, path_output_figures + "small_element.vtk");
        }

        // Stabilization of the bulk
        // convdiff.mat_.clear();
        convdiff.addFaceStabilization(+innerProduct(1. / h * tau0 * jump(u), jump(v)) +
                                            innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)),
                                        Thi, TimeMacro);

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
                                Thi, INTEGRAL_BOUNDARY);

        convdiff.addLinear(+innerProduct(g.expr(), 0.5 * fabs(vel * n) * v) -
                                innerProduct(g.expr(), 0.5 * (vel * n) * v) -
                                innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                            Thi, INTEGRAL_BOUNDARY);

        convdiff.addBilinear(+innerProduct(u, 0.5 * fabs(vel * n) * v) + innerProduct(u, 0.5 * (vel * n) * v) -
                                    innerProduct(D * grad(u) * n, v)   // from IBP
                                    - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                    + innerProduct(u, lambda / h * v)  // added penalty
                                ,
                                interface);

        convdiff.addLinear(+innerProduct(g.expr(), 0.5 * fabs(vel * n) * v) -
                                innerProduct(g.expr(), 0.5 * (vel * n) * v) -
                                innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
                            interface);

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
        convdiff.addLinear(+innerProduct(f.expr(), v), Thi);

        
        matlab::Export(convdiff.mat_[0], path_output_data + "mat_h_dg" + std::to_string(h) + "_" +
                                                    std::to_string(j + 1) + ".dat");

        // Solve linear system
        convdiff.solve("mumps");

        double errBulk = L2normCut(b0h, fun_uBulk, 0, 1);
        std::cout << "|| u-uex||_2 = " << errBulk << '\n';

        errors.at(j) = errBulk;

        // #ifdef USE_MPI
        //          if ((iterations == 1) && MPIcf::IamMaster()) {
        // #else
        if (MPIcf::IamMaster() && (iterations == 1)) {
            Paraview<Mesh> writer(Thi, path_output_figures + "bulk_dg.vtk");
            writer.add(b0h, "bulk", 0, 1);
            Fun_h uBex(Wh, fun_uBulk);
            Fun_h fB(Wh, fun_rhsBulk);
            Expression2 uuh(b0h, 0, op_id);
            Expression2 uuex(uBex, 0, op_id);
            writer.add(uBex, "bulk_exact", 0, 1);
            writer.add(fB, "bulk_rhs", 0, 1);
            writer.add(fabs(b0h.expr() - uBex.expr()), "bulk_error");
            // writer.add(fabs(uuh - uuex), "bulk_error");
            writer.add(levelSet, "levelSet0", 0, 1);
        }


// Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
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

    return 0;
}
