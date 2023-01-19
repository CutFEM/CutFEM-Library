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

using namespace globalVariable; // to access some globally defined constants

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

// Initial solution surface
R fun_uSurfInit(double *P, const int i) { return 0.; }

// Exact solution surface
R fun_uSurf(double *P, const int i, const R t) {
    double r0 = 1., x = P[0], y = P[1];
    return 0;
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

// RHS fB bulk
R fun_rhsSurf(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0;
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

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// Note: standard for this program is to solve the equations on the inner domain
// Omega 2. To solve on Omega 1, define the word "omega1" for the pre-processor
// (only works for example 1)

//* Set numerical example (options: "example1", "lehrenfeld")
#define example1
//* Set parameter D (options: "convection_dominated" means D=0.01, else D=1)
#define convection_dominated

//* Choose domain to solve on (options: "omega1", "omega2")
#define omega1
// If "omega1":
// Set type of BCs on outer boundary (options: "dirichlet1" or "neumann1")
#define neumann1
// Set type of BCs on interface (options: "dirichlet2" or "neumann2")
#define neumann2

// If "omega2":
// Set type of BCs on interface (options: "dirichlet", "neumann")
#define dirichlet

//* Set scheme for the method (options: "classical", "conservative".)
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

// [0.0128738, 0.00319005, 0.00074366]

// Do not touch
#ifdef example1
// #ifdef omega1
// using namespace Example1_Omega1;
// #else
using namespace Example1; // on Omega 2
// #endif
#endif
#if defined(lehrenfeld)
using namespace Lehrenfeld_Convection_Dominated;
#endif

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);
    // Mesh settings and data objects
    const size_t iterations = 2; // number of mesh refinements   (set to 1 to run
                                 // only once and plot to paraview)
    int nx = 15, ny = 15;        // starting mesh size
    double h  = 0.0125;          // starting mesh size
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

    // Data file to hold problem data

    std::ofstream outputData(pathOutputFolder + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errors; // array to hold bulk errors
    std::array<double, iterations> errors_surface; // array to hold surface errors
    std::array<double, iterations> hs;     // array to hold mesh sizes
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
        const int divisionMeshSize = 2;

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

        const double D      = 0.01;
        const double DGamma = 1.;

        double lambda = 1.; // Nitsche's method penalty parameter

        // CG stabilization parameter
        const double tau1 = 0.1;

        FESpace2 Vh(Th, DataFE<Mesh>::P1); // continuous basis functions

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
        double q0_0, q0_1, qp_0, qp_1; // integral values to be computed
        double intF = 0, intG = 0;     // hold integrals of rhs and Neumann bcs

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

            ls.begin()->swap(ls[nbTime - 1]);

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

#ifdef omega1
            Thi.truncate(interface, -1); // remove part with negative sign of level
#elif defined(omega2)
            Thi.truncate(interface, 1); // remove part with positive sign of level
                                        // set to get inner domain
#endif
            // Cut FE space
            CutSpace Wh(Thi, Vh);

            // Surface active mesh
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            // Surface FE space
            CutSpace WhGamma(ThGamma, Vh);

            // Add time slab to cut space
            convdiff.initSpace(Wh, In);

            // Add time slab to surface space on top of the data
            convdiff.add(WhGamma, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Vh, In, fun_rhsBulk);
            Fun_h g(Vh, In, fun_uBulk); // create an FE-function of the exact bulk solution

            Fun_h fS(Vh, In, fun_rhsSurf);
            Fun_h gS(Vh, In, fun_uSurf);

            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);         // computed Neumann BC
            Fun_h g_Neumann_left(Wh, In, fun_neumann_left);     // label 1
            Fun_h g_Neumann_bottom(Wh, In, fun_neumann_bottom); // label 4
            Fun_h g_Neumann_right(Wh, In, fun_neumann_right);   // label 2
            Fun_h g_Neumann_top(Wh, In, fun_neumann_top);       // label 3

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);             // bulk
            FunTest uS(WhGamma, 1), vS(WhGamma, 1); // bulk

            // Data for initial solution

            // variable are [vB, vS]
            int idx_s0  = Wh.NbDoF() * In.NbDoF(); // where v0 start in the solution array
            int tot_dof = (int)convdiff.get_nb_dof();

            // Data for initial solution
            Rn data_u0(tot_dof, 0.); // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));           // initial data bulk
            KN_<R> data_S0(data_u0(SubArray(WhGamma.NbDoF(), idx_s0))); // initial data surface

            if (iter == 0) {
                interpolate(Wh, data_B0, fun_uBulkInit);
                interpolate(WhGamma, data_S0, fun_uSurfInit);
            }

            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);
            Fun_h s0h(WhGamma, data_S0);

            // Plot initial solution in paraview
            // #ifdef USE_MPI
            // if (iter == 0 && MPIcf::IamMaster()) {
            // #else
            if (iter == 0) {
                // #endif
                Paraview<Mesh> writerInitial(Thi, pathOutputFigures + "BulkInitial.vtk");
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

            // Bulk
            convdiff.addBilinear(+innerProduct(dt(u), v), Thi, In);
            // Surface
            convdiff.addBilinear(+innerProduct(dt(uS), vS), interface, In);

            // Time penalty term bulk LHS
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);

            // Time penalty term bulk RHS
            convdiff.addLinear(+innerProduct(b0h.expr(), v), Thi, 0, In);

            // Time penalty term surface LHS
            convdiff.addBilinear(+innerProduct(uS, vS), *interface(0), In, 0);

            // Time penalty term surface RHS
            convdiff.addLinear(+innerProduct(s0h.expr(), vS), *interface(0), In, 0);

#endif

            //* Diffusion terms

            // Bulk
            convdiff.addBilinear(+innerProduct(D * grad(u), grad(v)), Thi, In);

            // Surface
            convdiff.addBilinear(+innerProduct(DGamma * gradS(uS), gradS(vS)), interface, In);

            //* Convection terms

#if defined(classical)
            // Bulk
            convdiff.addBilinear(+innerProduct((vel.exprList() * grad(u)), v), Thi, In);

            // Surface
            convdiff.addBilinear(+innerProduct((vel.exprList() * grad(uS)), vS) +
                                     innerProduct(uS * divS(vel), vS),
                                 interface, In);

#elif defined(conservative)
            convdiff.addBilinear(-innerProduct(u, (vel.exprList() * grad(v))), Thi, In);
#endif

            //* Stabilization

#if defined(fullstab)

            convdiff.addFaceStabilization(+innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)), Thi, In);
            convdiff.addFaceStabilization(+innerProduct(tau1 * jump(grad(uS) * n), jump(grad(vS) * n)), ThGamma, In);

#elif defined(macro)

            TimeMacroElement<Mesh> TimeMacro(Thi, qTime, 0.125);

            // Visualize macro elements
            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, pathOutputFigures + "Th" + std::to_string(iter + 1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);

                // domain = 0,

                writerMacro.writeFaceStab(Thi, 0,
                                          pathOutputFigures + "FullStabilization" + std::to_string(iter + 1) + ".vtk");
                writerMacro.writeActiveMesh(Thi, pathOutputFigures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
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
            convdiff.addFaceStabilization(+innerProduct(h * tau1 * jump(grad(u) * n), jump(grad(v) * n)), Thi, In,
                                          TimeMacro);

            // matlab::Export(convdiff.mat_, "mat.dat");
            // getchar();

#endif

            //* Boundary conditions

// If Omega1
#ifdef omega1

//* Dirichlet on both outer and inner boundary
#if defined(dirichlet1) && defined(dirichlet2)
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

            // Dirichlet inner
            convdiff.addBilinear(-innerProduct(D * grad(u) * n, v)      // from IBP
                                     - innerProduct(u, D * grad(v) * n) // added to make symmetric
                                     + innerProduct(u, lambda / h * v)  // added penalty
                                 ,
                                 interface, In);

            convdiff.addLinear(-innerProduct(g.expr(), D * grad(v) * n) + innerProduct(g.expr(), lambda / h * v),
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
            // TODO: The below functions do not work, I believe
            // convdiff.addLinear(innerProduct(g_Neumann_left.expr(), v), Thi, INTEGRAL_BOUNDARY, In,
            // (std::list<int>){1}); convdiff.addLinear(innerProduct(g_Neumann_bottom.expr(), v), Thi,
            // INTEGRAL_BOUNDARY, In, (std::list<int>){4}); convdiff.addLinear(innerProduct(g_Neumann_right.expr(), v),
            // Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){2}); convdiff.addLinear(innerProduct(g_Neumann_top.expr(),
            // v), Thi, INTEGRAL_BOUNDARY, In, (std::list<int>){3});

#if defined(conservative)
            convdiff.addBilinear(innerProduct(vel.exprList() * u * n, v), Thi, INTEGRAL_BOUNDARY, In);
#endif
            // Neumann inner
            //convdiff.addLinear(-innerProduct(g_Neumann.expr(), v), interface, In);
            //convdiff.addBilinear(innerProduct(jump(u,uS), jump(v,vS)), interface, In);
            convdiff.addBilinear(innerProduct(u-uS, v-vS), interface, In);

#endif

// If Omega2
#elif defined(omega2)
#ifdef neumann
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
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

            // Add RHS in bulk
            convdiff.addLinear(+innerProduct(f.expr(), v), Thi, In);
            convdiff.addLinear(+innerProduct(fS.expr(), vS), interface, In);

            if (iter == total_number_iteration - 1)
                matlab::Export(convdiff.mat_[0],
                               pathOutputFolder + "mat_h" + std::to_string(h) + "_" + std::to_string(j + 1) + ".dat");

            // Solve linear system
            convdiff.solve("umfpack");

            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            // Compute conservation error
            if (iterations == 1) {

                intF = integral(Thi, In, f, 0, qTime);
#ifdef neumann
                intG = integral(g_Neumann, In, interface, 0);
#endif

                Fun_h funuh(Wh, data_u0);

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
                           << ((q_1 - qp_1) - intF - intG) << '\n';
                qp_1 = q_1;
            }


            Rn sol(Wh.get_nb_dof(), 0.);
            sol += data_u0(SubArray(Wh.get_nb_dof(), 0));       // add first time dof
            sol += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));     // add second time dof
            Fun_h funuh(Wh, sol);
            double errBulk = L2normCut(funuh, fun_uBulkD, current_time + dT, 0, 1);
            std::cout << " t_n -> || u-uex||_2 = " << errBulk << '\n';

            Rn solsurf(WhGamma.get_nb_dof(), 0.);
            solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0));  
            solsurf += data_u0(SubArray(WhGamma.get_nb_dof(), idx_s0 + WhGamma.get_nb_dof()));
            Fun_h funuh_surf(WhGamma, solsurf);
            double errSurf = L2normSurf(funuh_surf, fun_uSurf, *interface(lastQuadTime), current_time+dT, 0, 1);
            std::cout << " t_n -> || uS-uSex||_2 = " << errSurf << '\n';

            errors.at(j) = errBulk;
            errors_surface.at(j) = errSurf;

            // #ifdef USE_MPI
            //          if ((iterations == 1) && MPIcf::IamMaster()) {
            // #else
            if ((iterations == 1)) {
                Fun_h sol(Wh, data_u0);
                Paraview<Mesh> writer(Thi, pathOutputFigures + "Bulk" + std::to_string(iter + 1) + "DG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, current_time);
                Fun_h fB(Wh, fun_rhsBulk, current_time);
                Expression2 uuh(sol, 0, op_id);
                Expression2 uuex(uBex, 0, op_id);
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB, "bulk_rhs", 0, 1);
                // writer.add(fabs(uuh - uuex), "bulk_error");
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

    std::cout << '\n';
    std::cout << "Errors Surface = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_surface.at(i);
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
