#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
// #include "cfmpi.hpp"
#include "finiteElement.hpp"
// #include "baseCutProblem.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"

using namespace globalVariable;

double fun_boundary(double *P, int elementComp) { return 2 * sin(2 * pi * P[0]) * sin(4 * pi * P[1]); }
double fun_rhs(double *P, int elementComp) { return 40 * pow(pi, 2) * sin(2 * pi * P[0]) * sin(4 * pi * P[1]); }
double fun_exact(double *P, int elementComp) { return 2 * sin(2 * pi * P[0]) * sin(4 * pi * P[1]); }
double fun_zero(double *P, int elementComp) { return 0.; }
double fun_test(double *P, int elementComp) { return P[0] * P[1]; }

#define quad

int main(int argc, char **argv) {

    const size_t iterations = 6;

    // 2 DIMENSIONAL PROBLEM
    // =====================================================
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

    // INITIALIZE MPI
    // =====================================================
    // MPIcf cfMPI(argc, argv);

    // MESH AND PROBLEM PARAMETERS
    // =====================================================
    double h = 0.1;
    int lx = 1., ly = 1.;

    std::array<double, iterations> errors; // array to hold bulk errors
    std::array<double, iterations> hs;     // array to hold mesh sizes

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {
        int nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;

        // CONSTRUCTION OF THE MESH AND FINITE ELEMENT SPACE
        // =====================================================
        Mesh Th(nx, ny, 0., 0., lx, ly);
        FESpace Vh(Th, DataFE<Mesh>::P3);

        hs.at(j)      = h;
        double lambda = 100.;

        Th.info();
        std::cout << "h = " << h << "\n";
        std::cout << "nx = " << nx << "\n";
        std::cout << "ny = " << ny << "\n";

        // OBJECTS NEEDED FOR THE PROBLEM
        // =====================================================

        ProblemOption option;
        option.order_space_element_quadrature_ = 6;
        FEM<Mesh> poisson(Vh, option);

        Fun_h fh(Vh, fun_rhs);
        Fun_h gh(Vh, fun_exact);
        Fun_h funTest(Vh, fun_test); // computer Neumann BC
        Fun_h gh_zero(Vh, fun_zero);
        FunTest u(Vh, 1), v(Vh, 1);
        Normal n;

        // ASSEMBLY OF THE LINEAR SYSTEM
        // =====================================================
        poisson.addBilinear(innerProduct(grad(u), grad(v)), Th);
        // poisson.addBilinear(innerProduct(u, v), Th);
        poisson.addLinear(innerProduct(fh.expr(), v), Th);

        poisson.setDirichlet(gh_zero, Th);

        // poisson.addBilinear(-innerProduct(grad(u) * n, v) - innerProduct(u, grad(v) * n) +
        //                         innerProduct(u, lambda / h * v) // added penalty
        //                     ,
        //                     Th, INTEGRAL_BOUNDARY);

        // poisson.addLinear(- innerProduct(gh.expr(), grad(v) * n) + innerProduct(gh.expr(), lambda / h * v), Th,
        //                   INTEGRAL_BOUNDARY);

        // RESOLUTION OF THE LINEAR SYSTEM
        // =====================================================
        poisson.solve("umfpack");

        // COMPUTE THE L2 ERROR
        // =====================================================
        Fun_h femSolh(Vh, poisson.rhs_);

        R errU = L2norm(femSolh, fun_exact, 0, 1);
        // R errU = L2norm(femSolh, fun_rhs, 0, 1);
        std::cout << "errU = " << errU << "\n";

        errors.at(j) = errU;

        // PRINT THE SOLUTION TO PARAVIEW
        // =====================================================
        // Paraview<Mesh> writer(Th, "poisson" + std::to_string(j) + ".vtk");
        // writer.add(femSolh, "poisson", 0, 1);
        // Fun_h uBex(Vh, fun_exact);
        // Fun_h fB(Vh, fun_rhs);
        // auto test_dx(dx(funTest.expr()));
        // auto test_dy(dy(funTest.expr()));
        // writer.add(uBex, "bulk_exact", 0, 1);
        // // writer.add(fabs(femSolh.expr() - uBex.expr()), "bulk_error");
        // writer.add(fabs(femSolh.expr() - fB.expr()), "bulk_error");
        // writer.add(fB, "bulk_rhs", 0, 1);
        // writer.add(funTest, "test", 0, 1);
        // writer.add(test_dx, "test_dx");
        // writer.add(test_dy, "test_dy");

        h *= 0.5;
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

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
}
