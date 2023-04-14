/*
   We consider a time-dependent bulk problem on Omega2.
   We consider problems with both Neumann and Dirichlet boundary conditions,
   and we consider scheme II, III and the Reynold scheme.

   // Problem:
   Find u in Omega2 such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    on Omega2.

    // Numerical method:s
    A space-time Cutfem, using the level-set method,
    which allows for both dg and cg.

    Classical : Integration by parts on full convection term, no term added to
   make anti-symmetric Conservative: Reynold's transport theorem is used to make
   the bilinear form fulfill a conservation law
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
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

using namespace globalVariable;

namespace Diffusion {

        // RHS for surface variable
        double fun_rhs0(const R2 P, int elementComp) {
                R x = P.x, y = P.y;
                return (9*y*(3*x*x - y*y))/(x*x + y*y);

        }

        // Exact solution on surface
        double fun_uSurface(const R2 P,  int elementComp) {
                return 3.0*P.x*P.x*P.y - pow(P.y,3);
        }

        // Exact solution on surface for time-step t
        double fun_uSurfaceT(const R2 P,  int elementComp, double t) {
                return 3.0*P.x*P.x*P.y - pow(P.y,3);
        }

        // Level-set function
        double fun_levelSet(const R2 P, const int i) {
                return (P.x)*(P.x) + (P.y)*(P.y) - 1.0 - Epsilon;
        }

        double fun_one(const R2 P, const int i) { return 1.; }

}

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

using namespace Diffusion;

#define dg

#define use_h


int main(int argc, char **argv) {

    verbose = 0;    // 2 for more info

    MPIcf cfMPI(argc, argv);

    std::cout << std::setprecision(16);

    // Mesh settings and data objects
    const size_t iterations = 1; // number of mesh refinements   (set to 1 to
                                 // run only once and plot to paraview)
    int nx = 20, ny = 15;        // starting mesh size (only apply if use_n is defined)
    double h  = 0.1;             // starting mesh size
    
	std::array<double, iterations> errors;             // array to hold L2 errors vs h
    std::array<double, iterations> gamma_length_h;
    std::array<double, iterations> nxs; // array to hold mesh sizes
    std::array<double, iterations> nys; // array to hold mesh sizes
    std::array<double, iterations> hs;  // array to hold mesh sizes
    
    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Mesh
        double lx = 3., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        
        Mesh Th(nx, ny, -1.5, -1.5, lx, ly);

        std::cout << "Number of elements in background mesh: " << Th.nbElements() << "\n";

        // Paths to store data
        const std::string path_output_data = "../output_files/surface/laplace_beltrami/data/";
        const std::string path_figures     = "../output_files/surface/laplace_beltrami/paraview/";

        // Create directory if not already existent
        if (MPIcf::IamMaster()) {
            std::filesystem::create_directories(path_output_data);
            std::filesystem::create_directories(path_figures);
        }

        // Data file to hold problem data
        std::ofstream output_data(path_output_data + "data.dat", std::ofstream::out);
        
        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;
        
        std::cout << "h = " << h << "\n";
        std::cout << "nx = " << nx << "\n";
        std::cout << "ny = " << ny << "\n";
        
	#if defined(cg)
        // CG stabilization parameter
        double tau1 = 1e-6;
	#elif defined(dg)
		// DG penalty and stabilization parameters
		double lambda = 50;
		double tau1 = 1, tau2 = 1, tau3 = 0.1;
    #endif

        // Background FE Space, Time FE Space & Space-Time Space
	#if defined(cg)
		GFESpace<Mesh> Vh(Th, DataFE<Mesh>::P1);  // continuous basis functions
	#elif defined(dg)
		GFESpace<Mesh> Vh(Th, DataFE<Mesh>::P1dc);  // continuous basis functions
	#endif
        // Set up level-set function
        GFESpace<Mesh> Lh(Th, DataFE<Mesh>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Th, levelSet);

        // Create active meshes
        ActiveMesh<Mesh> ThGamma(Th);
        ThGamma.createSurfaceMesh(interface);
        std::cout << "Number of elements in active mesh: " << ThGamma.get_nb_element() << "\n";

        CutSpace Wh(ThGamma, Vh);
        ProblemOption problem_option;
        problem_option.order_space_element_quadrature_ = 5;     // 5 = default
        CutFEM<Mesh> surfactant(Wh, problem_option);

        double errL2 = 0., intGamma = 0.; 

        Fun_h funrhs(Wh, fun_rhs0);
		Fun_h funone(Wh, fun_one);

        Normal n;
        Tangent t;

        gnuplot::save(Th);
        //gnuplot::save(interface, "interface.dat");
        gnuplot::save(interface);

        // Test and Trial functions
        FunTest u(Wh, 1), v(Wh, 1); 

        // Scheme for diffusion
        surfactant.addBilinear(+innerProduct(gradS(u), gradS(v)), interface);

	#if defined(dg)
		surfactant.addBilinear(
            - innerProduct(average(gradS(u)*t,0.5,-0.5), jump(v)) 
            - innerProduct(jump(u), average(gradS(v)*t,0.5,-0.5))    
            + innerProduct(lambda/h*jump(u), jump(v))                             
            , interface
            , INTEGRAL_INNER_NODE_2D
        );

		// Stabilization
        surfactant.addFaceStabilization(
			+ innerProduct(tau1 * 1./h/h * jump(u), jump(v))
			+ innerProduct(tau2 * jump(grad(u) * n), jump(grad(v) * n)), ThGamma);
		
		surfactant.addBilinear(
                innerProduct(tau3*h*h*grad(u)*n, grad(v)*n)
                , interface
        );

	#elif defined(cg)

		// Stabilization
        surfactant.addFaceStabilization(
			+ innerProduct(tau1 * jump(grad(u) * n), jump(grad(v) * n)), ThGamma);

	#endif

        // Add RHS on surface
        surfactant.addLinear(+innerProduct(funrhs.expr(), v), interface);

		// Add Lagrange multiplier such that the averages of the test functions are zero
        surfactant.addLagrangeMultiplier(innerProduct(1.,v), 0., interface);
        
        matlab::Export(surfactant.mat_[0], "mat_triangle_" + std::to_string(j) + ".dat");
        
        // Solve linear system
        surfactant.solve("mumps");

        KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
        Fun_h us(Wh, dw);

        errL2 = L2normSurf(us, fun_uSurface, interface, 0, 1);

        intGamma = integral(funone, interface, 0);
        errors.at(j)         = errL2;
        gamma_length_h.at(j) = fabs(intGamma - 2*pi);

        if (iterations == 1) {

            Paraview<Mesh> writer(ThGamma, path_figures + "surfactant.vtk");

            Fun_h uS_ex(Wh, fun_uSurface);
            writer.add(us, "surfactant", 0, 1);
            writer.add(uS_ex, "surfactant_exact", 0, 1);
            writer.add(fabs(us.expr() - uS_ex.expr()), "surfactant_error");
            writer.add(levelSet, "levelSet", 0, 1);
            // writer.add(ls[2], "levelSet2", 0, 1);
        }


        // Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_h)
        h *= 0.5;
#endif
    }

    std::cout << "\n";
    std::cout << "Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    std::cout << "\n";

	std::cout << "\n";
    std::cout << "Length Gamma = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << gamma_length_h.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    std::cout << "\n";

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";
	
	std::cout << "\n";

    std::cout << "nx = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nxs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    return 0;
}
