/*  
   We consider a time-dependent bulk problem on Omega2. 
   We consider problems with both Neumann and Dirichlet boundary conditions,
   and we consider scheme II, III and the Reynold scheme.

   // Problem:
   Find u in Omega2 such that

    dt(u) + beta*grad(u) - D*laplace(u) = f,    on Omega2.

    // Numerical method:
    A space-time Cutfem, using the level-set method,
    which allows for both dg and cg.

    Scheme II  : Integration by parts on half of the convection term, term added to make anti-symmetric
    Scheme III : Integration by parts on full convection term, no term added to make anti-symmetric
*/

// TODO: Allow product between double and CutFEM_R2 
// TODO: Allow fabs on Normal and CutFEM_R2 

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
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"

// Numerical examples

namespace Frachon1 {

    double c0 = 0.5;

    R fun_levelSet(const R2 P, const int i, double t) {
        return -P.x - P.y + c0;
    }

    R fun_levelSet(const R2 P, const int i) {
        return -P.x - P.y + c0;
    }

    R fun_boundary(const R2 P, int elementComp, double t) {
        return sin(pi*(P.x + P.y- 4*t));
    }

    R fun_uBulk(const R2 P, int elementComp, int domain, double t) {
        // return 2*sin(pi*(P.x+P.y-3*t));
        if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
        else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
    }

    R fun_uBulkInit(const R2 P, int elementComp, int domain) {
        if(domain == 0) return sin(pi*(P.x + P.y));
        else return 4./3*sin(4./3*pi*(P.x + P.y - c0/4));
    }

    // Velocity field
    R fun_velocity(const R2 P, const int i, int domain) {
        if (domain == 0) {
            if (i == 0) return 3.;
            else return 1.;
        }
        else {
            if (i == 0) return 2.;
            else return 1.;
        }
    }

}

namespace Frachon2 {

    double c0 = 0.25;

    R one(const R2 P, int elementComp, int domain, double t) {return 1.;}

    R fun_levelSet(const R2 P, const int i, double t) {
        return -P.x - P.y + c0;
    }

    R fun_levelSet(const R2 P, const int i) {
        return -P.x - P.y + c0;
    }

    R fun_boundary(const R2 P, int elementComp, double t) {
        return 0.;
    }

    R fun_uBulk(const R2 P, int elementComp, int domain, double t) {
        // return 2*sin(pi*(P.x+P.y-3*t));
        if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
        else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
    }

    R fun_uBulkInit(const R2 P, int elementComp, int domain) {

        if(domain == 1) return 0;
        double xs = -0.3, ys = -0.3;
        double r = 0.3;
        double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
        if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
        else return 0.;
    }

    // Velocity field
    R fun_velocity(const R2 P, const int i, int domain) {
        if (domain == 0) {
            if (i == 0) return 3.;
            else return 1.;
        }
        else {
            if (i == 0) return 1.;
            else return 2.;
        }
    }

}

// Setup two-dimensional class types
const int d = 2;
typedef Mesh2 Mesh;
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<Mesh2> Fun_h;

// Note: standard for this program is to solve the equations on the inner domain Omega 2.
// To solve on Omega 1, define the word "omega1" for the pre-processor (only works for example 1)


// Choose Discontinuous or Continuous Galerkin method (options: "dg", "cg")
#define dg
// Set numerical example (options: "frachon1", "flower", "example1")
#define frachon2
// Set scheme for the dg method (options: "conservative", "classical")
#define conservative
// Set stabilization method (options: "fullstab", "macro") 
#define fullstab       

#define use_vel   // (options: "use_vel", "use_beta")

#ifdef frachon1
    using namespace Frachon1;
#elif defined(frachon2)
    using namespace Frachon2;
#endif

int main(int argc, char** argv) {
    
    // Mesh settings and data objects
    const size_t iterations = 1;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 11, ny = 11;       // starting mesh size
    const double lx = 2., ly = 2.;    // domain length

#ifdef frachon1    
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeCoupled/Frachon1/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeCoupled/Frachon1/paraview/";
#elif defined(frachon2)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeCoupled/Frachon2/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeCoupled/Frachon2/paraview/";
#endif

    // Initialize MPI
    MPIcf cfMPI(argc,argv);

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::experimental::filesystem::create_directories(pathOutputFolder);
        std::experimental::filesystem::create_directories(pathOutputFigures);
    }

    // Data file to hold problem data
    std::ofstream outputData(pathOutputFolder + "data.dat", std::ofstream::out);

    // Arrays to hold data
    std::array<double, iterations> errorsBulk;          // array to hold bulk errors
    std::array<double, iterations> hs;                  // array to hold mesh sizes

    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {

        // Define background mesh
        Mesh Th(nx, ny, -1., -1., lx, ly);

    //// Parameters
        
        // Mesh size
        double h = lx/(nx);
        hs.at(j) = h;
        int divisionMeshSize = 2*3*sqrt(10);
        //int divisionMeshSize = 2;

        // Time
        double dT = h / 3 / sqrt(10) * 0.5; //h/divisionMeshSize; // Time step size
        double tfinal = .5;            // Final time
        GTime::total_number_iteration = (int)(tfinal/dT);
        dT = tfinal / GTime::total_number_iteration;
        GTime::time_step = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "Iteration " << j + 1 << "/" << iterations << std::endl;
            std::cout << "h = " << h << std::endl;
            std::cout << "nx = " << nx << std::endl;
            std::cout << "dT = " << dT << std::endl;
        }

        // Penalty parameter for outer boundary
        double lambdaB = 3;    // coefficient on the outer boundary

        CutFEM_R2 beta({R2(3,1), R2(2,1)});

        CutFEMParameter lambda(0., 1.);

        CutFEMParameter lambdaE(3, 2);  //  ||beta||_inf in Omega_1/Omega_2

        // DG stabilization parameters
        double tau20 = 1e-2, tau21 = 1e-2;      // bulk


        // Background FE Space, Time FE Space & Space-Time Space

        FESpace2 Vh(Th, DataFE<Mesh>::P1dc);        // discontinuous basis functions
        // 1D Time mesh
        Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        // Quadrature data
        const QuadratureFormular1d& qTime(*Lobatto(3));
        const Uint nbTime = qTime.n;
        const Uint ndfTime = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        double dt_levelSet = dT/(nbTime-1);
        vector<Fun_h> ls(nbTime);

        for (int i=0; i<nbTime; i++) ls[i].init(Lh, fun_levelSet, 0.);

        // Declare time dependent interface (although static interface)
        TimeInterface<Mesh> interface(qTime);
        
        // Bulk-Surface Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;
        
        // computation of the interface in each of the three quadrature points
        for (int i = 0; i < nbTime; ++i) {
        
            interface.init(i, Th, ls[i]);

        }

        int iter = 0;
        double q0_0, q0_1, qp_0, qp_1;
        double intF = 0, intG = 0;      // hold integrals of rhs and Neumann bcs
        // Iterate over time-slabs
        while (iter < GTime::total_number_iteration) {
            
            
            GTime::current_iteration = iter;
            const TimeSlab &In(Ih[iter]);

            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << GTime::total_number_iteration << std::endl;
            std::cout << " Time      \t : \t" << GTime::current_time() << std::endl;
    
            ActiveMesh<Mesh> Khi(Th, interface);
            
            // Velocity field

            Lagrange2 FEvelocity(0);
            CutSpace Wh(Khi, Vh);

            Fun_h one_fun(Wh, In, one);

            Expression2 fun(one_fun, 0, 0, 0);

            std::cout << integral(Khi, In, fun, INTEGRAL_BOUNDARY, qTime, {1,2,4}) - 6*dT << std::endl;
            return 0;
        }
    }
}
