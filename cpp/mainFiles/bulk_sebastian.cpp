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
namespace Circle {

    // RHS for bulk variables
    double fun_rhsBulk(const R2 P, int elementComp, const R t) {
        R x = P.x, y = P.y;
        return -(exp(- x*x - y*y - t/4 + 1)*(120*x*x*x*x*y + 96*x*x*x*y + 80*x*x*y*y*y + 72*x*x*y*y - 465*x*x*y - 36*x*x - 32*x*y*y*y - 96*x*y - 40*y*y*y*y*y - 24*y*y*y*y + 155*y*y*y + 36*y*y))/10;
    }

    double fun_uBulk(const R2 P,  int elementComp, const R t) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3))*exp(-0.25*t);
    }

    double fun_uBulkD(const R2 P,  int elementComp, int domain, const R t) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3))*exp(-0.25*t);
    }

    double fun_uBulkInit(const R2 P, int elementComp) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3));
    }

    // Level-set function
    double fun_levelSet(const R2 P, const int i) {
        R r0 = 0.5, r1 = 0.15;
        return -(sqrt((P.x)*(P.x) + (P.y)*(P.y)) - 1);
    }

    double fun_velocity(const R2 P, const int i) {
        if (i == 0) return 0.8;
        else return 0.6;
    }
}

namespace Flower {

    // RHS for bulk variables
    double fun_rhsBulk(const R2 P, int elementComp, const R t) {
        R x = P.x, y = P.y;
        return -(exp(- x*x - y*y - t/4 + 1)*(120*x*x*x*x*y + 96*x*x*x*y + 80*x*x*y*y*y + 72*x*x*y*y - 465*x*x*y - 36*x*x - 32*x*y*y*y - 96*x*y - 40*y*y*y*y*y - 24*y*y*y*y + 155*y*y*y + 36*y*y))/10;
    }

    double fun_uBulk(const R2 P,  int elementComp, const R t) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3))*exp(-0.25*t);
    }

    double fun_uBulkD(const R2 P,  int elementComp, int domain, const R t) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3))*exp(-0.25*t);
    }

    double fun_uBulkInit(const R2 P, int elementComp) {
        return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*(3.0*P.x*P.x*P.y - pow(P.y,3));
    }

    // Level-set function
    double fun_levelSet(const R2 P, const int i) {
        R r0 = 0.5, r1 = 0.15;
        return -(sqrt((P.x)*(P.x) + (P.y)*(P.y)) - r0 - r1*cos(5*atan2(P.y, P.x)));
    }

    double fun_velocity(const R2 P, const int i) {
        if (i == 0) return 0.8;
        else return 0.6;
    }
}

namespace Example1 {
    /* This works for running Test – i.e. a pure bulk problem on Omega_2. */

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) - 0.17;
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return (pi*cos(pi*y)*sin(pi*x)*(2*cos(pi*t)*cos(pi*t) - 1)*((7*sin(pi*t))/25 - x + 0.5))
                /(250*sqrt(pow(y + (7*cos(t*pi))/25 - 0.5,2) + pow((7*sin(t*pi))/25 - x + 0.5,2)))
                - (pi*cos(pi*x)*sin(pi*y)*(2*cos(pi*t)*cos(pi*t) - 1)*(y + (7*cos(pi*t))/25 - 0.5))
                /(250*sqrt(pow(y + (7*cos(t*pi))/25 - 0.5,2) + pow((7*sin(t*pi))/25 - x + 0.5,2)));

    }

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return M_PI*(0.5-P.y);
        else return M_PI*(P.x-0.5);
    }

    // Normal x-direction
    R n1(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.x - xc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Normal y-direction
    R n2(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.y - yc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y);
    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return (M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*cos(M_PI*y))/125
        - (4*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(2*M_PI*t))/5
        - (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*sin(M_PI*y)*(x - 0.5))/5
        + (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*y)*sin(M_PI*x)*(y - 0.5))/5;

    }
}

namespace Example1_Omega1 {

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) - 0.17;
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return (pi*cos(pi*y)*sin(pi*x)*(2*cos(pi*t)*cos(pi*t) - 1)*((7*sin(pi*t))/25 - x + 0.5))
                /(250*sqrt(pow(y + (7*cos(t*pi))/25 - 0.5,2) + pow((7*sin(t*pi))/25 - x + 0.5,2)))
                - (pi*cos(pi*x)*sin(pi*y)*(2*cos(pi*t)*cos(pi*t) - 1)*(y + (7*cos(pi*t))/25 - 0.5))
                /(250*sqrt(pow(y + (7*cos(t*pi))/25 - 0.5,2) + pow((7*sin(t*pi))/25 - x + 0.5,2)));

    }

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return M_PI*(0.5-P.y);
        else return M_PI*(P.x-0.5);
    }

    // Normal x-direction
    R n1(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.x - xc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Normal y-direction
    R n2(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.y - yc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y);
    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return (M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*cos(M_PI*y))/125
        - (4*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(2*M_PI*t))/5
        - (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*sin(M_PI*y)*(x - 0.5))/5
        + (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*y)*sin(M_PI*x)*(y - 0.5))/5;

    }
}

namespace Lehrenfeld {

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        double r0 = 1.;
        double x = P.x, y = P.y;

        return -(sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y) - r0);
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        double r0 = 1.;
        return -(sqrt(P.x*P.x + P.y*P.y) - r0);
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return 0.;

    }

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return 1-P.y*P.y;
        else return 0.;
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.;
    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        double r0 = 1., x = P.x, y = P.y;
        return cos(M_PI*sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/r0)*sin(M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        double r0 = 1., x = P.x, y = P.y;
        return cos(M_PI*sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/r0)*sin(M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return M_PI*cos(M_PI*t)*cos(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            + (M_PI*sin(M_PI*t)*sin(M_PI*sqrtf((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)))/
            sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)
            + (M_PI*M_PI*sin(M_PI*t)*cos(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            *(2*x + 2*t*(y*y - 1))*(2*x + 2*t*(y*y - 1)))/(4*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            + (M_PI*M_PI*sin(M_PI*t)*cos(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            *(2*y + 4*t*y*(x + t*(y*y - 1)))*(2*y + 4*t*y*(x + t*(y*y - 1))))
            /(4*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            + (M_PI*sin(M_PI*t)*sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            *(8*t*t*y*y + 4*t*(x + t*(y*y - 1)) + 2))
            /(2*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            - (M_PI*sin(M_PI*t)*sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            *(2*x + 2*t*(y*y - 1))*(2*x + 2*t*(y*y - 1)))/(4*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)
            *sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            - (M_PI*sin(M_PI*t)*sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            *(2*y + 4*t*y*(x + t*(y*y - 1)))*(2*y + 4*t*y*(x + t*(y*y - 1))))
            /(4*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)*sqrt((x + t*(y*y - 1))
            *(x + t*(y*y - 1)) + y*y)) - (M_PI*sin(M_PI*t)
            *sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(x + t*(y*y - 1)))
            /sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y) + (M_PI*sin(M_PI*t)
            *sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(2*x + 2*t*(y*y - 1)))
            /(2*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y));

    }

}

namespace Convection {
    /* This works for running Test – i.e. a pure bulk problem on Omega_2. */

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return -sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) + 0.17;
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        return -sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) + 0.17;
    }

    // Velocity field
    R fun_velocity(const R2 P, const int i) {
        if (i == 0) return M_PI*(0.5-P.y);
        else return M_PI*(P.x-0.5);
    }

    // Normal x-direction
    R n1(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.x - xc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Normal y-direction
    R n2(const R2 P, const R t) {
        R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
        return (P.y - yc)/(sqrt((P.y-yc)*(P.y-yc) + (P.x-xc)*(P.x-xc)));
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y);
    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        return 0.5 + 0.4*cos(M_PI*P.x)*cos(M_PI*P.y)*cos(2*M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*y)*sin(M_PI*x)*(y - 0.5))/5
        - (2*M_PI*M_PI*cos(2*M_PI*t)*cos(M_PI*x)*sin(M_PI*y)*(x - 0.5))/5
        - (4*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(2*M_PI*t))/5;

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
// Set numerical example (options: "circle", "flower", "example1", "lehrenfeld")
#define example1
// Set boundary condition type on Omega 2 (options: "dirichlet", "neumann" – note: neumann only works combined with example1)
#define neumann
// Set scheme for the dg method (options: "schemeII", "schemeIII", "reynold" see thesis. Irrelevant if "cg" is defined instead of "dg")
#define reynold
// Set stabilization method (options: "fullstab", "macro")
#define fullstab
// Decide whether to solve for level set function, or to use exact (options: "levelsetsolve", "levelsetexact")
#define levelsetexact
// Solve on Omega_1 (options: "omega1" or anything else to solve on Omega 2)
#define omega2
// If "omega1" is defined, set type of BCs on outer boundary (options: "dirichlet1" or "neumann1")
#define dirichlet1



#ifdef circle
    using namespace Circle;
#elif defined(flower)
    using namespace Flower;
#elif defined(example1)
    #ifdef omega1
        using namespace Example1_Omega1;
    #else
        using namespace Example1;   // on Omega 2
    #endif
#elif defined(lehrenfeld)
    using namespace Lehrenfeld;   // on Omega 2
#elif defined(convection)
    using namespace Convection;
#endif

int main(int argc, char** argv) {

    // Mesh settings and data objects
    const size_t iterations = 1;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 10, ny = 10;       // starting mesh size

#ifdef example1
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Example1/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Example1/paraview/";

#elif defined(lehrenfeld)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Lehrenfeld/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Lehrenfeld/paraview/";
#elif defined(convection)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Convection/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Convection/paraview/";
#else
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Other/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Other/paraview/";
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
    #if defined(example1) || defined(convection)
        const double lx = 1., ly = 1.;
        Mesh Th(nx, ny, 0., 0., lx, ly);
    #elif defined(lehrenfeld)
        const double lx = 7., ly = 3.;
        Mesh Th(nx, ny, -3.5, -1.5, lx, ly);
    #else
        const double lx = 3., ly = 3.;
        Mesh Th(nx, ny, -1.5, -1.5, lx, ly);
    #endif

    //// Parameters

        // Mesh size
        double h = lx/(nx);
        //double h = sqrt(lx*lx/(nx*nx) + ly*ly/(ny*ny));
        hs.at(j) = h;
        int divisionMeshSize = 3;

        // Time
        double dT = h/divisionMeshSize; // Time step size
        double tfinal = .45;            // Final time
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

    #if defined(flower) || defined(circle) // For Circle and Flower domains, we use different problem parameters
        // General problem parameters
        double A2 = 0.5;        // diffusion coefficients
        double kappa2 = 0.5, kappa02 = 2;
        double kappaTilde2 = kappa2/kappa02;

        // Constants for penalty terms
        double tau_a2 = 300;        // diffusion penalty scaling
        double tau_b2 = 0;        // convection penalty scaling
        // Bulk penalties
        double lambdaA = kappaTilde2*A2*tau_a2/h;      // diffusion term
        double lambdaB = kappaTilde2*A2*tau_b2;        // convection term

        // Penalty parameter for outer boundary
        double lambda = A2*10/h;    // coefficient on the outer boundary

        #ifdef dg
        // DG stabilization parameters
        double tau20 = A2, tau21 = 0.1*A2;      // bulk
        #elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 0.1*A2;
        #endif
    #elif defined(lehrenfeld)
        double A2 = 1;
        double kappaTilde2 = 1;

        // Constants for penalty terms
        double tau_a2 = 500;        // diffusion penalty scaling
        double tau_b2 = 10;        // convection penalty scaling

        // Bulk penalties
        double lambdaA = tau_a2/h;      // diffusion term
        double lambdaB = tau_b2;            // convection term

        // Penalty parameter for outer boundary
        double lambda = 500/h;    // only used when Dirichlet BCs apply

        #ifdef dg
        // DG stabilization parameters
        double tau20 = 1e-1, tau21 = 1e-1;      // bulk
        #elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 1e-1;
        #endif
    #elif defined(convection)
        double A2 = 0;
        double kappaTilde2 = 1;

        // Constants for penalty terms
        double tau_a2 = 500;        // diffusion penalty scaling
        double tau_b2 = 10;        // convection penalty scaling

        // Bulk penalties
        double lambdaA = tau_a2*A2/h;      // diffusion term
        double lambdaB = tau_b2;            // convection term

        // Penalty parameter for outer boundary
        double lambda = 10;    // only used when Dirichlet BCs apply

        #ifdef dg
        // DG stabilization parameters
        double tau20 = 1e-1, tau21 = 1e-1;      // bulk
        #elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 1e-1;
        #endif
    #else   // For Example 1
        double A2 = 0.01;
        double kappaTilde2 = 1;

        // Constants for penalty terms
        double tau_a2 = 500;        // diffusion penalty scaling
        double tau_b2 = 10;        // convection penalty scaling

        // Bulk penalties
        double lambdaA = tau_a2/h;      // diffusion term
        double lambdaB = tau_b2;            // convection term

        // Penalty parameter for outer boundary
        double lambda = 500/h;    // coefficient on the outer boundary

        #ifdef dg
        // DG stabilization parameters
        double tau20 = 1e-1, tau21 = 1e-1;      // bulk
        #elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 1e-1;
        #endif
    #endif

        // =============================================================== //
        // ------------------------ Spaces ------------------------------- //
        // =============================================================== //

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space

    #ifdef dg
        FESpace2 Vh(Th, DataFE<Mesh>::P1dc);        // discontinuous basis functions
    #elif defined(cg)
        FESpace2 Vh(Th, DataFE<Mesh>::P1);          // continuous basis functions
    #endif
        // 1D Time mesh
        Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        // Quadrature data
        const QuadratureFormular1d& qTime(*Lobatto(3));
        const Uint nbTime = qTime.n;
        const Uint ndfTime = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Velocity field
        Lagrange2 FEvelocity(2);
        FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);

        double dt_levelSet = dT/(nbTime-1);
        vector<Fun_h> ls(nbTime);
    #if defined(levelsetexact)
        for (int i=0; i<nbTime; i++) ls[i].init(Lh, fun_levelSet, 0.);
    #elif defined(levelsetsolve)
        for (int i=0; i<nbTime; i++) ls[i].init(Lh, fun_levelSet);
    #endif

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        // Bulk-Surface Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;

        int iter = 0;
        double q0_0, q0_1, qp_0, qp_1;
        double intF = 0, intG = 0;      // hold integrals of rhs and Neumann bcs

        // Iterate over time-slabs
        while (iter < GTime::total_number_iteration) {


            GTime::current_iteration = iter;
            const TimeSlab& In(Ih[iter]);

            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << GTime::total_number_iteration << std::endl;
            std::cout << " Time      \t : \t" << GTime::current_time() << std::endl;

            ls.begin()->swap(ls[nbTime - 1]);

            // computation of the interface in each of the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

            #if defined(levelsetexact)
                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);
            #endif
                // FIXME: Below gives a zsh trace trap error after in the second time step
                interface.init(i, Th, ls[i]);

                // We solve for the level-set using Crank-Nicholson in time
            #if defined(levelsetsolve)
                if (i<lastQuadTime) {
                    LevelSet::move((ls[i], vel, vel, dt_levelSet, ls[i+1]));
                }
            #endif

          }

            // Create active meshes
            ActiveMesh<Mesh> Kh2(Th);
            Kh2.truncate(interface, 1);


            // Cut FE space
            CutSpace Wh(Kh2, Vh);
            Wh.info();
            convdiff.initSpace(Wh, In);

            // gnuplot::save(Th, "Th.dat");
            // gnuplot::save(cutThTime,  "cutThTime.dat");
            // gnuplot::save(*interface[0],  "interface.dat");

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Wh, In, fun_rhsBulk);
            Fun_h g(Wh, In, fun_uBulk);    // create an FE-function of the exact bulk solution Omega1

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1);

            // Data for initial solution
            Rn data_u0;                           // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));        // initial data bulk

            if (iter == 0) interpolate(Wh, data_B0, fun_uBulkInit);
            Rn uh(data_u0);

            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);

            // Plot initial solution in paraview
            if (iter == 0 && MPIcf::IamMaster()) {
                Paraview<Mesh> writerInitial(Kh2, pathOutputFigures + "BulkInitial.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);

                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);

                writerInitial.add(uBex, "bulk_exact", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);

            }


        //// Assembling linear and bilinear forms

        // Partial derivative in time, bulk

        // reynold scheme
        #ifdef reynold


            convdiff.addBilinear(
                - innerProduct(u, kappaTilde2*dt(v))
                , Kh2
                , In
            );

            double cc0 = 1./In.T.mesure()*1./qTime[0].a;

            convdiff.addBilinear(
                + innerProduct(u, cc0*kappaTilde2*v)
                , Kh2
                , (int)lastQuadTime
                , In
            );

            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), cc0*kappaTilde2*v)
                , Kh2
                , 0
                , In
            );

        // classical scheme
        #else

            convdiff.addBilinear(
                + innerProduct(dt(u), kappaTilde2*v)
                , In
            );

            // Time penalty term bulk LHS
            double cc0 = 1./In.T.mesure()*1./qTime[0].a;

            convdiff.addBilinear(
                + innerProduct(u, cc0*kappaTilde2*v)
                , 0
                , In
            );

            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), cc0*kappaTilde2*v)
                , 0
                , In
            );
        #endif

        //// Scheme for diffusion

            convdiff.addBilinear(
                + innerProduct(kappaTilde2*A2*grad(u), grad(v))
                , Kh2
                , In
            );

        #ifdef dg
            // integral on inner edges (E_h)
            // convdiff.addBilinear(
            //     - innerProduct(kappaTilde2*A2*average(grad(u)*n), jump(v))
            //     - innerProduct(kappaTilde2*A2*jump(u), average(grad(v)*n))
            //     + innerProduct(lambdaA*jump(u), jump(v))
            //     , Kh2
            //     , INTEGRAL_INNER_EDGE_2D
            //     , In
            // );

        #endif

        //// Schemes for convection

        #if defined(schemeII) && defined(dg)
            // Bulk convection
            convdiff.addBilinear(
                + innerProduct((vel.expression(),grad(u)), kappaTilde2*v)*0.5
                - innerProduct(kappaTilde2*u, (vel.expression(),grad(v)))*0.5
                , In
            );

            // Inner edges bulk convection
            convdiff.addEdgeIntegral(
                + innerProduct(average((vel*n)*u), kappaTilde2*jump(v))*0.5
                - innerProduct(jump((vel*n)*u), kappaTilde2*average(v))*0.5
                + innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))
                , In
            );

        #elif defined(schemeIII) && defined(dg)
            // Bulk convection
            convdiff.addBilinear(
                - innerProduct(kappaTilde2*u, (vel.expression(),grad(v)))
                , In
            );

            // Inner edges bulk convection
            convdiff.addEdgeIntegral(
                + innerProduct(average((vel*n)*u), kappaTilde2*jump(v))
                + innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))
                , In
            );

        #elif defined(reynold)

            // Convection term
            convdiff.addBilinear(
                - innerProduct(kappaTilde2*u, (vel.expression()*grad(v)))
                , Kh2
                , In
            );


            // Added terms

            convdiff.addBilinear(
              innerProduct(jump(u), jump(v))
                // + innerProduct(average(vel*n*u), jump(v))
                // + innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))
                , Kh2
                , INTEGRAL_INNER_EDGE_2D
                , In
                , 1
            );
            getchar();



        #else   // classic CG scheme
            convdiff.addBilinear(
                + innerProduct((vel.expression(),grad(u)), kappaTilde2*v)
                , In
            );
        #endif


        // Stabilization

        #ifdef macro
            TimeMacroElement TimeMacro(Wh, In, qTime, 0.125);

            // gnuplot::save(Th, "Th.dat");
            // gnuplot::save(Wh,0, "Vh1.dat");
            // gnuplot::save(cutThTime, "cutThTime.dat");
            // gnuplot::save(*interface[0], "interface0.dat");
            // gnuplot::save(*interface[1], "interface1.dat");
            // gnuplot::save(*interface[2], "interface2.dat");
            // gnuplot::save(TimeMacro);

            // Stabilization of the bulk
            convdiff.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), kappaTilde2*jump(v))
                + innerProduct(h*tau21*jump(grad(u)), kappaTilde2*jump(grad(v)))
                , In
                , TimeMacro
            );

        #elif defined(fullstab)
            convdiff.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), kappaTilde2*jump(v))
                + innerProduct(h*tau21*jump(grad(u)), kappaTilde2*jump(grad(v)))
                , Kh2
                , In
            );
        #endif


        // Boundary conditions on interface

    #ifdef neumann
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
            convdiff.addLinear(
                + innerProduct(g_Neumann.expression(), kappaTilde2*v)
                , interface
                , In
            );
        #ifdef schemeII
            convdiff.addBilinear(
                + innerProduct((vel*n)*u, kappaTilde2*v)*0.5
                , interface
                , In
            );
        #elif defined(schemeIII)
            convdiff.addBilinear(
                + innerProduct((vel*n)*u, kappaTilde2*v)
                , interface
                , In
            );
        #endif

    #elif defined(dirichlet)
            convdiff.addBilinear(
                - innerProduct(kappaTilde2*A2*grad(u)*n, v)   // from IBP
                - innerProduct(u, kappaTilde2*A2*grad(v)*n)   // added to make symmetric
                + innerProduct(u, kappaTilde2*lambda*v)       // added penalty
                , interface
                , In
            );


            // RHS on Gamma
        #if defined(schemeII) && defined(dg)
            convdiff.addLinear(
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)*0.5
                , interface
                , In
            );

        #elif defined(schemeIII) && defined(dg)
            convdiff.addLinear(
                - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                , interface
                , In
            );

        #elif defined(reynold) || defined(cg)   // cg reynold, dg reynold and cg classic
            convdiff.addLinear(
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                , interface
                , In
            );

        #endif

    #endif

    // Boundary conditions on outer boundary
    #if defined(omega1) && defined(example1)

        #ifdef dirichlet1
        // solve with Dirichlet BCs

            convdiff.addBilinearFormBorder(
                - innerProduct(kappaTilde2*A2*grad(u)*n, v)   // from IBP
                - innerProduct(u, kappaTilde2*A2*grad(v)*n)   // added to make symmetric
                + innerProduct(u, kappaTilde2*lambda*v)       // added penalty
                , In
            );

               // RHS on external boundary
        #if defined(schemeII) && defined(dg)
            convdiff.addLinearFormBorder(
                - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)*0.5
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                , In
            );

        #elif defined(schemeIII) && defined(dg)
            // NOTE: This gives errors
            convdiff.addLinearFormBorder(
                - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                , In
            );

        #elif defined(reynold) || defined(cg)   // cg reynold, dg reynold and cg classic
            convdiff.addLinearFormBorder(
                - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)
                + innerProduct(g.expression(), kappaTilde2*lambda*v)
                , In
            );

        #endif

        #elif defined(neumann1)
        // solve with Neumann BCs
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
        #ifdef schemeII
            convdiff.addBilinearFormBorder(
                + innerProduct((vel*n)*u, kappaTilde2*v)*0.5
                , In
            );
        #elif defined(schemeIII)
            convdiff.addBilinearFormBorder(
                + innerProduct((vel*n)*u, kappaTilde2*v)
                , In
            );
        #endif

            convdiff.addLinearFormBorder(
                + innerProduct(g_Neumann.expression(), kappaTilde2*v)
                , In
            );
        #endif

    #endif

            // Add RHS on bulk
            convdiff.addLinear(
                + innerProduct(f.expression(), kappaTilde2*v)
                , Kh2
                , In
            );

            // Compute integrals
            intF = integral(Kh2, In, f, 0, qTime);
        #ifdef neumann
            intG = integral(g_Neumann, In, interface, 0);
        #elif defined(dirichlet)
            double intu = integralSurf(b0h, In, qTime);
            double intg = integralSurf(g, In, qTime);
            double intDiff = lambda*(intu-intg);
            double intGrad = integral<Mesh>((A2, grad(u)*n), b0h, In, qTime);

            //std::cout << "u - g = " << intu-intg << std::endl;
            //std::cout << "intDiff = " << intDiff << std::endl;
        #endif



            // Solve linear system
            convdiff.solve("mumps");

            uh = convdiff.rhs_;
            convdiff.saveSolution(uh);

            // Compute conservation error
            if (iterations == 1) {

                Fun_h funuh(Wh, data_u0);

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_u0(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Kh2, funsol, 0, 0);
                sol2 += data_u0(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Kh2, funsol, 0, lastQuadTime);

                if(iter==0) {
                    q0_0 = q_0; q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Kh2, b0h, 0, 0);
                }
                // std::cout << q_0 - q0_1 << std::endl;
                // std::cout << q_0 << std::endl;
                // std::cout << q0_1 << std::endl;
                // getchar();

                // cout << " q0 - qt_1 \t" << q_1-qp_1 << std::endl;
                // cout << " int F \t" << intF << std::endl;

                outputData << setprecision(10);
                outputData << GTime::current_time() << ","
                        << (q_1-qp_1) << ","
                        << intF << ","
                        << intG << ","
                    #ifdef neumann
                        << ((q_1 -qp_1) - intF - intG) << std::endl;
                    #elif defined(dirichlet)
                        << ((q_1 -qp_1) - intF + intGrad + intDiff) << "," << std::endl;
                    #endif
                        //    << q0_1 - q_1 <<  std::endl;
                qp_1 = q_1;

            }


            R errBulk =  L2normCut(b0h, fun_uBulkD, GTime::current_time(), 0, 1);
            std::cout << std::endl;
            std::cout << " L2 Error \t : \t" << errBulk << std::endl;

            //gnuplot::save(cutThTime,"cutThTime.dat");
            //gnuplot::save(Th, "Th.dat");

            errorsBulk.at(j) = errBulk;

            if ((iterations == 1) && MPIcf::IamMaster()) {

                Paraview<Mesh> writer(Kh2, pathOutputFigures + "Bulk" + to_string(iter + 1)+"DG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, GTime::current_time());
                Fun_h fB(Wh, fun_rhsBulk, GTime::current_time());
                writer.add(uBex,"bulk_exact", 0, 1);
                writer.add(fB,"bulk_rhs", 0, 1);
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);
            }


            iter++;

        }

        // Refine mesh
        nx *= 2;
        ny *= 2;

    }

    std::cout << std::endl;
    std::cout << "Errors Bulk = [";
    for (int i=0; i<iterations; i++) {

        std::cout << errorsBulk.at(i);
        if (i < iterations-1) {
            std::cout << ", ";
        }

    }
    std::cout << "]" << std::endl;

    std::cout << "h = [";
    for (int i=0; i<iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations-1) {
            std::cout << ", ";
        }

    }
    std::cout << "]" << std::endl;

    return 0;
}