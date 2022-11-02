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

    Classical : Integration by parts on full convection term, no term added to make anti-symmetric
    Conservative: Reynold's transport theorem is used to make the bilinear form fulfill a conservation law
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
    R fun_velocity(const R2 P, const int i, const R t) {
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
        //return cos(M_PI*((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/(r0*r0))*sin(M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        double r0 = 1., x = P.x, y = P.y;
        return cos(M_PI*sqrt((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/r0)*sin(M_PI*t);
        //return cos(M_PI*((x - (1-y*y)*t)*(x - (1-y*y)*t) + y*y)/(r0*r0))*sin(M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return M_PI*cos(M_PI*t)*cos(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))
            + (M_PI*sin(M_PI*t)*sin(M_PI*sqrt((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)))/
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

        //return M_PI*cos(M_PI*t)*cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)) + 2*M_PI*sin(M_PI*t)*sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y)) + M_PI*sin(M_PI*t)*sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(8*t*t*y*y + 4*t*(x + t*(y*y - 1)) + 2) + M_PI*M_PI*cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*sin(M_PI*t)*(2*x + 2*t*(y*y - 1))*(2*x + 2*t*(y*y - 1)) + M_PI*M_PI*cos(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*sin(M_PI*t)*(2*y + 4*t*y*(x + t*(y*y - 1)))*(2*y + 4*t*y*(x + t*(y*y - 1))) + M_PI*sin(M_PI*t)*sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(2*x + 2*t*(y*y - 1)) - 2*M_PI*sin(M_PI*t)*sin(M_PI*((x + t*(y*y - 1))*(x + t*(y*y - 1)) + y*y))*(y*y - 1)*(x + t*(y*y - 1));

    }
    
}

namespace Lehrenfeld_6_2 {

    // Level-set function
    R fun_levelSet(const R2 P, const int i, const R t) {
        double r0 = .5 + Epsilon;
        double x = P.x, y = P.y;
        
        //double rho(R x, R y, R t) return 1./M_PI*sin(2*M_PI*t);
        //double r(R x, R y, R t) return sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y);

        return -(sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y) - r0);
    }

    // Level-set function initial
    R fun_levelSet(const R2 P, const int i) {
        double r0 = .5 + Epsilon;
        double x = P.x, y = P.y;
        
        return -(sqrt(x*x + y*y) - r0);
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
    
        return 0.;
        
    }

    // Velocity field
    R fun_velocity(const R2 P, const int i, const R t) {
        if (i == 0) return 2*cos(2*M_PI*t);
        else return 0.;
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.;
    }
}

namespace Lehrenfeld_6_2_Convection_Dominated {

    // Level-set function
    R fun_levelSet(const R2 P, const int i, const R t) {
        double r0 = .5 + Epsilon;
        double x = P.x, y = P.y;
        
        //double rho(R x, R y, R t) return 1./M_PI*sin(2*M_PI*t);
        //double r(R x, R y, R t) return sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y);

        return -(sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y) - r0);
    }

    // Level-set function initial
    R fun_levelSet(const R2 P, const int i) {
        double r0 = .5 + Epsilon;
        double x = P.x, y = P.y;
        
        return -(sqrt(x*x + y*y) - r0);
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
    
        return 0.;
        
    }


    // Velocity field
    R fun_velocity(const R2 P, const int i, const R t) {
        if (i == 0) return 2*cos(2*M_PI*t);
        else return 0.;
    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        return 0.;
    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        double r0 = .5, x = P.x, y = P.y;
        
        return cos(M_PI*((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y)/(r0*r0))*sin(M_PI*t);
        //return cos(M_PI*sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y)/r0)*sin(M_PI*t);
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        double r0 = .5, x = P.x, y = P.y;
        
        return cos(M_PI*((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y)/(r0*r0))*sin(M_PI*t);
        //return cos(M_PI*sqrt((x-1./M_PI*sin(2*M_PI*t))*(x-1./M_PI*sin(2*M_PI*t)) + y*y)/r0)*sin(M_PI*t);
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
        return 3.1415927*cos(3.1415927*t)*cos(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y) + 0.50265482*sin(3.1415927*t)*sin(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y) + 6.3165468*y*y*sin(3.1415927*t)*cos(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y) + 1.5791367*sin(3.1415927*t)*cos(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y)*(2.0*x - 0.63661977*sin(6.2831853*t))*(2.0*x - 0.63661977*sin(6.2831853*t)) - 25.132741*cos(6.2831853*t)*sin(3.1415927*t)*sin(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y)*(2.0*x - 0.63661977*sin(6.2831853*t)) + 50.265482*cos(6.2831853*t)*sin(3.1415927*t)*sin(12.566371*(x - 0.31830989*sin(6.2831853*t))*(x - 0.31830989*sin(6.2831853*t)) + 12.566371*y*y)*(x - 0.31830989*sin(6.2831853*t));
    }


    
}

namespace Lehrenfeld_6_3 {

    // Level-set function
    double fun_levelSet(const R2 P, const int i, const R t) {
        double x = P.x, y = P.y;

        return -(min(sqrt(x*x + (y-t+3./4)*(y-t+3./4)), sqrt(x*x + (y-t+3./4)*(y+t-3./4))) - 0.5);
    }

    // Level-set function initial
    double fun_levelSet(const R2 P, const int i) {
        
        double x = P.x, y = P.y;

        return -(min(sqrt(x*x + (y+3./4)*(y+3./4)), sqrt(x*x + (y+3./4)*(y-3./4))) - 0.5);
    }

    // The rhs Neumann boundary condition
    R fun_neumann_Gamma(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;
    
        return 0.;
        
    }

    // Velocity field
    R fun_velocity(const R2 P, const int i, const R t) {
        
        if ((P.y > 0 && t <= 3./4) || (P.y < 0 && t > 3./4)) {
            if (i == 0) {
                return 0.;
            }
            else return -1.;
        }

        else if ((P.y <= 0 && t <= 3./4) || (P.y > 0 && t > 3./4)) {
            if (i == 0) {
                return 0.;
            }
            else return 1.;
        }

    }

    // Initial solution bulk
    R fun_uBulkInit(const R2 P, const int i) {
        
        // Sign(y)
        if (P.y > 0) return 1.;
        else if (P.y == 0) return 0.;
        else return -1.;

    }

    // Exact solution bulk
    R fun_uBulk(const R2 P, const int i, const R t) {
        return 0.;
    }

    R fun_uBulkD(const R2 P, const int i, const int d, const R t) {
        return 0.;
    }

    // RHS fB bulk
    R fun_rhsBulk(const R2 P, const int i, const R t) {
        R x = P.x, y = P.y;

        return 0.;

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


// Choose Discontinuous or Continuous Galerkin method (options: "dg", "cg")
#define dg
// Set numerical example (options: "lehrenfeld_6_2", "lehrenfeld_6_3")
#define lehrenfeld_6_2

#define convection_dominated_not
// Set boundary condition type on Omega 2 (options: "dirichlet", "neumann" – note: neumann only works combined with example1)
#define neumann
// Set scheme for the dg method (options: "conservative", "classical" see thesis. Irrelevant if "cg" is defined instead of "dg")
#define classical  
// Set stabilization method (options: "fullstab", "macro") 
#define fullstab       
// Decide whether to solve for level set function, or to use exact (options: "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h
#define use_t

#if defined(lehrenfeld_6_2)
    #ifdef convection_dominated
    using namespace Lehrenfeld_6_2_Convection_Dominated;
    #else
    using namespace Lehrenfeld_6_2;
    #endif
#elif defined(lehrenfeld_6_3)
    using namespace Lehrenfeld_6_3;
#endif

int main(int argc, char** argv) {
    
    // Mesh settings and data objects
    const size_t iterations = 5;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 10, ny = 10;       // starting mesh size
    double h = 0.05;             // starting mesh size
    double dT = 0.25;
#if defined(lehrenfeld_6_2)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Lehrenfeld_6_2/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Lehrenfeld_6_2/paraview/";
#elif defined(lehrenfeld_6_3)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeBulk/Lehrenfeld_6_3/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeBulk/Lehrenfeld_6_3/paraview/";
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
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {
        
    #ifdef lehrenfeld_6_2
        const double lx = 2., ly = 1.2; 
        #ifdef use_h
        nx = (int)(lx/h)+1, ny = (int)(ly/h)+1;
        #elif defined(use_n)
        h = lx/(nx-1);
        #endif
        Mesh Th(nx, ny, -1., -0.6, lx, ly);
    #elif defined(lehrenfeld_6_3)
        const double lx = 2.7, ly = 2.7;
        #ifdef use_h
        nx = (int)(lx/h)+1, ny = (int)(ly/h)+1;
        #elif defined(use_n)
        h = lx/(nx-1);
        #endif
        Mesh Th(nx, ny, -1.35, -1.35, lx, ly);
    #endif

        // Mesh size
        int divisionMeshSize = 2;
        //double dT = h/divisionMeshSize;
        hs.at(j) = h;
        dts.at(j) = dT;
        
        double tfinal = 0.25;            // Final time
    #ifdef use_t
        GTime::total_number_iteration = ceil(tfinal/dT);
    #else
        GTime::total_number_iteration = int(tfinal/dT);
    #endif
        dT = tfinal / GTime::total_number_iteration;
        GTime::time_step = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "Iteration " << j + 1 << "/" << iterations << std::endl;
        }
        std::cout << "h = " << h << std::endl;
        std::cout << "nx = " << nx << std::endl;
        std::cout << "nx = " << ny << std::endl;
        std::cout << "dT = " << dT << std::endl;

        #ifdef convection_dominated
            double A2 = 0.01;
        #else 
            double A2 = 1;
        #endif
        
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
        double tau20 = 5e-2, tau21 = 5e-2;      // bulk
        #elif defined(cg)
        // CG stabilization parameters
        double tau20 = 0, tau21 = 1e-1;
        #endif        


        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        FESpace2 Vh2(Th, DataFE<Mesh>::P2);
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
        vector<Fun_h> vel(nbTime);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        double dt_levelSet = dT/(nbTime-1);
        vector<Fun_h> ls(nbTime);

        for (int i=0; i<nbTime; i++) ls[i].init(Lh, fun_levelSet, 0.);

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);
        
        // Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;
        
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
    
            ls.begin()->swap(ls[nbTime - 1]);

            // computation of the interface in each of the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);

                vel[i].init(VelVh, fun_velocity, tt);

                interface.init(i, Th, ls[i]);
            }

            // Create active meshes
            ActiveMesh<Mesh> Kh2(Th);
            Kh2.truncate(interface, -1);

            // Cut FE space
            CutSpace Wh(Kh2, Vh); 
            
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Vh2, In, fun_rhsBulk);

            Fun_h g(Vh2, In, fun_uBulk);    // create an FE-function of the exact bulk solution Omega1

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega2                 
            
            // Data for initial solution
            Rn data_u0;                           // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));        // initial data bulk
            
            if (iter == 0) interpolate(Wh, data_B0, fun_uBulkInit);
            
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
        #ifdef conservative  

            convdiff.addBilinear(
                - innerProduct(u, dt(v))
                , Kh2
                , In
            );

            convdiff.addBilinear(
                + innerProduct(u, v)
                , Kh2
                , (int)lastQuadTime
                , In
            );

            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), v)
                , Kh2
                , 0
                , In
            );

        // classical scheme
        #elif defined(classical)   

            convdiff.addBilinear(
                + innerProduct(dt(u), v)
                , Kh2
                , In
            );

            convdiff.addBilinear(
                + innerProduct(u, v)
                , Kh2
                , 0
                , In
            );

            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), v)
                , Kh2
                , 0
                , In
            );

        #endif
            
        //// Scheme for diffusion

            convdiff.addBilinear(
                + innerProduct(A2*grad(u), grad(v))
                , Kh2
                , In
            );
            
        #ifdef dg
            // integral on inner edges (E_h)
            convdiff.addBilinear(
                - innerProduct(A2*average(grad(u)*n), jump(v))      
                - innerProduct(A2*jump(u), average(grad(v)*n))     
                + innerProduct(lambdaA*jump(u), jump(v))          
                , Kh2
                , INTEGRAL_INNER_EDGE_2D
                , In
            );
        #endif
            
        //// Schemes for convection

        #if defined(dg) && (defined(conservative) || defined(classical))

            // Convection terms

            for (int i=0; i<nbTime; ++i) {
                
                // Convection term
                convdiff.addBilinear(
                - innerProduct(u, (vel[i].expression()*grad(v)))
                , Kh2
                , In
                , i
                );

                // Flux terms 
                convdiff.addBilinear(
                + innerProduct(average(vel[i]*n*u), jump(v))
                //+ innerProduct(lambdaB*fabs(vel[i]*n)*jump(u), jump(v))
                + innerProduct(0.5*fabs(vel[i]*n)*jump(u), jump(v))
                , Kh2
                , INTEGRAL_INNER_EDGE_2D
                , In
                , i
                );
            }

        #endif

        // Stabilization

        #ifdef macro    
            TimeMacroElement2<Mesh> TimeMacro(Kh2, qTime, 0.125);

            // Stabilization of the bulk 
            convdiff.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), jump(v))
                + innerProduct(h*tau21*jump(grad(u)), jump(grad(v)))
                , Kh2
                , In
                , TimeMacro
            );

        #elif defined(fullstab)
            convdiff.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), jump(v))
                + innerProduct(h*tau21*jump(grad(u)), jump(grad(v)))
                , Kh2
                , In
            );
        #endif
            
               
            // Boundary conditions on interface
            Fun_h g_Neumann(Wh, In, fun_neumann_Gamma);
            convdiff.addLinear(
                + innerProduct(g_Neumann.expression(), v)
                , interface
                , In
            );

            #if defined(classical)
            
            for (int i=0; i<nbTime; ++i) {
                convdiff.addBilinear(
                    + innerProduct((vel[i]*n)*u, v)
                    , interface
                    , In
                    , i
                );
            }
            #endif

            // Add RHS on bulk
            convdiff.addLinear(
                + innerProduct(f.expression(), v)
                , Kh2
                , In
            );

            // Compute integrals
            intF = integral(Kh2, In, f, 0, qTime);
            intG = integral(g_Neumann, In, interface, 0);
            
            if (iter == GTime::total_number_iteration-1) matlab::Export(convdiff.mat_, pathOutputFolder + "mat_h" + to_string(h) + "_" + to_string(j+1) + ".dat");
            // Solve linear system
            convdiff.solve();
            
            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

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
            
                outputData << setprecision(10);
                outputData << GTime::current_time() << ","
                        << (q_1-qp_1) << ","
                        << intF << ","
                        << intG << ","
                    #ifdef neumann
                        << ((q_1 -qp_1) - intF - intG) << std::endl;
                    #elif defined(dirichlet)
                        //<< ((q_1 -qp_1) - intF + intGrad + intDiff) << "," << std::endl;
                        << "," << std::endl;
                    #endif
                qp_1 = q_1;
                
            }

            {
                Rn sol(Wh.get_nb_dof(), 0.);
                sol += data_u0(SubArray(Wh.get_nb_dof(), 0));
                Fun_h funuh(Wh, sol);
                double errL2 =  L2normCut(funuh, fun_uBulkD,  GTime::current_time(), 0, 1);
                std::cout << " t_{n-1} -> || u-uex||_2 = " << errL2 << std::endl;

                sol  += data_u0(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
                errL2 =  L2normCut(funuh, fun_uBulkD, GTime::current_time()+dT, 0, 1);
                std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;
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
                
                writer.add(uBex, "bulk_exact", 0, 1);
                writer.add(fB,"bulk_rhs", 0, 1);
                writer.add(ls[0], "levelSet0", 0, 1);
                writer.add(ls[1], "levelSet1", 0, 1);
                writer.add(ls[2], "levelSet2", 0, 1);
                writer.add(vel[0], "velocity0", 0, 1);
                writer.add(vel[1], "velocity1", 0, 1);
                writer.add(vel[2], "velocity2", 0, 1);
            }

            if (iterations > 1 && iter == GTime::total_number_iteration-1) outputData << h << "," << dT << "," << errBulk << std::endl;

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
    
    std::cout << "dT = [";
    for (int i=0; i<iterations; i++) {

        std::cout << dts.at(i);
        if (i < iterations-1) {
            std::cout << ", ";
        }

    }
    std::cout << "]" << std::endl;
    return 0;
}

