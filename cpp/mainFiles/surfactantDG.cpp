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


namespace NumericSurfactantEllipse2D {

    R fun_init_surfactant(const R2 P, const int i) {return P.x*P.y ;} //2*... = 4*pi
    R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-4*t);}

    R fun_velocity(const R2 P, const int i, const R t) {
        R a = (1+0.25*sin(2*M_PI*t));
        R b = M_PI/4*(cos(2*M_PI*t))/a;
        return (i==0) ? b*P.x : 0;
    }

    R fun_levelSet(const R2 P, const int i, const R t) {
        R a = (1+0.25*sin(2*M_PI*t));
        // a = 1;
        return sqrt(P.x*P.x/a + P.y*P.y) - 1 - Epsilon;
    }

  R fun_rhs(const R2 P, const int cc, const R t) {
    R x = P.x,  y = P.y;
    R r = (4*x*y*exp(-4*t)*(sin(2*M_PI*t) + 4)*
           (2048*x*x*y*y + 336*y*y*y*y*pow(sin(2*M_PI*t),2) +
            52*y*y*y*y*pow(sin(2*M_PI*t),3) + 3*y*y*y*y*pow(sin(2*M_PI*t),4)
            + 64*x*x*x*x*sin(2*M_PI*t) + 960*y*y*y*y*sin(2*M_PI*t) + 1024*x*x*x*x
            + 1024*y*y*y*y + 144*x*x*y*y*pow(sin(2*M_PI*t),2)
            + 4*x*x*y*y*pow(sin(2*M_PI*t),3) + 1024*x*x*y*y*sin(2*M_PI*t)))
      /pow((y*y*pow(sin(2*t*M_PI),2) + 8*y*y*sin(2*t*M_PI) + 16*x*x + 16*y*y),3)
      - 4*x*y*exp(-4*t) + (x*y*M_PI*exp(-4*t)*cos(2*M_PI*t))/(sin(2*M_PI*t) + 4)
      + (x*y*y*y*M_PI*exp(-4*t)*cos(2*M_PI*t)*(sin(2*M_PI*t) + 4))/
      (y*y*pow(sin(2*M_PI*t),2) + 8*y*y*sin(2*M_PI*t) + 16*x*x + 16*y*y);

    return r;
  }
}

namespace NumericSurfactantEllipse2DVersion2 {
  R fun_init_surfactant(const R2 P, const int i) {return P.x*P.y ;} //2*... = 4*pi
  R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-t/4);}
  R fun_velocity(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    R b = M_PI/4*(cos(2*M_PI*t))/a;
    return (i==0) ? b*P.x : 0;
  }
  R fun_levelSet(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    // a = 1;
    return sqrt(P.x*P.x/a + P.y*P.y) - 1 - Epsilon;
  }
  R fun_rhs(const R2 P, const int cc, const R t) {
    R x = P.x,  y = P.y;
    // R r = (4*x*y*exp(-4*t)*(sin(2*M_PI*t) + 4)*
    //        (2048*x*x*y*y + 336*y*y*y*y*pow(sin(2*M_PI*t),2) +
    //         52*y*y*y*y*pow(sin(2*M_PI*t),3) + 3*y*y*y*y*pow(sin(2*M_PI*t),4)
    //         + 64*x*x*x*x*sin(2*M_PI*t) + 960*y*y*y*y*sin(2*M_PI*t) + 1024*x*x*x*x
    //         + 1024*y*y*y*y + 144*x*x*y*y*pow(sin(2*M_PI*t),2)
    //         + 4*x*x*y*y*pow(sin(2*M_PI*t),3) + 1024*x*x*y*y*sin(2*M_PI*t)))
    //   /pow((y*y*pow(sin(2*t*M_PI),2) + 8*y*y*sin(2*t*M_PI) + 16*x*x + 16*y*y),3)
    //   - 4*x*y*exp(-4*t) + (x*y*M_PI*exp(-4*t)*cos(2*M_PI*t))/(sin(2*M_PI*t) + 4)
    //   + (x*y*y*y*M_PI*exp(-4*t)*cos(2*M_PI*t)*(sin(2*M_PI*t) + 4))/
    //   (y*y*pow(sin(2*M_PI*t),2) + 8*y*y*sin(2*M_PI*t) + 16*x*x + 16*y*y);
    R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t)
    + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi)
    + 16*x*x + 16*y*y, 2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t))/(sin(2*pi*t) + 4)
    + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
    return r;
  }
  // R fun_rhs(const R2 P, const int cc, const R t) {
  //   R x = P.x,  y = P.y;
  //   R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t) + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi) + 16*x*x + 16*y*y,2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t)*(y*y*sin(2*pi*t) - 4*x*x + 4*y*y))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y) + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
  //   return r;
  // }
}

namespace Zahedi1 {

    R fun_init_surfactant(const R2 P, const int i) {return P.x*P.y + P.x*P.x*P.x*P.y*P.y;} 
    R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-4*t) + P.x*P.x*P.x*P.y*P.y;}
    R fun_velocity(const R2 P, const int i, const R t) {
        R a = (1+0.25*sin(2*M_PI*t));
        R b = M_PI/4*(cos(2*M_PI*t))/a;
        return (i==0) ? b*P.x : 0;
    }
    R fun_levelSet(const R2 P, const int i, const R t) {
        R a = (1+0.25*sin(2*M_PI*t));
        // a = 1;
        return sqrt(P.x*P.x/a + P.y*P.y) - 1 - Epsilon;
    }

    R fun_rhs(const R2 P, const int cc, const R t) {
        double x = P.x, y = P.y;

        return (exp(-4*t)*(1024*x*y*y*y - 512*x*x*x*x*x*x*x*exp(4*t) + 1024*x*x*x*y + 704*x*y*y*y*sin(2*pi*t) + 320*x*x*x*y*sin(2*pi*t) + 2816*x*x*x*y*y*y*y*exp(4*t) + 3840*x*x*x*x*x*y*y*exp(4*t) + 160*x*y*y*y*sin(2*pi*t)*sin(2*pi*t) + 16*x*x*x*y*sin(2*pi*t)*sin(2*pi*t) + 12*x*y*y*y*sin(2*pi*t)*sin(2*pi*t)*sin(2*pi*t) - 1536*x*y*y*y*y*y*y*exp(4*t) + 624*x*x*x*y*y*y*y*exp(4*t)*sin(2*pi*t)*sin(2*pi*t) + 56*x*x*x*y*y*y*y*exp(4*t)*sin(2*pi*t)*sin(2*pi*t)*sin(2*pi*t) - 1536*x*y*y*y*y*y*y*exp(4*t)*sin(2*pi*t) + 2304*x*x*x*y*y*y*y*exp(4*t)*sin(2*pi*t) + 960*x*x*x*x*x*y*y*exp(4*t)*sin(2*pi*t) - 576*x*y*y*y*y*y*y*exp(4*t)*sin(2*pi*t)*sin(2*pi*t) - 96*x*y*y*y*y*y*y*exp(4*t)*sin(2*pi*t)*sin(2*pi*t)*sin(2*pi*t) - 6*x*y*y*y*y*y*y*exp(4*t)*sin(2*pi*t)*sin(2*pi*t)*sin(2*pi*t)*sin(2*pi*t)))/pow((y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi) + 16*x*x + 16*y*y),2) - 4*x*y*exp(-4*t) + (x*pi*cos(2*pi*t)*(3*x*x*y*y + exp(-4*t)*y))/(sin(2*pi*t) + 4) + (y*y*pi*cos(2*pi*t)*(sin(2*pi*t) + 4)*(x*x*x*y*y + exp(-4*t)*x*y))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
    }
}

namespace Shi1 {

    R fun_init_surfactant(const R2 P, const int i) {return P.y/sqrt(P.x*P.x+P.y*P.y) + 2;} 
    R fun_sol_surfactant(const R2 P,  const int i, const R t) {
        R x = P.x, y = P.y;
        return exp(-t/4)*(y/sqrt((x-t)*(x-t) + y*y)) + 2;
    }
    R fun_velocity(const R2 P, const int i, const R t) {
        return (i==0) ? 1. : 0.;
    }
    R fun_levelSet(const R2 P, const int i, const R t) {
        return sqrt((P.x-t)*(P.x-t) + P.y*P.y) - 2 - Epsilon;
    }

    R fun_rhs(const R2 P, const int cc, const R t) {
        R x = P.x, y = P.y;

        return -(y*exp(-t/4)*(t*t - 2*t*x + x*x + y*y - 4))
        /(4*(t*t - 2*t*x + x*x + y*y)*sqrt(t*t - 2*t*x + x*x + y*y));
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
#define cg
// Set numerical example (options: "frachon1", "zahedi1")
#define frachon1
// Set scheme for the dg method (options: "conservative", "classical" see thesis. Irrelevant if "cg" is defined instead of "dg")
#define classical
// Set stabilization method (options: "fullstab", "macro") 
#define macro       
// Decide whether to solve for level set function, or to use exact (options: "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h
#define use_t

#if defined(frachon1)
    using namespace NumericSurfactantEllipse2D;
#elif defined(zahedi1)
    using namespace Zahedi1;
#elif defined(shi1)
    using namespace Shi1;
#endif

// Shi1
// fullstab: [0.000683471, 0.000170117, 4.33821e-05, 1.08361e-05]
// [0.000683471, 0.000170117, 4.33821e-05, 1.08361e-05]
// macro: [0.000674731, 0.000165793, 4.24598e-05, 1.06673e-05]
//[0.000655723, 0.000158013, 3.99581e-05, 1.02663e-05]

// Frachon1 
// fullstab: [0.00245676, 0.000604105, 0.00014715, 3.82835e-05]
// macro: [0.00195667, 0.00046839, 0.000111216, 2.96682e-05]

// fullstab [0.000470812, 0.000513443, 0.000608677, 0.000633784, 0.000639416]
// macro [0.000516745, 0.000389518, 0.000474983, 0.000503638, 0.000512846]

// fullstab: [0.000666664, 0.000406623, 0.000421387
// macro: [0.000164092, 0.000380699, 0.000495177

int main(int argc, char** argv) {
    
    // Mesh settings and data objects
    const size_t iterations = 5;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 20, ny = 20;       // starting mesh size
    double h = 0.1;             // starting mesh size
    double dT = 0.125;

#if defined(frachon1)
    // Paths to store data
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeSurface/frachon1/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeSurface/frachon1/paraview/";  
#elif defined(zahedi1)  
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeSurface/zahedi1/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeSurface/zahedi1/paraview/";  
#elif defined(shi1)
    const std::string pathOutputFolder = "../outputFiles/SpaceTimeSurface/shi1/data/";
    const std::string pathOutputFigures = "../outputFiles/SpaceTimeSurface/shi1/paraview/";  
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
    std::array<double, iterations> errorsSurface;          // array to hold bulk errors
    std::array<int, iterations> number_of_stabilized_edges;          // array to count stabilized edges
    std::array<double, iterations> hs;                  // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {
        
    #if defined(frachon1) or defined(zahedi1)
        const double lx = 4., ly = 4.; 
        #ifdef use_h
        nx = (int)(lx/h)+1, ny = (int)(ly/h)+1;
        #elif defined(use_n)
        h = lx/(nx-1);
        #endif
        Mesh Th(nx, ny, -2, -2, lx, ly);
    #elif defined(shi1)
        const double lx = 8., ly = 6.; 
        #ifdef use_h
        nx = (int)(lx/h)+1, ny = (int)(ly/h)+1;
        #elif defined(use_n)
        h = lx/(nx-1);
        #endif
        Mesh Th(nx, ny, -3, -3, lx, ly);

    #endif


        // [0.00483141, 0.00322363, 0.00183045, 0.00112077, 0.00112941]
        // Number of stabilized edges = [12, 14, 8, 8, 8]

        // Errors Bulk = [0.00714964, 0.00323403, 0.00183045, 0.00112077, 0.00112941]
        // Number of stabilized edges = [21, 14, 8, 8, 8]

        //// Parameters
        //double tfinal = 1.25;            // Final time
        double tfinal = 0.25;            // Final time
        
    #ifdef use_t
        GTime::total_number_iteration = int(tfinal/dT);
    #else
        int divisionMeshSize = 3;
        //int divisionMeshSize = 2*3*pi;
        //int divisionMeshSize = 18;

        double dT = h/divisionMeshSize;
        //double dT = 3*h;
        GTime::total_number_iteration = int(tfinal/dT);
    #endif
        dT = tfinal / GTime::total_number_iteration;
        GTime::time_step = dT;
        
        hs.at(j) = h;
        dts.at(j) = dT;

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

        
        double D = 1.;
  
        
        // Constants for penalty terms
        double tau_a0 = 50;        // diffusion penalty scaling
        double tau_b0 = 1;        // convection penalty scaling
        
        // Bulk penalties
        double lambdaA = tau_a0/h;      // diffusion term
        double lambdaB = tau_b0;            // convection term


        #ifdef dg
        // DG stabilization parameters
        double tau0 = 1e-1, tau1 = 1e-1, tau2 = 1e-1;
        #elif defined(cg)
        // CG stabilization parameters
        double tau0 = 0, tau1 = 1e-1, tau2 = 1e-1;
        #endif        


        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
    #ifdef dg
        FESpace2 Vh(Th, DataFE<Mesh>::P1dc);        // discontinuous basis functions
    #elif defined(cg)
        FESpace2 Vh(Th, DataFE<Mesh>::P1);          // continuous basis functions
    #endif
        FESpace2 Vh2(Th, DataFE<Mesh>::P1);          // higher order space
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
        #ifdef shi1 
        Lagrange2 FEvelocity(2);
        #else
        Lagrange2 FEvelocity(2);
        #endif
        FESpace2 VelVh(Th, FEvelocity);
        vector<Fun_h> vel(nbTime);

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
        
        // Convection-Diffusion Problem Object
        CutFEM<Mesh> surfactant(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;
        
        int iter = 0;
        double q0_0, q0_1, qp_0, qp_1;
        double intF = 0, intG = 0;      // hold integrals of rhs and Neumann bcs
        double errL2 = 0.;

        // Iterate over time-slabs
        while (iter < GTime::total_number_iteration) {
            
            GTime::current_iteration = iter;
            const TimeSlab &In(Ih[iter]);

            double tid = GTime::current_time();

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
                vel[i].init(VelVh, fun_velocity, tt);
            #endif
                
                // FIXME: This version of the code has a memory leak (it doesn't delete the interface pointers)
                interface.init(i, Th, ls[i]);
                
            #if defined(levelsetsolve)
                // We solve for the level-set using Crank-Nicholson in time
                if (i<lastQuadTime) {
                    LevelSet::move(ls[i], vel, vel, dt_levelSet, ls[i+1]);
                }
            #endif
            }
            
            
            // Create active meshes
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);
            CutSpace Wh(ThGamma, Vh);
            
            // Data for initial solution
            Rn datau0;
            surfactant.initSpace(Wh, In);
            surfactant.initialSolution(datau0);
            KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
            if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
            Rn uh(datau0);
            Fun_h u0(Wh, datau0);


            // Objects needed for the weak form
            Normal n;
            Conormal t;

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega2                 
            
        //// Assembling linear and bilinear forms
        
        // Partial derivative in time, bulk
            
        // reynold scheme
        #ifdef conservative  

            surfactant.addBilinear(
                - innerProduct(u, dt(v))
                , interface
                , In
            );

            surfactant.addBilinear(
                + innerProduct(u, v)
                , *interface(lastQuadTime)
                , In
                , (int)lastQuadTime
            );

            // Time penalty term bulk RHS
            surfactant.addLinear(
                + innerProduct(u0.expression(), v)
                , *interface(0)
                , In
                , 0
            );

            surfactant.addBilinear(
                innerProduct(u, v)
                , *interface(0)
                , In
                , 0
            );
            surfactant.addLinear(
                innerProduct(u0.expression(), v)
                , *interface(0)
                , In
                , 0
            );

        // classical scheme
        #elif defined(classical)   

            surfactant.addBilinear(
                + innerProduct(dt(u), v)
                , interface
                , In
            );

            surfactant.addBilinear(
                innerProduct(u, v)
                , *interface(0)
                , In
                , 0
            );
            surfactant.addLinear(
                innerProduct(u0.expression(), v)
                , *interface(0)
                , In
                , 0
            );

        #endif
            
        //// Scheme for diffusion

            surfactant.addBilinear(
                + innerProduct(D*gradS(u), gradS(v))
                , interface
                , In
            );
            
        #ifdef dg
            // integral on inner edges (E_h)
            surfactant.addBilinear(
                - innerProduct(D*average(gradS(u)*t, 0.5, -0.5), jump(v))      
                - innerProduct(D*jump(u), average(gradS(v)*t, 0.5, -0.5))     
                + innerProduct(lambdaA*jump(u), jump(v))          
                , interface
                , INTEGRAL_INNER_NODE_2D
                , In
            );
        #endif
            
        //// Schemes for convection

        #if defined(dg) && (defined(conservative) || defined(classical))

            // Convection terms
                       
            for (int i=0; i<nbTime; ++i) {
                
                // Convection term
                surfactant.addBilinear(
                    - innerProduct(u, (vel[i].expression()*grad(v)))
                    //- innerProduct(u*dxS(vel[i]) + u*dyS(vel[i]), v)
                    , interface
                    , In
                    , i
                );

                // Flux terms 
                surfactant.addBilinear(
                    + innerProduct(average(vel[i]*t*u, 0.5, -0.5), jump(v))
                    + innerProduct(lambdaB*jump(fabs(vel[i]*t)*u), jump(v))
                    // + innerProduct(average((vel.at(i) * t) * u), jump(v))
                    // - innerProduct(jump(u), average((vel.at(i) * t) * v))
                    // + innerProduct(1 * jump(u * fabs(vel.at(i) * n)), jump(v))
                    , interface
                    , INTEGRAL_INNER_NODE_2D
                    , In
                    , i
                );
            }
        
        #elif defined(cg)

            #if defined(classical)
            for(int i=0;i<nbTime;++i) {      // computation of the curvature
                ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
                ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
                surfactant.addBilinear(
                    + innerProduct(dx(u)*vx + dy(u)*vy, v)
                    + innerProduct(u*dxS(vel[i]) + u*dyS(vel[i]), v)
                    , interface
                    , In
                    , i
                );
            }

            #elif defined(conservative)
            for(int i=0;i<nbTime;++i) {      // computation of the curvature
                ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
                ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
                surfactant.addBilinear(
                    - innerProduct(u, dx(v)*vx + dy(v)*vy)
                    , interface
                    , In
                    , i
                );
            }
            #endif

        #endif

        // Stabilization

        #ifdef macro    
            TimeMacroElementSurface<Mesh> TimeMacro(ThGamma, interface, qTime, 0.125);
            number_of_stabilized_edges.at(j) = TimeMacro.number_of_inner_edges();

            if (iterations == 1 && h > 0.01) {
                Paraview<Mesh> writerMacro(Th, pathOutputFigures + "Th" + to_string(iter+1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);  
                std::cout << writerMacro.get_nb_stab_elems(ThGamma, 0) << std::endl;
                // domain = 0, 

                writerMacro.writeFaceStab(ThGamma, 0, pathOutputFigures + "FullStabilization" + to_string(iter+1) + ".vtk");
                writerMacro.writeActiveMesh(ThGamma, pathOutputFigures + "ActiveMesh" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0, pathOutputFigures + "macro" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroInnerEdge(TimeMacro, 0, pathOutputFigures + "macro_inner_edge" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroOutterEdge(TimeMacro, 0, pathOutputFigures + "macro_outer_edge" + to_string(iter+1) + ".vtk");
                writerMacro.writeSmallElements(TimeMacro, 0, pathOutputFigures + "small_element" + to_string(iter+1) + ".vtk");
            }

            surfactant.addFaceStabilization(
                + innerProduct(tau0/h/h*jump(u), jump(v))
                + innerProduct(tau1*jump(grad(u)), jump(grad(v)))
                , ThGamma
                , In
                , TimeMacro
            );

        #elif defined(fullstab)

            surfactant.addFaceStabilization(
                + innerProduct(tau0/h/h*jump(u), jump(v))
                + innerProduct(tau1*jump(grad(u)), jump(grad(v)))
                , ThGamma
                , In
            );
        #endif
            
            surfactant.addBilinear(
                + innerProduct(tau2*h*h*grad(u)*n, grad(v)*n)
                , interface
                , In
            );
            
            
            Fun_h funrhs(Vh2, In, fun_rhs);

            // Add RHS on surface
            surfactant.addLinear(
                + innerProduct(funrhs.expression(), v)
                , interface
                , In
            );

            // Compute integrals
            intF = integral(funrhs, In, interface, 0);
            
    #ifdef conservative
        #ifdef fullstab
            #ifdef use_t
            if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_conservative_fullstab_h" + to_string(h) + "_" + to_string(j+1) + ".dat");
            #else
            if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_conservative_fullstab_j" + to_string(j+1) + ".dat");
           #endif
        #elif defined(macro)
           #ifdef use_t
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_conservative_macro_h" + to_string(h) + "_" + to_string(j+1) + ".dat");
           #else
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_conservative_macro_j" + to_string(j+1) + ".dat");
           #endif
        #endif
    #elif defined(classical)
        #ifdef fullstab
           #ifdef use_t
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_classical_fullstab_h" + to_string(h) + "_" + to_string(j+1) + ".dat");
           #else
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_classical_fullstab_j" + to_string(j+1) + ".dat");
           #endif
        #elif defined(macro)
           #ifdef use_t
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_classical_macro_h" + to_string(h) + "_" + to_string(j+1) + ".dat");
           #else
           if (iter == GTime::total_number_iteration-1) matlab::Export(surfactant.mat_, pathOutputFolder + "mat_classical_macro_j" + to_string(j+1) + ".dat");
           #endif
        #endif
    #endif

            // Solve linear system
            surfactant.solve();

            KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
            uh = dw;
            surfactant.saveSolution(uh);

            // Compute conservation error

            // Compute error
            {
                Rn sol(Wh.get_nb_dof(), 0.);
                sol += uh(SubArray(Wh.get_nb_dof(), 0));
                Fun_h funuh(Wh, sol);
                errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(0), tid,0,1);
                std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

                sol  += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
                errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(nbTime-1), tid+dT,0,1);
                std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
            }   
            
            errorsSurface.at(j) = errL2;

            if (iterations == 1) {
                Fun_h sol(Wh, uh);
            
                Paraview<Mesh> writer(ThGamma, pathOutputFigures + "surfactant_"+to_string(iter+1)+".vtk");
            
                Fun_h uS_ex(Wh, fun_sol_surfactant, GTime::current_time());
                writer.add(sol , "surfactant", 0, 1);
                writer.add(uS_ex , "surfactant_exact", 0, 1);
                writer.add(ls[0], "levelSet",  0, 1);
                // writer.add(ls[2], "levelSet2", 0, 1);
            
            }

            if (iterations > 1 && iter == GTime::total_number_iteration-1) outputData << h << "," << dT << "," << errL2 << std::endl;

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
    std::cout << "Errors Surface = [";
    for (int i=0; i<iterations; i++) {

        std::cout << errorsSurface.at(i);
        if (i < iterations-1) {
            std::cout << ", ";
        }

    }
    std::cout << "]" << std::endl;

    std::cout << std::endl;
    std::cout << "Number of stabilized edges = [";
    for (int i=0; i<iterations; i++) {

        std::cout << number_of_stabilized_edges.at(i);
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

