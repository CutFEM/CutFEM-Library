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
        return exp(-t/4)*sin(asin(y/sqrt((x-t)*(x-t) + y*y))) + 2;
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
#define shi1
// Set scheme for the dg method (options: "conservative", "classical" see thesis. Irrelevant if "cg" is defined instead of "dg")
#define conservative
// Set stabilization method (options: "fullstab", "macro") 
#define fullstab       
// Decide whether to solve for level set function, or to use exact (options: "levelsetsolve", "levelsetexact")
#define levelsetexact

#define use_h

#if defined(frachon1)
    using namespace NumericSurfactantEllipse2D;
#elif defined(zahedi1)
    using namespace Zahedi1;
#elif defined(shi1)
    using namespace Shi1;
#endif

int main(int argc, char** argv) {
    
    // Mesh settings and data objects
    const size_t iterations = 4;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 20, ny = 20;       // starting mesh size
    double h = 0.1;             // starting mesh size

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

        // Mesh size
        int divisionMeshSize = 3;
        double dT = h/divisionMeshSize;
        hs.at(j) = h;
        dts.at(j) = dT;
        
        double tfinal = 0.25;            // Final time
        GTime::total_number_iteration = (int)(tfinal/dT);
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
        double tau0 = 0, tau1 = 1e-2, tau2 = 1e-2;
        #endif        


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
        #ifdef shi1 
        Lagrange2 FEvelocity(0);
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
            TimeMacroElement<Mesh> TimeMacro(Kh2, qTime, 0.125);

            surfactant.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), jump(v))
                + innerProduct(h*tau21*jump(grad(u)), jump(grad(v)))
                , Kh2
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
            
            
            Fun_h funrhs(Wh, In, fun_rhs);

            // Add RHS on surface
            surfactant.addLinear(
                + innerProduct(funrhs.expression(), v)
                , interface
                , In
            );

            // Compute integrals
            intF = integral(funrhs, In, interface, 0);
            
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


// #include <cassert>
// #include <cmath>
// #include <cstdlib>
// #include <fstream>
// #include <iostream>
// #include <experimental/filesystem>

// #ifdef USE_MPI
// #  include "cfmpi.hpp"
// #endif
// #include "baseProblem.hpp"
// #include "../num/matlab.hpp"
// #include "paraview.hpp"

// namespace NumericSurfactant2D0 {
//   R fun_levelSet(const R2 P, const int i, const double t ) {
//     R x = P.x,  y = P.y;
//     return sqrt((x-t)*(x-t) + y*y) - 1;
//     // return sqrt((x)*(x) + y*y) - 1;
//   }
//   R fun_init_surfactant(const R2 P, int i) {    // = 8PI
//     return P.y/1+2;
//   }
//   R fun_velocity(const R2 P, int i,const double t) {
//     return (i==0)?1:0;
//   }
//   // R fun_sol_surfactant(const R2 P,  const int i, const R t) {
//   //   return P.y/1+2;;
//   // }
//   R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-0.25*t);}

//   R fun_rhs(const R2 P, const int cc, const R t) {
//     R x = P.x,  y = P.y;
//     R r = y*exp(-t/4) - (x*y*exp(-t/4))/4 - (y*exp(-t/4)*(3*t - 4*x))/(t*t - 2*t*x + x*x + y*y);
//     return r;
//   }
//   // R fun0(const R2 P, const R t) {return P.x*P.y*exp(-4*t);}
//   // R fun_sol_surfactant(const R2 P, const R t) {
//   //   R x = P.x,  y = P.y;
//   //   return exp(-t/4)*y/(sqrt((x-t)*(x-t)+y*y)) + 2;
//   // }
//   // R2 fun_velocity_field(const R2 P, int i) {
//   //   //return R2((P.y+2)*(P.y+2)/3,0);
//   //   return (i==0);
//   // }
//   // R fun_div_surfactant(const R2 P, const R t) {
//   //   return 0;
//   // }
//   // R fun_rhs_surfactant(const R2 P, const R t) {
//   //   return 0;
//   // }

//   // R fun_levelSet(const R2 P) {
//   //   R x = P.x,  y = P.y;
//   //   return sqrt(x*x + y*y) - 2;
//   // }
//   // Diff<R,2> fun_levelSet_timeDA(Diff<R,2> P[2], const R t) {
//   //   return sqrt((P[0]-t)*(P[0]-t) + P[1]*P[1]) - 2;
//   // }

// }
// namespace NumericSurfactantEllipse2D {
//   R fun_init_surfactant(const R2 P, const int i) {return P.x*P.y ;} //2*... = 4*pi
//   R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-t/4);}
//   R fun_velocity(const R2 P, const int i, const R t) {
//     R a = (1+0.25*sin(2*M_PI*t));
//     R b = M_PI/4*(cos(2*M_PI*t))/a;
//     return (i==0) ? b*P.x : 0;
//   }
//   R fun_levelSet(const R2 P, const int i, const R t) {
//     R a = (1+0.25*sin(2*M_PI*t));
//     // a = 1;
//     return sqrt(P.x*P.x/a + P.y*P.y) - 1.0000079;
//   }
//   R fun_rhs(const R2 P, const int cc, const R t) {
//     R x = P.x,  y = P.y;
//     // R r = (4*x*y*exp(-4*t)*(sin(2*M_PI*t) + 4)*
//     //        (2048*x*x*y*y + 336*y*y*y*y*pow(sin(2*M_PI*t),2) +
//     //         52*y*y*y*y*pow(sin(2*M_PI*t),3) + 3*y*y*y*y*pow(sin(2*M_PI*t),4)
//     //         + 64*x*x*x*x*sin(2*M_PI*t) + 960*y*y*y*y*sin(2*M_PI*t) + 1024*x*x*x*x
//     //         + 1024*y*y*y*y + 144*x*x*y*y*pow(sin(2*M_PI*t),2)
//     //         + 4*x*x*y*y*pow(sin(2*M_PI*t),3) + 1024*x*x*y*y*sin(2*M_PI*t)))
//     //   /pow((y*y*pow(sin(2*t*M_PI),2) + 8*y*y*sin(2*t*M_PI) + 16*x*x + 16*y*y),3)
//     //   - 4*x*y*exp(-4*t) + (x*y*M_PI*exp(-4*t)*cos(2*M_PI*t))/(sin(2*M_PI*t) + 4)
//     //   + (x*y*y*y*M_PI*exp(-4*t)*cos(2*M_PI*t)*(sin(2*M_PI*t) + 4))/
//     //   (y*y*pow(sin(2*M_PI*t),2) + 8*y*y*sin(2*M_PI*t) + 16*x*x + 16*y*y);
//     R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t)
//     + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi)
//     + 16*x*x + 16*y*y, 2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t))/(sin(2*pi*t) + 4)
//     + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
//     return r;
//   }
//   // R fun_rhs(const R2 P, const int cc, const R t) {
//   //   R x = P.x,  y = P.y;
//   //   R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t) + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi) + 16*x*x + 16*y*y,2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t)*(y*y*sin(2*pi*t) - 4*x*x + 4*y*y))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y) + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
//   //   return r;
//   // }
// }

// #define FORMULATION1;
// // #define FORMULATION2;

// using namespace NumericSurfactantEllipse2D ;
// // using namespace NumericSurfactant2D0 ;
// typedef Mesh2 Mesh;
// typedef FESpace2 Space;
// typedef CutFESpace<Mesh> CutSpace;
// typedef TestFunction<2> FunTest;
// typedef FunFEM<Mesh2> Fun_h;

// double solveCG_1(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
//   const int d = 2;

//   const double cpubegin = CPUtime();

//   double errL2 = 0.;

//   // BACKGROUND MESH
//   // Mesh Kh(nx, nx, -2, -2, 4., 4.);
//   Space Vh(Kh, DataFE<Mesh2>::P1);
//   double meshSize = hi;//(4./(nx-1));

//   // TIME PARAMETER
//   int divisionMeshsize = 4;
//   // double dT = dt;//meshSize/divisionMeshsize;
//   GTime::total_number_iteration = int(tfinal/dT);
//   dT = tfinal / GTime::total_number_iteration;
//   GTime::time_step = dT;
//   Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
//   FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
//   const QuadratureFormular1d& qTime(*Lobatto(3));
//   const Uint nbTime = qTime.n;
//   const Uint ndfTime = Ih[0].NbDoF();
//   const Uint lastQuadTime = nbTime-1;


//   //CREATE SPACE AND VELOCITY FUNCTION
//   Lagrange2 FEvelocity(2);
//   Space VelVh  (Kh, FEvelocity);
//   vector<Fun_h> vel(nbTime);


//   // CREATE THE LEVELSET
//   Space Lh  (Kh, DataFE<Mesh2>::P1);
//   double dt_levelSet = dT/(nbTime-1);
//   vector<Fun_h> ls(nbTime);
//   for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
//   // LevelSet2 levelSet(Lh_k);

//   // TIME INTERFACE AT ALL QUADRATURE TIME
//   TimeInterface<Mesh> interface(qTime);

//   // INITIALIZE THE PROBLEM
//   CutFEM<Mesh2> surfactant(qTime);
//   const R epsilon_S = 1.;

//   int iter = 0, iterfig = 0;
//   while( iter < GTime::total_number_iteration ) {

//     GTime::current_iteration = iter;
//     const TimeSlab& In(Ih[iter]);

//     std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
//     double tid = GTime::current_time();
//     // INITIALIZE LEVELSET, VELOCITY, INTERFACE
//     ls.begin()->swap(ls[nbTime-1]);
//     for(int i=0;i<nbTime;++i) {
//       R tt = In.Pt(R1(qTime(i).x));
//       ls[i].init(Lh, fun_levelSet, tt);
//       vel[i].init(VelVh, fun_velocity, tt);
//       interface.init(i,Kh,ls[i]);
//     }

//     // CREATE THE ACTIVE MESH AND CUT SPACE
//     ActiveMesh<Mesh> Kh0(Kh);
//     Kh0.createSurfaceMesh(interface);
//     CutSpace Wh(Kh0, Vh);

//     // INITIALIZE THE PROBLEM
//     Rn datau0;
//     surfactant.initSpace(Wh, In);
//     surfactant.initialSolution(datau0);
//     KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
//     if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
//     Rn uh(datau0);
//     Fun_h u0(Wh, datau0);

//     Normal n;
//     FunTest s(Wh,1), r(Wh,1);


//     std::cout << " Begin assembly " << std::endl;
//     double tt0 = MPIcf::Wtime();

//     surfactant.addBilinear(
//       innerProduct(dt(s), r)
//       + innerProduct(epsilon_S*gradS(s), gradS(r))
//       , interface
//       , In
//     );

//     for(int i=0;i<nbTime;++i) {      // computation of the curvature
//       ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
//       ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
//       surfactant.addBilinear(
//         innerProduct(dx(s)*vx + dy(s)*vy, r)
//         + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
//         , interface
//         ,In, i
//       );
//     }
//     surfactant.addFaceStabilization(
//       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
//       , Kh0 
//       , In
//     );

//     //     // surfactant.addBilinearFormInterface(
//     //     //   innerProduct(1e-2*grad(s)*n, grad(r)*n)
//     //     //   ,In
//     //     // );
//     //
//         surfactant.addBilinear(
//           innerProduct(s, r)
//           , *interface(0)
//           , In
//           , 0
//         );
//         surfactant.addLinear(
//           innerProduct(u0.expression(), r)
//           , *interface(0)
//           , In
//           , 0
//         );

//         Fun_h funrhs(Wh, In, fun_rhs);
//         // Fun_h funrhsp(cutVh, In);
//         // projection(funrhs, funrhsp, In, interface ,0);
//         surfactant.addLinear(
//           innerProduct(funrhs.expression(), r)
//           , interface
//           , In
//         );
//         // double  intF = integralSurf(funrhsp, In, interface) ;
//     //     intFnp = integralSurf(funrhs, In, qTime);
//     //
//     //
//     // matlab::Export(surfactant.mat_, "matA.dat");
//     // surfactant.cleanMatrix();
//     surfactant.solve();

//     KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
//     uh = dw;
//     surfactant.saveSolution(uh);
//     //



//     //  -----------------------------------------------------
//     //                COMPUTE ERROR
//     //  -----------------------------------------------------
//     {
//       Rn sol(Wh.get_nb_dof(), 0.);
//       sol += uh(SubArray(Wh.get_nb_dof(), 0));
//       Fun_h funuh(Wh, sol);
//       errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(0), tid,0,1);
//       std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

//       sol  += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
//       errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(nbTime-1), tid+dT,0,1);
//       std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
//     }



//     //  -----------------------------------------------------
//     //                     PLOTTING
//     //  -----------------------------------------------------
//     // if(MPIcf::IamMaster() ) {
//     //     // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
//     //     Fun_h sol(Wh, uh);
//     //
//     //     // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
//     //     Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");
//     //
//     //     writer.add(sol , "surfactant", 0, 1);
//     //     writer.add(ls[0], "levelSet",  0, 1);
//     //     // writer.add(ls[2], "levelSet2", 0, 1);
//     //   }


//     // return 0.;
//     iter += 1;
//   }

//   return errL2;
// }
// double solveCG_2(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
//   const int d = 2;
//   // typedef Mesh2 Mesh;
//   // typedef FESpace2 Space;
//   // typedef CutFESpace<Mesh> CutSpace;
//   // typedef TestFunction<2> FunTest;
//   // typedef FunFEM<Mesh2> Fun_h;

//   const double cpubegin = CPUtime();

//   double errL2 = 0.;

//   // BACKGROUND MESH
//   // Mesh Kh(nx, nx, -2, -2, 4., 4.);
//   Space Vh(Kh, DataFE<Mesh2>::P1);
//   double meshSize = hi;//(4./(nx-1));

//   // TIME PARAMETER
//   int divisionMeshsize = 4;
//   // double dT = 0.01;//meshSize/divisionMeshsize;
//   GTime::total_number_iteration = int(tfinal/dT);
//   dT = tfinal / GTime::total_number_iteration;
//   GTime::time_step = dT;
//   Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
//   FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
//   const QuadratureFormular1d& qTime(*Lobatto(3));
//   const Uint nbTime = qTime.n;
//   const Uint ndfTime = Ih[0].NbDoF();
//   const Uint lastQuadTime = nbTime-1;



//   //  -----------------------------------------------------
//   //             Create files to save results
//   //  -----------------------------------------------------
//   // std::string pathOutpuFolder = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
//   // std::string pathOutpuFigure = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
//   // CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
//   // cout << " path to the output files : "+pathOutpuFolder << std::endl;
//   // cout << " path to the vtk files : "+pathOutpuFigure << std::endl;
//   std::ofstream outputData("dataDrop"+to_string(i)+".dat", std::ofstream::out);
//   // cout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
//   // cout << " Creating the file \"output.txt\" to save cout" << std::endl;


//   //CREATE SPACE AND VELOCITY FUNCTION
//   Lagrange2 FEvelocity(2);
//   Space VelVh  (Kh, FEvelocity);
//   vector<Fun_h> vel(nbTime);


//   // CREATE THE LEVELSET
//   Space Lh  (Kh, DataFE<Mesh2>::P1);
//   double dt_levelSet = dT/(nbTime-1);
//   vector<Fun_h> ls(nbTime);
//   for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
//   // LevelSet2 levelSet(Lh_k);

//   // TIME INTERFACE AT ALL QUADRATURE TIME
//   TimeInterface<Mesh> interface(qTime);

//   // INITIALIZE THE PROBLEM
//   CutFEM<Mesh2> surfactant(qTime);
//   const R epsilon_S = 1.;

//   int iter = 0, iterfig = 0;
//   double q_init0, q_init1, qp1;
//   while( iter < GTime::total_number_iteration ) {

//     GTime::current_iteration = iter;
//     const TimeSlab& In(Ih[iter]);

//     std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
//     double tid = GTime::current_time();
//     // INITIALIZE LEVELSET, VELOCITY, INTERFACE
//     ls.begin()->swap(ls[nbTime-1]);
//     for(int i=0;i<nbTime;++i) {
//       R tt = In.Pt(R1(qTime(i).x));
//       ls[i].init(Lh, fun_levelSet, tt);
//       vel[i].init(VelVh, fun_velocity, tt);
//       interface.init(i,Kh,ls[i]);
//     }

//     // CREATE THE ACTIVE MESH AND CUT SPACE
//     ActiveMesh<Mesh> Kh0(Kh);
//     Kh0.createSurfaceMesh(interface);
//     CutSpace Wh(Kh0, Vh);

//     // INITIALIZE THE PROBLEM
//     Rn datau0;
//     surfactant.initSpace(Wh, In);
//     surfactant.initialSolution(datau0);
//     KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
//     if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
//     Rn uh(datau0);
//     Fun_h u0(Wh, datau0);

//     Normal n;
//     FunTest s(Wh,1), r(Wh,1);


//     std::cout << " Begin assembly " << std::endl;
//     double tt0 = MPIcf::Wtime();

//     surfactant.addBilinear(
//       -innerProduct(s, dt(r))
//       + innerProduct(epsilon_S*gradS(s), gradS(r))
//       , interface
//       , In
//     );

//     for(int i=0;i<nbTime;++i) {
//       ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
//       ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
//       surfactant.addBilinear(
//         - innerProduct(s, dx(r)*vx)
//         - innerProduct(s, dy(r)*vy)
//         , interface
//         , In
//         , i
//       );
//     }
//     // surfactant.addEdgeIntegral(
//     //   innerProduct(1e-2*jump(grad(s).t()*n), jump(grad(r).t()*n)),
//     //   In
//     // );


//     Fun_h funrhs (Wh, In, fun_rhs);
//     // Fun_h funrhsp(cutVh, In);
//     // projection(funrhs, funrhsp, In, interface ,0);
//     surfactant.addLinear (
//       innerProduct(funrhs.expression(), r)
//       , interface
//       , In
//     );
//     double intF   = integral(funrhs, In, interface, 0) ;
//     // intFnp = integralSurf(funrhs , In, qTime);

//     surfactant.addBilinear(
//       innerProduct(s, r)
//       , *interface(lastQuadTime)
//       , In
//       , lastQuadTime
//     );
//     surfactant.addLinear(
//       innerProduct(u0.expression(), r)
//       , *interface(0)
//       , In
//       , 0
//     );

//     // FunTest D2un = grad(grad(r)*n)*n;
//     surfactant.addBilinear(
//       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
//       , Kh0
//       , INTEGRAL_INNER_EDGE_2D
//       , In
//     );
//     // if stab interface, can have h in both stabilization
//     // surfactant.addBilinear(
//     //   innerProduct(1e-1*grad(s)*n, grad(r)*n)
//     //   , interface
//     //   , In
//     // );


//     // matlab::Export(surfactant.mat_, "matA.dat");
//     // surfactant.cleanMatrix();
//     surfactant.solve();

//     KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
//     uh = dw;
//     surfactant.saveSolution(uh);



//     Rn sol(Wh.get_nb_dof(), 0.);
//     sol += uh(SubArray(Wh.get_nb_dof(), 0));
//     sol += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
//     Fun_h funuh_0(Wh, uh);
//     Fun_h funuh_1(Wh, sol);

//     //  -----------------------------------------------------
//     //                COMPUTE ERROR
//     //  -----------------------------------------------------
//     {
//       Fun_h funuh(Wh, sol);
//       errL2 =  L2normSurf(funuh_0, fun_sol_surfactant, interface(0), tid,0,1);
//       std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

//       errL2 =  L2normSurf(funuh_1, fun_sol_surfactant, interface(nbTime-1), tid+dT,0,1);
//       std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
//     }
//     //  -----------------------------------------------------
//     // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
//     //  -----------------------------------------------------
//     {
//       double q0 = integral(funuh_0, interface(0), 0);
//       double q1 = integral(funuh_1, interface(lastQuadTime), 0);
//       if(iter==0){
//         q_init0 = q0;
//         q_init1 = q1;
//         qp1     = q1;
//         q_init1 = integral(u0, interface(0), 0);
//       }

//       outputData << setprecision(10);
//       outputData << std::setw(10) << std::setfill(' ') << GTime::current_time()
//                  << std::setw(20) << std::setfill(' ') << q0
//                  << std::setw(20) << std::setfill(' ') << q1
//                  << std::setw(20) << std::setfill(' ') << intF
//                  << std::setw(20) << std::setfill(' ') << fabs(q1 - qp1 - intF)
//                  << std::endl;
//                  qp1 = q1;
//       }



//     //  -----------------------------------------------------
//     //                     PLOTTING
//     //  -----------------------------------------------------
//     // if(MPIcf::IamMaster() ) {
//     //     // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
//     //     Fun_h sol(Wh, uh);
//     //
//     //     // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
//     //     Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");
//     //
//     //     writer.add(sol , "surfactant", 0, 1);
//     //     writer.add(ls[0], "levelSet",  0, 1);
//     //     // writer.add(ls[2], "levelSet2", 0, 1);
//     //   }


//     // return 0.;
//     iter += 1;
//   }


//   outputData.close();
//   return errL2;
// }

// double solveDG_1(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
//   const int d = 2;
//   // typedef Mesh2 Mesh;
//   // typedef FESpace2 Space;
//   // typedef CutFESpace<Mesh> CutSpace;
//   // typedef TestFunction<2> FunTest;
//   // typedef FunFEM<Mesh2> Fun_h;

//   const double cpubegin = CPUtime();

//   double errL2 = 0.;

//   // BACKGROUND MESH
//   // Mesh Kh(2*nx, nx, -2, -2, 8., 4.);
//   Space Vh(Kh, DataFE<Mesh2>::P1);
//   double meshSize = hi;//(4./(nx-1));
//   // double hi = meshSize;

//   // TIME PARAMETER
//   int divisionMeshsize = 4;
//   // double dT = hi/8;//meshSize/divisionMeshsize;
//   GTime::total_number_iteration = int(tfinal/dT);
//   dT = tfinal / GTime::total_number_iteration;
//   GTime::time_step = dT;
//   Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
//   FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
//   const QuadratureFormular1d& qTime(*Lobatto(3));
//   const Uint nbTime = qTime.n;
//   const Uint ndfTime = Ih[0].NbDoF();
//   const Uint lastQuadTime = nbTime-1;


//   //CREATE SPACE AND VELOCITY FUNCTION
//   Lagrange2 FEvelocity(2);
//   Space VelVh  (Kh, FEvelocity);
//   vector<Fun_h> vel(nbTime);

//   // CREATE THE LEVELSET
//   Space Lh  (Kh, DataFE<Mesh2>::P1);
//   double dt_levelSet = dT/(nbTime-1);
//   vector<Fun_h> ls(nbTime);
//   for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
//   // LevelSet2 levelSet(Lh_k);

//   // TIME INTERFACE AT ALL QUADRATURE TIME
//   TimeInterface<Mesh> interface(qTime);

//   // INITIALIZE THE PROBLEM
//   CutFEM<Mesh2> surfactant(qTime);
//   const R epsilon_S = 1.;

//   int iter = 0, iterfig = 0;
//   while( iter < GTime::total_number_iteration ) {

//     GTime::current_iteration = iter;
//     const TimeSlab& In(Ih[iter]);

//     std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
//     double tid = GTime::current_time();
//     // INITIALIZE LEVELSET, VELOCITY, INTERFACE
//     ls.begin()->swap(ls[nbTime-1]);
//     for(int i=0;i<nbTime;++i) {
//       R tt = In.Pt(R1(qTime(i).x));
//       ls[i].init(Lh, fun_levelSet, tt);
//       vel[i].init(VelVh, fun_velocity, tt);
//       interface.init(i,Kh,ls[i]);
//     }

//     // CREATE THE ACTIVE MESH AND CUT SPACE
//     ActiveMesh<Mesh> Kh0(Kh);
//     Kh0.createSurfaceMesh(interface);
//     CutSpace Wh(Kh0, Vh);

//     // INITIALIZE THE PROBLEM
//     Rn datau0;
//     surfactant.initSpace(Wh, In);
//     surfactant.initialSolution(datau0);
//     KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
//     if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
//     Rn uh(datau0);
//     Fun_h u0(Wh, datau0);

//     Normal n;
//     Conormal conormal;
//     FunTest s(Wh,1), r(Wh,1);


//     std::cout << " Begin assembly " << std::endl;
//     double tt0 = MPIcf::Wtime();

//     surfactant.addBilinear(
//       + innerProduct(dt(s), r)
//       + innerProduct(epsilon_S*gradS(s), gradS(r))
//       , interface
//       , In
//     );


//     for(int i=0;i<nbTime;++i) {      // computation of the curvature
//         ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
//         ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
//         // Interior edges surface convective term
//         // surfactant.addBilinear(
//         //   - innerProduct(s, dx(r)*vx + dy(r)*vy)
//         //   //+ innerProduct(dx(s)*vx + dy(s)*vy, r)  // TEST
//         //   + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
//         //   , interface
//         //   , In , i
//         // );
//         // surfactant.addBilinear(
//         //   + innerProduct(average((vel[i] * conormal)*s, 0.5, -0.5), jump(r))
//         //   //+ innerProduct(average((vel.at(i) * t)*s, 1, 1), average(r)) // TEST
//         //   //- innerProduct(jump(s), average((vel.at(i) * t) * r, 0.5, -0.5))
//         //   //+ innerProduct(10 * jump((vel.at(i)*t)*s), jump(r))
//         //   + innerProduct(0*average(fabs(vel[i]*conormal)*s, 1, 1), jump(r))
//         //   , interface
//         //   , innerRidge
//         //   , In , i
//         // );
//         surfactant.addBilinear(
//             + innerProduct(dx(s)*vx + dy(s)*vy, r)*0.5
//             - innerProduct(s, dx(r)*vx + dy(r)*vy)*0.5
//             + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
//             , interface
//             , In , i
//         );

//         // "Variant 2"
//         surfactant.addBilinear(
//             + innerProduct(average((vel[i]*conormal)*s, 0.5, -0.5), jump(r))*0.5
//             - innerProduct(jump(s), average((vel[i]*conormal)*r, 0.5, -0.5))*0.5
//             + innerProduct(10*jump(s), jump(r))
//             //+ innerProduct(lambdaBSdiv*average((vel[i]*t)*s, 1, 1), average((vel[i]*t)*r, 1, 1))                  // (3.16)
//             , interface
//             , INTEGRAL_INNER_NODE_2D
//             , In , i
//         );

//     }

//     // Stabilization
//     surfactant.addFaceStabilization(
//         + innerProduct(1./hi/hi*jump(s), jump(s))
//         + innerProduct(1.*jump(grad(s)*n), jump(grad(r)*n))
//         , Kh0
//         , In
//     );
    
//     surfactant.addBilinear(
//       innerProduct(1e-2*hi*hi*grad(s)*n, grad(r)*n)
//       , interface
//       , In
//     );
      
//     // Time penalty
//     surfactant.addBilinear(
//         innerProduct(s, r)
//         , *interface(0)
//         , In
//         , 0
//     );
//     surfactant.addLinear(
//         innerProduct(u0.expression(), r)
//         , *interface(0)
//         , In
//         , 0
//     );

//     Fun_h funrhs(Wh, In, fun_rhs);
//     // Fun_h funrhsp(cutVh, In);
//     // projection(funrhs, funrhsp, In, interface ,0);
//     surfactant.addLinear(
//         innerProduct(funrhs.expression(), r)
//         , interface
//         , In
//     );
//         // double  intF = integralSurf(funrhsp, In, interface) ;
//     //     intFnp = integralSurf(funrhs, In, qTime);
//     //
//     //
//     // matlab::Export(surfactant.mat_, "matA.dat");
//     // surfactant.cleanMatrix();
//     surfactant.solve();

//     KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
//     uh = dw;
//     surfactant.saveSolution(uh);
//     //



//     //  -----------------------------------------------------
//     //                COMPUTE ERROR
//     //  -----------------------------------------------------
//     {
//       Rn sol(Wh.get_nb_dof(), 0.);
//       sol += uh(SubArray(Wh.get_nb_dof(), 0));
//       Fun_h funuh(Wh, sol);
//       errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(0), tid,0,1);
//       std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

//       sol  += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
//       errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(nbTime-1), tid+dT,0,1);
//       std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
//     }



//     //  -----------------------------------------------------
//     //                     PLOTTING
//     //  -----------------------------------------------------
//     if(MPIcf::IamMaster() ) {
//       // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
//       Fun_h sol(Wh, uh);
//       Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");

//       writer.add(sol , "surfactant", 0, 1);
//       writer.add(ls[0], "levelSet",  0, 1);
//       // writer.add(ls[2], "levelSet2", 0, 1);
//     }


//     // return 0.;
//     iter += 1;
//   }

//   return errL2;
// }


// int main(int argc, char** argv ) {


//   MPIcf cfMPI(argc,argv);
//   int nx = 20;
//   for(int iter=0; iter<4;++iter) {

//     Mesh Kh(nx, nx, -2, -2, 4., 4.);
//     double h = 4./(nx-1);
//     double dt = h/4;
//     double err = solveDG_1(iter, Kh, nx, 0.25, h, dt);

//     nx = 2*nx-1;
//   }
//   return 0;
// }


//   //
//   //
//   //
//   // cout << " nx        \t" << nx << std::endl;
//   // cout << " Mesh size \t" << meshSize << std::endl;
//   // cout << " Time Step \t" << dT << std::endl;
//   // cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
//   // cout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
//   // cout << " number of quadrature points in time : \t" << nbTime << std::endl;
//   // cout << " number of dof in time per time slab : \t" << ndfTime << std::endl;
//   //
//   //
//   // // Set parameters for paraview PLOTTING
//   // const bool writeVTKFiles = true;
//   // const bool saveSurfactantVTK = true;
//   // const bool saveVelocityVTK = false;
//   // const bool saveLevelSetVTK = false;
//   // const int frequencyPlottingTheSolution = 10;




// //   const CutFEM_Parameter& h(Parameter::h);
// //   const CutFEM_Parameter& h_E(Parameter::meas);
// //
// //   double q0_0, q0_1, qp_0, qp_1;
// //   double intF = 0, intFnp = 0;
// //   int iter = 0, iterfig = 0;
// //   while( iter < GTime::total_number_iteration ) {
// //
// //     GTime::current_iteration = iter;
// //     const TimeSlab& In(Ih[iter]);
// //
// //     std::cout << " ------------------------------------------------------------- "<< std::endl;
// //     std::cout << " ------------------------------------------------------------- "<< std::endl;
// //     std::cout << " ITERATION \t : \t" << iter << std::endl;
// //     std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;
// //
// //
// //     ls_k.begin()->swap(ls_k[nbTime-1]);
// //     // computation of the nbTime interfaces
// //     for(int i=0;i<nbTime;++i) {
// //
// //       R tt = In.Pt(R1(qTime(i).x));
// //       ls_k[i].init(Lh_k, fun_levelSet, tt);
// //       ls[i].init(Lh, fun_levelSet, tt);
// //
// //       mapping[i] = new Mapping2(VelVh, ls_k[i]);
// //
// //       // projection(vel[i], vel[i]);
// //       vel[i].init(VelVh, fun_velocity, tt);
// //
// //       // projection(ls_k[i], ls[i]);
// //       interface.init(i,Th,ls[i].v);
// //
// //       // gnuplot::save(*interface[i], "interface"+to_string(i)+".dat");
// //       // if(i<nbTime-1) {
// //       //   LevelSet levelSet(ls_k[i], vel[i], vel[i+1], dt_levelSet);
// //       //   ls_k[i+1].init(levelSet.rhs);
// //       //   if(iter%frequencyReinitialization == 0 && i == 0)
// //       //    reinitialization.perform(ls_k[i], ls_k[i], *interface[i]);
// //       // }
// //     }
// //
// //     // vel[i].init(VelVh, fun_velocity, In,tt);
// //
// //
// //     // Create the Active Cut Mesh for insoluble surfactant
// //     Mesh2 cutTh(interface);
// // 	  FESpace2 cutVh(cutTh, interface, DataFE<Mesh2>::P1);
// //     cutVh.backSpace = &Lh_k;  // save backSpace to save solution
// //
// //     // FESpace2 cutVh2(cutTh, interface, FEvelocity);
// //     // cutVh2.backSpace = &VelVh;
// //     // FESpace2 VelVh  (Th, FEvelocity);
// //     // Fun_h vel2(cutVh2, In, fun_velocity);
// //
// //
// //     surfactant.initSpace(cutVh, In);
// //
// //     Rn datau0;
// //     surfactant.initialSolution(datau0);
// //
// //     KN_<double> datas0(datau0(SubArray(cutVh.NbDoF(),0)));
// //     if(iter == 0)interpolate(cutVh, datas0, fun_init_surfactant);
// //     Rn uh(datau0);
// //     Fun_h u0(cutVh, datau0);
// //     Normal n;
// //     FunTest s(cutVh,1), r(cutVh,1);
// //     cout << " Number of dof of the problem \t" << surfactant.nDoF << std::endl;
// //
// //
// //
// //     std::cout << " Begin assembly " << std::endl;
// //     double tt0 = CPUtime();
// //
// //
// //
// // #ifdef FORMULATION1
// //
// //     surfactant.addBilinear(
// //           // innerProduct(dt(s), r)
// //         + innerProduct(epsilon_S*gradS(s), gradS(r))
// //         , interface
// //         , In
// //     );
// //
// //     // for(int i=0;i<nbTime;++i) {      // computation of the curvature
// //     //   ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
// //     //   ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
// //     //   // ExpressionFunFEM<Mesh2> vx(vel,0,op_id);
// //     //   // ExpressionFunFEM<Mesh2> vy(vel,1,op_id);
// //     //   surfactant.addBilinear(
// //     //         innerProduct(dx(s)*vx + dy(s)*vy, r)
// //     //       + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
// //     //       , interface
// //     //       ,i , In
// //     //   );
// //     // }
// //
// //
// //     surfactant.addEdgeIntegral(
// //       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n)),
// //       In
// //     );
// //
// //     // surfactant.addBilinearFormInterface(
// //     //   innerProduct(1e-2*grad(s)*n, grad(r)*n)
// //     //   ,In
// //     // );
// //
// //     double cc = 1./In.T.mesure()*6.;
// //     surfactant.addBilinear(
// //       innerProduct(cc*s, r)
// //       , interface
// //       , 0
// //       , In
// //     );
// //     surfactant.addLinear(
// //       innerProduct(u0.expression(), cc*r)
// //       , interface
// //       , 0
// //       , In
// //     );
// //
// //     Fun_h funrhs(cutVh, In, fun_rhs);
// //     Fun_h funrhsp(cutVh, In);
// //     projection(funrhs, funrhsp, In, interface ,0);
// //     surfactant.addLinear(
// //       innerProduct(funrhs.expression(), r)
// //       , interface
// //       , In
// //     );
// //     intF = integralSurf(funrhsp, In, qTime) ;
// //     intFnp = integralSurf(funrhs, In, qTime);
// //
// //
// //     // for(int i=0;i<nbTime;++i) {
// //     //
// //     //   const R1 tq = In.map(qTime(i).x);
// //     //   Fun_h funrhs0(cutVh, fun_rhs, tq.x);
// //     //   Fun_h funrhsp(cutVh, fun_rhs, tq.x);
// //     //
// //     //   projection(funrhs0, funrhsp, *interface[i] );
// //     //
// //     //   // Fun_h funuh(cutVh, uh);
// //     //   // std::cout << integralSurf(funrhs0,0,i) << "\t"
// //     //   //           << integralSurf(funrhsp,0,i) << std::endl;
// //     //
// //     //   surfactant.addLinearFormInterface (
// //     //     innerProduct(funrhsp.expression(), r)
// //     //     , i, In
// //     //   );
// //     // }
// //
// // #else
// //
// //     surfactant.addBilinear(
// //         -innerProduct(s, dt(r))
// //         + innerProduct(epsilon_S*gradS(s), gradS(r))
// //         , interface
// //         , In
// //         // , mapping
// //     );
// //
// //     for(int i=0;i<nbTime;++i) {
// //       ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
// //       ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
// //       surfactant.addBilinear(
// //         - innerProduct(s, dx(r)*vx)
// //         - innerProduct(s, dy(r)*vy)
// //         , interface
// //         , i
// //         , In
// //         // , mapping
// //       );
// //     }
// //     // surfactant.addEdgeIntegral(
// //     //   innerProduct(1e-2*jump(grad(s).t()*n), jump(grad(r).t()*n)),
// //     //   In
// //     // );
// //
// //
// //     Fun_h funrhs (cutVh, In, fun_rhs);
// //     Fun_h funrhsp(cutVh, In);
// //     projection(funrhs, funrhsp, In, interface ,0);
// //     surfactant.addLinear (
// //       innerProduct(funrhs.expression(), r)
// //       , interface
// //       , In
// //       // , mapping
// //     );
// //     intF   = integralSurf(funrhsp, In, qTime) ;
// //     intFnp = integralSurf(funrhs , In, qTime);
// //
// //
// //     double ccend = 1./In.T.mesure()*1./qTime[lastQuadTime].a;;
// //
// //     surfactant.addBilinear(
// //       innerProduct(ccend*s, r)
// //       , interface
// //       , lastQuadTime
// //       , In
// //       // , mapping
// //     );
// //     double cc0 = 1./In.T.mesure()*1./qTime[0].a;;
// //     surfactant.addLinear(
// //       innerProduct(u0.expression(), cc0*r)
// //       , interface
// //       , 0
// //       , In
// //       // , mapping
// //     );
// //
// //     FunTest D2un = grad(grad(r)*n)*n;
// //
// //     surfactant.addEdgeIntegral(
// //       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
// //       // + innerProduct(1e-2*h_E*jump(D2un), h_E*jump(D2un)),
// //       , In
// //       // , mapping
// //     );
// //     // if stab interface, can have h in both stabilization
// //     // surfactant.addBilinear(
// //     //   innerProduct(1e-2*h_E*grad(s)*n, grad(r)*n)
// //     //   , interface
// //     //   , In
// //     // );
// //
// //
// //
// // #endif
// //     surfactant.solve();
// //
// //     KN_<double> dw(surfactant.rhs(SubArray(surfactant.nDoF, 0)));
// //     uh = dw;
// //     surfactant.saveSolution(uh);
// //
// //
// //
// //     // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
// //     {
// //       cout << "\n Features of the drop " << std::endl;
// //       Fun_h funuh(cutVh, uh);
// //
// //       Rn sol2(cutVh.NbDoF(), 0.);
// //       Fun_h funsol(cutVh, sol2);
// //       sol2 += uh(SubArray(cutVh.NbDoF(), 0));
// //       double q_0 = integralSurf(funsol, 0);
// //       sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
// //       double q_1 = integralSurf(funsol, 0, lastQuadTime);
// //       if(iter==0){ q0_0 = q_0;q0_1 = q_1;
// //                    qp_1 = q_1;
// //                    q0_1 = integralSurf(u0, 0);
// //                  }
// //
// //
// //       // cout << " q0 - qt_1 \t" << q_1-qp_1 << std::endl;
// //       // cout << " int F \t" << intF << std::endl;
// //
// //       outputData << setprecision(10);
// //       outputData << GTime::current_time() << "\t"
// //                  << (q_1-qp_1) << "\t"
// //                  << fabs(q_1-qp_1) << "\t"
// //                  << intF << "\t"
// //                  << intFnp << "\t"
// //                  << ((q_1 -qp_1) - intFnp) << "\t"
// //                  << q0_1 - q_1 << "\t"
// //                  << std::endl;
// //       qp_1 = q_1;
// //     }
// //
// //
// //     {
// //       Fun_h funsol(cutVh, fun_sol_surfactant, GTime::current_time());
// //       funsol.v -= uh(SubArray(cutVh.NbDoF(), 0));
// //
// //       Rn sol2(cutVh.NbDoF(), 0.);
// //       sol2 += uh(SubArray(cutVh.NbDoF(), 0));
// //       // sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
// //
// //       ExpressionFunFEM<Mesh> sol(funsol, 0, op_id);
// //       // cout << " error    ->    " << sqrt(integralSurf(sol*sol, cutVh, 0)) << std::endl;
// //       // Fun_h funuh(cutVh, uh);
// //       Fun_h funuh(cutVh, sol2);
// //       // cout << " error    ->    " << L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1) << std::endl;
// //       errL2 =  L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1);
// //       std::cout << " || u-uex||_2 = " << errL2 << std::endl;
// //     }
// //
// //
// //   //  -----------------------------------------------------
// //   //                     PLOTTING
// //   //  -----------------------------------------------------
// //   if(MPIcf::IamMaster() && iter%frequencyPlottingTheSolution == 0 && writeVTKFiles){
// //     std::cout << " Plotting " << std::endl;
// //
// //     if(saveSurfactantVTK) {
// //       KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
// //       Fun_h sol(cutVh, sols);
// //
// //       // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
// //       Paraview2 writer(cutVh, "surfactant_"+to_string(iterfig)+".vtk");
// //
// //       writer.add(sol , "surfactant", 0, 1);
// //       writer.add(ls[0], "levelSet",  0, 1);
// //       writer.add(ls[2], "levelSet2", 0, 1);
// //     }
// //     // if(saveLevelSetVTK) {
// //     //   Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
// //     //   writerLS.add(ls[0], "levelSet", 0, 1);
// //     // }
// //
// //     iterfig++;
// //   }
// //   for(int i=0;i<nbTime;++i) delete mapping[i];
// //   iter++;
// //
// //   }
// //
// //   myCout << 4./nx << "\t" << errL2 << std::endl;
// //   nx *= 2;
// //   ny *= 2;
// // }
// // return 0;
// // }
