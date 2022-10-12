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
        return -P.x - P.y + c0 + Epsilon;
    }

    R fun_levelSet(const R2 P, const int i) {
        return -P.x - P.y + c0;
    }

    R fun_boundary(const R2 P, int elementComp, double t) {
        return sin(pi*(P.x + P.y - 4*t));
    }

    R fun_uBulk(const R2 P, int elementComp, int domain, double t) {
        // return 2*sin(pi*(P.x+P.y-3*t));
        if(domain == 0) return sin(pi*(P.x + P.y - 4*t));
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
#define frachon1
// Set scheme for the dg method (options: "conservative", "classical")
#define conservative
// Set stabilization method (options: "fullstab", "macro") 
#define fullstab   

#define use_n


#ifdef frachon1
    using namespace Frachon1;
#elif defined(frachon2)
    using namespace Frachon2;
#endif
 
int main(int argc, char** argv) {
    
    // Mesh settings and data objects
    const size_t iterations = 4;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 10, ny = 10;       // starting mesh size
    double h = 0.125;             // starting mesh size
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
    std::array<double, iterations> dts;


    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {
        
        #ifdef use_h
        nx = (int)(lx/h)+1, ny = (int)(ly/h)+1;
        #elif defined(use_n)
        h = lx/(nx-1);
        #endif

        // Define background mesh
        Mesh Th(nx, ny, -1., -1., lx, ly);

    //// Parameters
        
        // Mesh size
        int divisionMeshSize = 2*3*sqrt(10);
        double dT = h/divisionMeshSize;
        hs.at(j) = h;
        dts.at(j) = dT;

        // Time
        double tfinal = 1.;            // Final time
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
        std::cout << "ny = " << ny << std::endl;
        std::cout << "dT = " << dT << std::endl;

        // Penalty parameter for outer boundary
        double lambdaB = 3;    // coefficient on the outer boundary

        CutFEMParameter lambda(0., 1.);     // lambda2 = 0, lambda1 = 1,  1+lambda2-lambda1 = 0
        //CutFEMParameter lambda(1./3, 1.);

        CutFEMParameter lambdaE(3, 2);  //  ||beta||_inf in Omega_1/Omega_2

        // DG stabilization parameters
        double tau20 = 5e-1, tau21 = 5e-1;      // bulk

        
        // Background FE Space, Time FE Space & Space-Time Space
        FESpace2 Vh2(Th, DataFE<Mesh>::P2);        // higher order space for interpolation
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
        double inflow = 0, outflow = 0;      // hold integrals of rhs and Neumann bcs
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

                // FIXME: This version of the code has a memory leak (it doesn't delete the interface pointers)
                interface.init(i, Th, ls[i]);
                
            }
            
            ActiveMesh<Mesh> Khi(Th, interface);
        
            CutSpace Wh(Khi, Vh);
            
            convdiff.initSpace(Wh, In);

            // Velocity field
            Lagrange2 FEvelocity(0);
            FESpace2 VelocitySpace(Th, FEvelocity);
            CutSpace VelSpace(Khi, VelocitySpace);
            Fun_h vel(VelSpace, fun_velocity);

            Normal n;

            // Exact solution
            Fun_h g(Vh2, In, fun_boundary);    // create an FE-function of the exact bulk solution Omega1

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega1 U Omega2                 
            
            // Data for initial solution
            Rn data_u0;                           // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));        // initial data bulk
            
            if (iter == 0) interpolate(Wh, data_B0, fun_uBulkInit);
            
            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);

        //// Assembling linear and bilinear forms
        
        #ifdef classical

            // Bulk terms
            convdiff.addBilinear(
                + innerProduct(dt(u), v)
                - innerProduct(u, (vel.expression()*grad(v)))
                , Khi
                , In
            );

            // Time penalty term bulk LHS
            convdiff.addBilinear(
                + innerProduct(u, v)
                , Khi
                , 0
                , In
            );

            
            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), v)
                , Khi
                , 0
                , In
            );


        #elif defined(conservative)
            
            // Bulk terms
            convdiff.addBilinear(
                - innerProduct(u, dt(v))
                - innerProduct(u, (vel.expression()*grad(v)))
                , Khi
                , In
            );

            // Time penalty term bulk LHS
            convdiff.addBilinear(
                + innerProduct(u, v)
                , Khi
                , (int)lastQuadTime
                , In
            );

            
            // Time penalty term bulk RHS
            convdiff.addLinear(
                + innerProduct(b0h.expression(), v)
                , Khi
                , 0
                , In
            );
        #endif

            
            // Inner edges bulk convection
            convdiff.addBilinear(
                + innerProduct(average((vel*n)*u), jump(v))         
                + innerProduct(0.5*fabs(vel*n)*jump(u), jump(v))  
                , Khi
                , INTEGRAL_INNER_EDGE_2D                 
                , In
            );

            // Interface terms q
            convdiff.addBilinear(
                + jump(innerProduct((vel*n)*u,v))
                + innerProduct(jump((vel*n)*u), jump(lambda*v))
                , interface
                , In
            );

            
            convdiff.addBilinear(
                + innerProduct((vel*n)*u, v)
                , Khi
                , INTEGRAL_BOUNDARY
                , In
                , {2, 3}          
            );  
            
           
            // convdiff.addBilinear(
            //     + innerProduct((vel*n)*u, 0.5*v) 
            //     + innerProduct(0.5*fabs(vel*n)*u, v)  
            //     , Khi
            //     , INTEGRAL_BOUNDARY
            //     , In
            //     , (std::list<int>){1,4}
            // );
            
            // convdiff.addLinear(
            //     - innerProduct(g.expression(), (vel*n)*(0.5*v))
            //     + innerProduct(g.expression(), fabs(vel*n)*0.5*v) 
            //     , Khi
            //     , INTEGRAL_BOUNDARY
            //     , In
            //     , (std::list<int>){1,4}          // 1 - bottom, 4 - left
            // );

            // NOTE: Below is equivalent to the above two when using upwind flux

            convdiff.addLinear(
                - innerProduct(g.expression(), (vel*n)*v) 
                , Khi
                , INTEGRAL_BOUNDARY
                , In
                , (std::list<int>){1,4}
            );
            
            // Stabilization
        #ifdef macro    
            TimeMacroElement<Mesh> TimeMacro(Khi, qTime, 0.125);

            if (iterations == 1) {
                Paraview<Mesh> writerMacro(Th, pathOutputFigures + "Th" + to_string(iter+1) + ".vtk");
                writerMacro.add(ls[0], "levelSet0.vtk", 0, 1);
                writerMacro.add(ls[1], "levelSet1.vtk", 0, 1);
                writerMacro.add(ls[2], "levelSet2.vtk", 0, 1);  

                // domain = 0, 

                writerMacro.writeFaceStab(Khi, 0, pathOutputFigures + "FullStabilization" + to_string(iter+1) + ".vtk");
                writerMacro.writeActiveMesh(Khi, pathOutputFigures + "ActiveMesh" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 0, pathOutputFigures + "macro_dom0_" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroInnerEdge(TimeMacro, 0, pathOutputFigures + "macro_inner_edge_dom0_" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroOutterEdge(TimeMacro, 0, pathOutputFigures + "macro_outer_edge_dom0_" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroElement(TimeMacro, 1, pathOutputFigures + "macro_dom1_" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroInnerEdge(TimeMacro, 1, pathOutputFigures + "macro_inner_edge_dom1_" + to_string(iter+1) + ".vtk");
                writerMacro.writeMacroOutterEdge(TimeMacro, 1, pathOutputFigures + "macro_outer_edge_dom1_" + to_string(iter+1) + ".vtk");
            }
            
            // Stabilization of the bulk 
            convdiff.addFaceStabilization(
                + innerProduct(tau20*jump(u), jump(v))
                + innerProduct(h*h*tau21*jump(grad(u)), jump(grad(v)))
                , Khi
                , In
                , TimeMacro
            );

        #elif defined(fullstab)
            convdiff.addFaceStabilization(
                + innerProduct(tau20*jump(u), jump(v))
                + innerProduct(h*h*tau21*jump(grad(u)), jump(grad(v)))
                , Khi
                , In
            );
        #endif
            
            // Solve linear system
            convdiff.solve("mumps");
            
            data_u0 = convdiff.rhs_;
            convdiff.saveSolution(data_u0);

            Expression2 bhexp(b0h, 0, op_id);
            Expression2 ghexp(g, 0, op_id);

            inflow = integral(Khi, In, (vel*n)*ghexp, INTEGRAL_BOUNDARY, qTime, {1,4});
            outflow = integral(Khi, In, (vel*n)*bhexp, INTEGRAL_BOUNDARY, qTime, {2,3});

            if (iterations == 1) {
                  
                Fun_h funuh(Wh, data_u0);

                Rn sol2(Wh.NbDoF(), 0.);
                Fun_h funsol(Wh, sol2);
                sol2 += data_u0(SubArray(Wh.NbDoF(), 0));
                double q_0 = integral(Khi, funsol, 0, 0);
                sol2 += data_u0(SubArray(Wh.NbDoF(), Wh.NbDoF()));
                double q_1 = integral(Khi, funsol, 0, lastQuadTime);
                
                if(iter==0) { 
                    q0_0 = q_0; q0_1 = q_1;
                    qp_1 = q_1;
                    q0_1 = integral(Khi, b0h, 0, 0);
                }
            
                outputData << setprecision(10);
                outputData << GTime::current_time() << ","
                        << (q_1 - qp_1) << ","
                        << inflow << ","
                        << outflow << ","
                        << ((q_1 - qp_1) + inflow + outflow)
                        << std::endl;

                qp_1 = q_1;
                
            }


            R errBulk =  L2normCut(b0h, fun_uBulk, GTime::current_time(), 0, 1);
            std::cout << std::endl;
            std::cout << " L2 Error \t : \t" << errBulk << std::endl;

            //gnuplot::save(cutThTime,"cutThTime.dat");
            //gnuplot::save(Th, "Th.dat");
            
            errorsBulk.at(j) = errBulk;

            if ((iterations == 1) && MPIcf::IamMaster()) {
                
                Paraview<Mesh> writer(Khi, pathOutputFigures + "Bulk" + to_string(iter + 1)+"DG.vtk");
                writer.add(b0h, "bulk", 0, 1);
                Fun_h uBex(Wh, fun_uBulk, GTime::current_time());
                writer.add(uBex,"bulk_exact", 0, 1);
                writer.add(ls[0], "levelSet", 0, 1);
            }

            if (iterations > 1 && iter == GTime::total_number_iteration-1) outputData << h << "," << dT << "," << errBulk << std::endl;

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

