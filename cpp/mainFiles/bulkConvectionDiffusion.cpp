/*  
   We consider the stationary problem solved in statConvectionDiffusion
   but with an added time component to the solution. We consider here a
   convection-diffusion problem.

   The main motivation for studying this example is that we have a working
   method for the corresponding stationary problem, and it has Dirichlet BCs.

   Problem:
   Find u on bulk domain Omega inside interface Gamma. 
   u(t,x,y) = 2*exp(1-x^2-y^2)*(3x^2y - y^3)*exp(-t/4)
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "../util/redirectOutput.hpp"
#include "baseCutProblem.hpp"
#include "../num/gnuplot.hpp"
#include "projection.hpp"
#include <typeinfo>


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

using namespace Flower;

// We consider a two-dimensional problem
const int d = 2;
typedef Mesh2 Mesh;
typedef FESpace2 FESpace;
typedef TestFunction<d> FunTest;
typedef FunFEM<Mesh2> Fun_h;

/*
CLASSIC1 : No integration by parts on convection term
CLASSIC2 : Integration by parts on full convection term, no term added to make anti-symmetric
CLASSIC3 : Integration by parts on half of the convection term, term added to make anti-symmetric
*/

#define CLASSIC2


/* Main DG Scheme */

void DG(int argc, char** argv) {
    
    // Mesh settings and data objects
    int iterations = 1;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 20, ny = 20;       // starting mesh size
    double lx = 6., ly = 6.;    // domain length

    // Initialize MPI
    MPIcf cfMPI(argc,argv);

    // Paths to store data
    std::string pathOutputFolder = "../outputFiles/BulkConvectionDiffusion/data/";
    std::string pathOutputFigures = "../outputFiles/BulkConvectionDiffusion/paraview/";

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::experimental::filesystem::create_directories(pathOutputFolder);
        std::experimental::filesystem::create_directories(pathOutputFigures);
    }

    // Data file to hold problem data
    std::ofstream outputData(pathOutputFolder+"dataDGClassic.dat", std::ofstream::out);

    // Arrays to hold data
    vector<double> errorsBulk(iterations);      // vector to hold bulk errors
    vector<double> hs(iterations);              // vector to hold mesh sizes

    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {

        // Define background mesh
        Mesh Th(nx, ny, -1.5, -1.5, lx, ly);

        // =============================================================== //
        // ------------------------ Parameters --------------------------- //
        // =============================================================== //
        
        // Mesh size
        double h = lx/nx;
        hs.at(j) = h;
        int divisionMeshSize = 2;

        // Time
        double dT = h/divisionMeshSize; // Time step size
        double tfinal = 1.5;           // Final time
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

        // General problem parameters
        double A2 = 0.5;        // diffusion coefficients
        double kappa2 = 0.5, kappa02 = 2;
        double kappaTilde2 = kappa2/kappa02;

        // Constants for penalty terms
        double tau_a2 = 300;        // diffusion penalty scaling
        double tau_b2 = 0;        // convection penalty scaling
        // Bulk penalties
        double lambdaA = kappaTilde2*A2*tau_a2/h;      // diffusion terms
        double lambdaB = kappaTilde2*A2*tau_b2;

        // Penalty parameter for outer boundary
        double lambda = A2*10/h;    // coefficient on the outer boundary

        // Stabilization parameters
        double tau20 = A2, tau21 = 0.1*A2;      // bulk


        // =============================================================== //
        // ------------------------ Spaces ------------------------------- //
        // =============================================================== //

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
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

        // Velocity field
        Lagrange2 FEvelocity(2);
        FESpace2 VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        FESpace2 Lh(Th, DataFE<Mesh2>::P1);
        double dt_levelSet = dT/(nbTime-1);
        vector<Fun_h> ls(nbTime);
        for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);

        // Declare time dependent interface
        TimeInterface2 interface((int)nbTime);
        
        // Bulk-Surface Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;
        
        int iter = 0;
        // Iterate over time-slabs
        while (iter < GTime::total_number_iteration) {
            
            GTime::current_iteration = iter;
            const TimeSlab &In(Ih[iter]);

            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ITERATION \t : \t" << iter + 1 << "/" << GTime::total_number_iteration << std::endl;
            std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;
    
            
            ls.begin()->swap(ls[nbTime - 1]);

            // computation of the interface in each of the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

                interface.init(i, Th, ls[i].v);

                if (i < lastQuadTime) {
                    LevelSet2 levelSet(ls[i], vel, vel, dt_levelSet);
                    ls[i+1].init(levelSet.rhs);
                }
            }

            // Interface mesh
            Mesh cutThTime(interface);

            // gnuplot::save(Th, "Th.dat");
            // gnuplot::save(cutThTime, "cutThTime.dat");

            // Cut FE space
            CutFESpace2 Wh(Vh, interface, {1});    // -1 denotes the inner bulk domain
    
            convdiff.initSpace(Wh, In);


            // gnuplot::save(Th, "Th.dat");
            // gnuplot::save(cutThTime,  "cutThTime.dat");
            // gnuplot::save(*interface[0],  "interface.dat");

            // =============================================================== //
            // ------------------ Variational Form Objects ------------------- //
            // =============================================================== //

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Wh, In, fun_rhsBulk);

            Fun_h g(Wh, In, fun_uBulk);    // create an FE-function of the exact bulk solution Omega1

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega2

            // =============================================================== //
            // ---------------------- Solution arrays ------------------------ //
            // =============================================================== //                    
            
            // Data for initial solution
            Rn data_u0;                           // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));        // initial data bulk
            
            if (iter == 0) interpolate(Wh, data_B0, fun_uBulkInit);
            
            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);
            
            // Plot initial solution in paraview
            if (iter == 0 && MPIcf::IamMaster()) {

                Paraview2 writerInitial(Wh, ls[0], pathOutputFigures + "BulkInitial.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
                
                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);
                
                writerInitial.add(uBex, "bulk_exact", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
                
            }
            
            // =============================================================== //
            // -------------------------- Assembly --------------------------- //
            // =============================================================== // 
            
            // ------------ TIME TERMS -------------- //

            // Partial derivative in time, bulk     (3rd term in (1.52))
            convdiff.addBilinear(
                + innerProduct(dt(u), kappaTilde2*v)
                , In
            );
            
            // ----- Bilinear form for diffusion ----- //

            // Integral on element for bulk variables (K_{h,1} and K_{h,2})
            convdiff.addBilinear(
                + innerProduct(kappaTilde2*A2*grad(u), grad(v))
                , In
            );

            // integral on Inner Edges for bulk variables (i.e. E_{h,1} and E_{h,2})
            convdiff.addEdgeIntegral(
                - innerProduct(kappaTilde2*A2*average(grad(u)*n), jump(v))      
                - innerProduct(kappaTilde2*A2*jump(u), average(grad(v)*n))     
                + innerProduct(lambdaA*jump(u),jump(v))                     // bulk penalty
                , In
            );

            // ----- Bilinear form for convection ------ //
#ifdef CLASSIC1
            // Bulk convection
            convdiff.addBilinear(
                + innerProduct((vel.expression(),grad(u)), kappaTilde2*v)
                , In
            );
#elif defined(CLASSIC2) 
            // Bulk convection
            convdiff.addBilinear(
                - innerProduct(kappaTilde2*u, (vel.expression(),grad(v)))
                , In
            );

            // Inner edges bulk convection
            convdiff.addEdgeIntegral(
                + innerProduct(average((vel*n)*u), kappaTilde2*jump(v))         
                //- innerProduct(jump((vel*n)*u), kappaTilde2*average(v))   
                + innerProduct(lambdaB*fabs(vel*n)*jump(u), jump(v))                   
                , In
            );

#elif defined(CLASSIC3)
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


#endif
            // Stabilization of the bulk 
            convdiff.addFaceStabilization(
                + innerProduct(1./h*tau20*jump(u), kappaTilde2*jump(v))
                + innerProduct(h*tau21*jump(grad(u)), kappaTilde2*jump(grad(v)))
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
            
            // ----- Add Dirichlet boundary conditions ------
            convdiff.addBilinear(
                    - innerProduct(kappaTilde2*A2*grad(u)*n, v)   // from IBP
                    - innerProduct(u, kappaTilde2*A2*grad(v)*n)   // added to make symmetric
                    + innerProduct(u, kappaTilde2*lambda*v)       // added penalty
                    , interface
                    , In
            );
            
    
            // RHS on Gamma
        #ifdef CLASSIC1
            convdiff.addLinear(
                    - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)   
                    + innerProduct(g.expression(), kappaTilde2*lambda*v)
                    , interface
                    , In
            );
        #elif defined(CLASSIC2)
            convdiff.addLinear(
                    - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)   
                    + innerProduct(g.expression(), kappaTilde2*lambda*v)
                    - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)
                    , interface
                    , In
            );
        #elif defined(CLASSIC3)
            convdiff.addLinear(
                    - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)   
                    + innerProduct(g.expression(), kappaTilde2*lambda*v)
                    - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)*0.5
                    , interface
                    , In
            );
        #endif

            // Add RHS on bulk
            convdiff.addLinear(
                    innerProduct(f.expression(), kappaTilde2*v)
                    , In
            );
            
            // Solve linear system
            convdiff.solve();
            
            data_u0 = convdiff.rhs;
            convdiff.saveSolution(data_u0);

            R errBulk =  L2normCut(b0h, fun_uBulkD, GTime::current_time(), 0, 1);
            std::cout << " || uBulk-uBulkex||_2 = " << errBulk << std::endl;

            gnuplot::save(cutThTime,"cutThTime.dat");
            gnuplot::save(Th, "Th.dat");
            
            errorsBulk.at(j) = errBulk;

            if ((iterations == 1) && MPIcf::IamMaster()) {
                
                Paraview2 writer(Wh, ls[0], pathOutputFigures + "Bulk" + to_string(iter + 1)+"DG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, GTime::current_time());
                Fun_h fB(Wh, fun_rhsBulk, GTime::current_time());
                writer.add(uBex,"bulk_exact", 0, 1);
                writer.add(fB,"bulk_rhs", 0, 1);
                writer.add(ls[0], "levelSet", 0, 1);
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
    

}


/* Main CG Scheme */

void CG(int argc, char** argv) {
    
    // Mesh settings and data objects
    int iterations = 5;         // number of mesh refinements   (set to 1 to run only once and plot to paraview)
    int nx = 15, ny = 15;       // starting mesh size
    double lx = 3., ly = 3.;    // domain length

    // Initialize MPI
    MPIcf cfMPI(argc,argv);

    // Paths to store data
    std::string pathOutputFolder = "../outputFiles/BulkConvectionDiffusion/data/";
    std::string pathOutputFigures = "../outputFiles/BulkConvectionDiffusion/paraview/";

    // Create directory if not already existent
    if (MPIcf::IamMaster()) {
        std::experimental::filesystem::create_directories(pathOutputFolder);
        std::experimental::filesystem::create_directories(pathOutputFigures);
    }

    // Data file to hold problem data
    //std::ofstream outputData(pathOutputFolder+"dataDGClassic.dat", std::ofstream::out);

    // Arrays to hold data
    vector<double> errorsBulk(iterations);      // vector to hold bulk errors
    vector<double> hs(iterations);              // vector to hold mesh sizes

    // Iterate over mesh sizes
    for (int j=0; j<iterations; ++j) {

        // Define background mesh
        Mesh Th(nx, ny, -1.5, -1.5, lx, ly);

        // =============================================================== //
        // ------------------------ Parameters --------------------------- //
        // =============================================================== //
        
        // Mesh size
        double h = lx/nx;
        hs.at(j) = h;
        int divisionMeshSize = 2;

        // Time
        double dT = h/divisionMeshSize; // Time step size
        double tfinal = 0.25;           // Final time
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

        // General problem parameters
        double A2 = 0.5;        // diffusion coefficients
        double kappa2 = 0.5, kappa02 = 2;
        double kappaTilde2 = kappa2/kappa02;


        // Penalty parameter for outer boundary
        double lambda = A2*1e1/h;    // coefficient on the outer boundary

        // Stabilization parameters
        double tau20 = A2, tau21 = 0.1*A2;      // bulk


        // =============================================================== //
        // ------------------------ Spaces ------------------------------- //
        // =============================================================== //

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        FESpace2 Vh(Th, DataFE<Mesh>::P1);        // discontinuous basis functions
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
        for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);

        // Declare time dependent interface
        TimeInterface2 interface((int)nbTime);
        
        // Bulk-Surface Convection-Diffusion Problem Object
        CutFEM<Mesh> convdiff(qTime);

        std::cout << "Number of time slabs \t : \t " << GTime::total_number_iteration << std::endl;
        
        int iter = 0;
        // Iterate over time-slabs
        while (iter < GTime::total_number_iteration) {
            
            GTime::current_iteration = iter;
            const TimeSlab &In(Ih[iter]);

            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ------------------------------------------------------------- " << std::endl;
            std::cout << " ITERATION \t : \t" << iter + 1 << "/" << GTime::total_number_iteration << std::endl;
            std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;
    
            
            ls.begin()->swap(ls[nbTime - 1]);

            // computation of the interface in each of the three quadrature points
            for (int i = 0; i < nbTime; ++i) {

                interface.init(i, Th, ls[i].v);

                if (i < lastQuadTime) {
                    LevelSet2 levelSet(ls[i], vel, vel, dt_levelSet);
                    ls[i+1].init(levelSet.rhs);
                }
            }

            // Interface mesh
            Mesh cutThTime(interface);

            // Cut FE space
            CutFESpace2 Wh(Vh, interface, {1});    // -1 denotes the inner bulk domain
    
            convdiff.initSpace(Wh, In);


            // gnuplot::save(Th, "Th.dat");
            // gnuplot::save(cutThTime,  "cutThTime.dat");
            // gnuplot::save(*interface[0],  "interface.dat");

            // =============================================================== //
            // ------------------ Variational Form Objects ------------------- //
            // =============================================================== //

            // Objects needed for the weak form
            Normal n;
            Tangent t;

            // Right hand side functions
            Fun_h f(Wh, In, fun_rhsBulk);

            Fun_h g(Wh, In, fun_uBulk);    // create an FE-function of the exact bulk solution Omega1

            // Test and Trial functions
            FunTest u(Wh, 1), v(Wh, 1); // Omega2

            // =============================================================== //
            // ---------------------- Solution arrays ------------------------ //
            // =============================================================== //                    
            
            // Data for initial solution
            Rn data_u0;                           // initial data total
            convdiff.initialSolution(data_u0);
            KN_<R> data_B0(data_u0(SubArray(Wh.NbDoF(), 0)));        // initial data bulk
            
            if (iter == 0) interpolate(Wh, data_B0, fun_uBulkInit);
            
            // Make function objects to use in innerProducts
            Fun_h b0h(Wh, data_B0);
            
            // Plot initial solution in paraview
            if (iter == 0 && MPIcf::IamMaster()) {

                Paraview2 writerInitial(Wh, ls[0], pathOutputFigures + "BulkInitial.vtk");
                writerInitial.add(b0h, "bulk", 0, 1);
                writerInitial.add(ls[0], "levelSet", 0, 1);
                
                // Add exact solutions
                Fun_h uBex(Wh, fun_uBulkD, 0.);
                
                writerInitial.add(uBex, "bulk_exact", 0, 1);
                
            }
            
            // =============================================================== //
            // -------------------------- Assembly --------------------------- //
            // =============================================================== // 
            
            // ------------ TIME TERMS -------------- //

            // Partial derivative in time, bulk     (3rd term in (1.52))
            convdiff.addBilinear(
                + innerProduct(dt(u), kappaTilde2*v)
                , In
            );
            
            // Integral on element for bulk variables (K_{h,1} and K_{h,2})
            convdiff.addBilinear(
                + innerProduct(kappaTilde2*A2*grad(u), grad(v))
                , In
            );

            // Bulk convection
            convdiff.addBilinear(
                - innerProduct(kappaTilde2*u, (vel.expression(),grad(v)))
                , In
            );

            // Stabilization of the bulk 
            convdiff.addFaceStabilization(
                + innerProduct(h*tau21*jump(grad(u)), kappaTilde2*jump(grad(v)))
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
            
            // ----- Add Dirichlet boundary conditions ------
            convdiff.addBilinear(
                    - innerProduct(kappaTilde2*A2*grad(u)*n, v)   // from IBP
                    - innerProduct(u, kappaTilde2*A2*grad(v)*n)   // added to make symmetric
                    + innerProduct(u, kappaTilde2*lambda*v)       // added penalty
                    , interface
                    , In
            );
            //getchar();
            // RHS on outer boundary
            convdiff.addLinear(
                    - innerProduct(g.expression(), kappaTilde2*A2*grad(v)*n)   
                    + innerProduct(g.expression(), kappaTilde2*lambda*v)
                    - innerProduct(g.expression(), (vel*n)*kappaTilde2*v)
                    , interface
                    , In
            );

            // convdiff.addLinear(
            //     + innerProduct(grad(g.expression())*n*A2, v)
            //     , interface
            //     , In
            // );

            // Add RHS on bulk
            convdiff.addLinear(
                    innerProduct(f.expression(), kappaTilde2*v)
                    , In
            );
            
            // Solve linear system
            convdiff.solve();
            
            data_u0 = convdiff.rhs;
            convdiff.saveSolution(data_u0);

            R errBulk =  L2normCut(b0h, fun_uBulkD, GTime::current_time(), 0, 1);
            std::cout << " || uBulk-uBulkex||_2 = " << errBulk << std::endl;

            
            errorsBulk.at(j) = errBulk;

            if ((iterations == 1) && MPIcf::IamMaster()) {
                
                Paraview2 writer(Wh, ls[0], pathOutputFigures + "Bulk" + to_string(iter + 1)+"CG.vtk");
                writer.add(b0h, "bulk", 0, 1);

                Fun_h uBex(Wh, fun_uBulk, GTime::current_time());
                Fun_h fB(Wh, fun_rhsBulk, GTime::current_time());
                writer.add(uBex,"bulk_exact", 0, 1);
                writer.add(fB,"bulk_rhs", 0, 1);
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
    

}


int main(int argc, char** argv) {

    DG(argc, argv);

    return 0;
}