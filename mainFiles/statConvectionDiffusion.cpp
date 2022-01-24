#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#  include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "baseCutProblem.hpp"
#include "../num/gnuplot.hpp"
#include "projection.hpp"
#include <typeinfo>

double fun_rhs1(const R2 P) {
    return -exp(- pow(P.x,2) - pow(P.y,2) + 1)*(12*pow(P.x,4)*P.y + 3*pow(P.x,3) + 8*pow(P.x,2)*pow(P.y,3) - 48*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 4*pow(P.y,5) + 16*pow(P.y,3));
}
double fun_rhs2(const R2 P) {
    return -2*exp(- pow(P.x,2) - pow(P.y,2) + 1)*(6*pow(P.x,4)*P.y + 3*pow(P.x,3) + 4*pow(P.x,2)*pow(P.y,3) - 24*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 2*pow(P.y,5) + 8*pow(P.y,3));
}

// need this function to compute on cutDomain
double fun_rhsBulk(const R2 P, int elementComp, int dom) {
    if(dom == 0) return fun_rhs1(P);
    else return fun_rhs2(P);
}

// ok for the interface problem
double fun_rhs0(const R2 P, int elementComp) {
    return 3*(- pow(P.x,5) + 2*pow(P.x,3)*pow(P.y,2) + 9*pow(P.x,2)*P.y + 3*P.x*pow(P.y,4) - 3*pow(P.y,3))/(pow(P.x,2) + pow(P.y,2));
}

// for bulk
double fun_uBulk(const R2 P,  int elementComp, int domain) {
    if(domain == 0) return     exp(1.0-P.x*P.x-P.y*P.y)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
    else            return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
}
// for surface
double fun_uSurface(const R2 P,  int elementComp) {
    return 3.0*P.x*P.x*P.y - pow(P.y,3);
}

double fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x)*(P.x) + (P.y)*(P.y)) - 1.0;
}

R2 fun_beta(const R2 P) {
    R2 beta(P.y, -P.x);
    return beta;
}
// Or, like in bulkSurface:
R fun_velocity(const R2 P, const int i){
    if(i == 0) return P.y;
    else       return -P.x;
}

// 2 DIMENSIONAL PROBLEM
// =====================================================
typedef Mesh2 Mesh;
typedef FESpace2 FESpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;

void assemble_forms(int nx, int ny, int argc, char** argv) {

}

int main(int argc, char** argv ) {

    // =====================================================
    // TODO: QUESTIONS:
    //=====================================================

    // INITIALIZE MPI
    // =====================================================
    MPIcf cfMPI(argc,argv);

    int nx = 28, ny = 28;
    double lx = 3., ly = 3.;
    double h = lx/nx;
    std::cout << "h = " << h << std::endl;

    // Problem parameters
    double A0 = 1.0, A1 = 1.0, A2 = 0.5;
    double kappa1 = 2.0, kappa2 = 0.5, kappa01 = 1, kappa02 = 2;
    double kappaTilde1 = kappa1/kappa01, kappaTilde2 = kappa2/kappa02;
    CutFEM_Parameter kappaTilde("kappaTilde", kappaTilde1, kappaTilde2);
    CutFEM_Parameter kappaTildeA("kappaTildeA", kappaTilde1*A1, kappaTilde2*A2);

    // Penalty parameters, arbitrarily chosen at this time.
    double tau_a0 = 10000, tau_b0 = 0;
    double tau_a1 = 10000, tau_b1 = 0;
    double tau_a2 = 10000, tau_b2 = 0;
    CutFEM_Parameter penalty1("penalty1", kappaTilde1*A1*tau_a1/h, kappaTilde2*A2*tau_a2/h);
    CutFEM_Parameter penalty2("penalty2", kappaTilde1*tau_b1, kappaTilde2*tau_b2);

    /*
     tau_b0 = 10
     || u  - u_ex ||_2  = 	0.0016111
     || u0 - u0_ex||_2  = 	0.00116648

      tau_b0 = 1
     || u  - u_ex ||_2  = 	0.0016111
     || u0 - u0_ex||_2  = 	0.00116648

     tau_b0 = 0
    || u  - u_ex ||_2  = 	0.00161297
    || u0 - u0_ex||_2  = 	0.00116649

     */

    // Stabilization parameters
    double tau11 = 0.1*A1, tau10 = A1, tau20 = A2, tau21 = 0.1*A2;
    double tau00 = A0, tau01 = A0, tau02 = 0.1*A0;
    CutFEM_Parameter tau_i0("taui0", tau10, tau20);
    CutFEM_Parameter tau_i1("taui1", tau11, tau21);


    // CONSTRUCTION OF THE MESH
    // =====================================================

    Mesh Th(nx, ny, -1.5, -1.5, lx, ly);


    // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
    // =====================================================
    FESpace Lh(Th, DataFE<Mesh>::P1); //dc are the discontinuous elements
    Fun_h levelSet(Lh, fun_levelSet);
    Interface2 interface(Th, levelSet.v);

    // Space for the velocity field ! THIS IS FROM BULK SURFACE,
    // ASK HOW VEL.EXPRESSION WORKS, AND WHAT SPACE TO USE
    Lagrange2 FEvelocity(1);
    FESpace2 VelVh(Th, FEvelocity);
    Fun2_h vel(VelVh, fun_velocity);

    // CONSTRUCTION OF THE FE SPACE AND THE CUT SPACE
    FESpace Vh(Th, DataFE<Mesh>::P1dc);
    CutFESpace2  Wh(Vh, interface, {1,-1});
    Mesh2 cutTh(interface);
    FESpace Sh(cutTh, DataFE<Mesh>::P1dc);
    Sh.backSpace = &Vh;

//    gnuplot::save(Th, "../../outputFiles/statConvectionDiffusion/Th.dat");
//    gnuplot::save(cutTh, "../../outputFiles/statConvectionDiffusion/cutTh.dat");
//    gnuplot::save(interface, "../../outputFiles/statConvectionDiffusion/interface.dat");
//    gnuplot::save(Wh,0,"../../outputFiles/statConvectionDiffusion/Vh1.dat");
//    gnuplot::save(Wh,1,"../../outputFiles/statConvectionDiffusion/Vh2.dat");


    // OBJECTS NEEDED FOR THE PROBLEM
    // =====================================================
    CutFEM<Mesh> convdiff(Wh);
    convdiff.add(Sh);
    Normal n;
    Tangent t;


    // variable are [v1,v2,v0]
    int idx_s0 = Wh.NbDoF(); // where v0 start in the solution array
    int idx_bulk = Sh.NbDoF();

    std::cout << "idx_s0 = " << idx_s0 << std::endl;
    std::cout << "idx_bulk = " << idx_bulk << std::endl;

    KN_<double> data_uh(convdiff.rhs(SubArray(Wh.NbDoF(),0)));
    KN_<double> data_sh(convdiff.rhs(SubArray(Sh.NbDoF(),idx_s0)));
    Fun_h uh(Wh, data_uh);  // fun representing bulk variables  [v1, v2]
    Fun_h us(Sh, data_sh);  // fun representing surface variable v0
    // Note that all this doesnt copy any data (so doesnt require time nor memory).
    // It just creates arrays on our existing data rhs. So uh and us will be modify
    // in the same time as rhs.
    // We basically just told the code which part of the array correspond
    // to which variable.



    Fun_h fh (Wh, fun_rhsBulk); // create a FE-function on the cutSpace [v1,v2]
    Fun_h fh0(Sh, fun_rhs0);  // create a function on the cutMesh [v0]
    Fun_h g (Wh, fun_uBulk);    // create an FE-function of the exact bulk solution

    //Fun_h betah(Wh, fun_beta);
    FunTest u(Wh,1), v(Wh,1); // Omega (both subdomains)
    FunTest u0(Sh,1), v0(Sh,1);   // Omega 0, i.e on the interface.
    FunTest u1(Wh,1,0,0), v1(Wh,1,0,0); // Omega1
    FunTest u2(Wh,1,0,1), v2(Wh,1,0,1); // Omega2

    // Assembly of the linear system

    // Integral on element for bulk variables (K_{h,1} and K_{h,2})
    convdiff.addBilinear(
            innerProduct(kappaTildeA*grad(u),grad(v))   // (3.12)
            + innerProduct(kappaTilde*(vel.expression(),grad(u)), v)*0.5 // (3.14)
            - innerProduct(kappaTilde*u,(vel.expression(),grad(v)))*0.5 // (3.14)
    );

    // integral on Edges for bulk variables (i.e. E_{h,1} and E_{h,2})
    convdiff.addBilinear(
            - innerProduct(kappaTildeA*average(grad(u)*n),jump(v))      // (3.12)
            - innerProduct(kappaTildeA*jump(u), average(grad(v)*n))         // (3.13)
            + innerProduct(average((vel*n)*u), kappaTilde*jump(v))*0.5         // (3.15)
            //- innerProduct(jump((vel*n)*u),    kappaTilde*average(v))*0.5      // (3.15)
              - innerProduct(jump(u),    kappaTilde*average((vel*n)*v))*0.5      // (3.15)
            , innerEdge
    );

    // CHANGED FROM LAST VERSION
    // Penalty terms on E_{h,1} and E_{h,2}
    convdiff.addBilinear(
            innerProduct(penalty1*jump(u),jump(v))                  // (3.13)
            + innerProduct(penalty2*jump(u*fabs(vel*n)),jump(v))    // (3.16)
            , innerEdge
    );

    // Integral on interface for surface variable (K_{h,0})
    convdiff.addBilinear(
            innerProduct(A0*gradS(u0),gradS(v0))            // (3.12)
            + innerProduct((vel.expression(),gradS(u0)),v0)*0.5     // (3.14)
            - innerProduct(u0,(vel.expression(),gradS(v0)))*0.5     // (3.14)
            , interface
    );

    // Point evaluation on E_{h,0}
    // TODO:
    //  Figure out which average operators to use

    // This one works best
//    convdiff.addBilinear(
//            - innerProduct(A0*average(gradS(u0)*t, 0.5,-0.5),jump(v0)) // (3.12)
//            - innerProduct(A0*jump(u0), average(gradS(v0)*t, 0.5,-0.5))    // (3.13)
//            + innerProduct(tau_a0*A0/h*jump(u0),jump(v0))                             // (3.13)
//            + innerProduct(average((vel*t)*u0, 0.5,-0.5), jump(v0))*0.5       // (3.15)
//            - innerProduct(jump((vel*t)*u0),    average(v0, 0.5,-0.5))*0.5    // (3.15)
//            + innerProduct(jump(u0*tau_b0*fabs(vel*t)),jump(v0))                  // (3.16)
//            , interface
//            , nodeEvaluation
//    );

    // NOTE: When jump and average operators are used for expressions involving the normal, we have to manually
    // take care of it by using jump = average(func, 1, 1) and average = average(func, 0.5, -0.5).
    convdiff.addBilinear(
            - innerProduct(average(A0*gradS(u0)*t, 0.5,-0.5),jump(v0)) // (3.12)/(3.81)
            - innerProduct(jump(u0), average(A0*gradS(v0)*t, 0.5,-0.5))    // (3.13)/(3.82) added for symmetry
            + innerProduct(tau_a0*A0/h*jump(u0),jump(v0))                             // (3.13)/(3.82)
            + innerProduct(average((vel*t)*u0, 0.5,-0.5), jump(v0))*0.5       // (3.15)/(3.84)
            - innerProduct(jump(u0),    average((vel*t)*v0, 0.5,-0.5))*0.5    // (3.15)/(3.84)
            // + innerProduct(average(tau_b0/fabs(vel*t)*(vel*t)*u0, 1,1),
            //                average((vel*t)*v0, 1, 1))                  // (3.16)
            + innerProduct(tau_b0*jump(u0), jump(v0))
            , interface
            , nodeEvaluation
    );

    // This seems good!
    // Mixed terms
    convdiff.addBilinear(
            + innerProduct(jump(kappa1*u1,kappa01*u0), jump(kappa1*v1, kappa01*v0))*(1.0/kappa01)
            + innerProduct(jump(kappa2*u2,kappa02*u0), jump(kappa2*v2, kappa02*v0))*(1.0/kappa02)
            , interface
    );


    // Stabilization of the bulk
    // Here you will need to get the macro elements and then
    // integrate. You don't have to implement anything here,
    // all is already somewhere. But just try to understand
    // what are those macro elements. Then we could see with
    // each other what to write

    MacroElement macro(Wh, 0.125);
    gnuplot::save(Th);
    gnuplot::save(interface);
    gnuplot::save(macro);

    convdiff.addFaceStabilization(
            innerProduct(1./h*tau_i0*jump(u)      , kappaTilde*jump(v))
            + innerProduct(h*tau_i1*jump(grad(u)), kappaTilde*jump(grad(v)))
            , macro
    );


    // Stabilization of the SURFACE

    // FIXME: This one does not yet work.
    // (3.22)
    MacroElementSurface macroInterface(interface, 0.5);
    gnuplot::save(macroInterface);

    // Full stabilization (3.22)
    convdiff.addBilinear(
            innerProduct(tau00*1./h/h*jump(u0), jump(v0))
            + innerProduct(tau01*jump(grad(u0)), jump(grad(v0)))
            , innerEdge
            // , macroInterface
    );

    // (3.22)
    convdiff.addBilinear(
            innerProduct(tau02*h*h*gradS(u0)*n, gradS(v0)*n)
            , interface
            // , macroInterface
    );

    //// ----------- Dirichlet boundary conditions ---------- ////

    // Nietsche terms to enforce BCs weakly (?)
    convdiff.addBilinear(
            - innerProduct(A1*grad(u1)*n,kappaTilde1*v1)
            - innerProduct(kappaTilde1*u1, A1*grad(v1)*n)
            + innerProduct(kappaTilde1*A1*tau_a1/h*u1, v1)
            , boundary
    );

    // RHS on outer boundary
    convdiff.addLinear(
            - innerProduct(g.expression(), (kappaTilde1*A1)*grad(v1)*n)
            + innerProduct(g.expression(), (kappaTilde1*A1*tau_a1/h)*v1)
            - innerProduct(g.expression(), (vel*n)*kappaTilde1*v1)*0.5,
            boundary
    );

    //// ----------------------------------------------------------------- ////

    // Add linear part on faces
    FunTest fhh0(fh0, 0, 1);

    convdiff.addLinear(
            // innerProduct(fh0.expression(),v0),
            innerProduct(fhh0,v0),
            interface
    );
    // Add RHS on bulk
    convdiff.addLinear(
            innerProduct(fh.expression(),kappaTilde*v)
    );

    //matlab::Export(convdiff.mat, "../../outputFiles/statConvectionDiffusion/matA_n_"+to_string(nx)+".dat");

    // Solve linear system
    convdiff.solve();

    // Compute L2 error
    R errU  = L2normCut(uh, fun_uBulk   , 0,1);
    R errU0 = L2norm(us, fun_uSurface, 0,1);
    //
    std::cout << " || u  - u_ex ||_2  = \t" << errU  << std::endl;
    std::cout << " || u0 - u0_ex||_2  = \t" << errU0 << std::endl;



    // PRINT THE SOLUTION TO PARAVIEW
    // =====================================================
    if(MPIcf::IamMaster()){
        Fun_h uex (Wh, fun_uBulk);
        Paraview2 writer(Wh, levelSet, "convDiff_Bulk.vtk");
        writer.add(uh, "convDiff", 0, 1);
        writer.add(uex,"convDiff_ex", 0, 1);

        Fun_h uexsurf (Sh, fun_uSurface);
        Paraview2 writerU0(Sh, "convDiff_Surface.vtk");
        writerU0.add(us, "convDiff", 0, 1);
        writerU0.add(levelSet, "levelSet", 0, 1);
        writerU0.add(uexsurf, "convDiff_ex_surface", 0, 1);
    }
}
