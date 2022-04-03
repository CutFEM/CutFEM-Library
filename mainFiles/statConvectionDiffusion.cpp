#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#  include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "projection.hpp"
#include <typeinfo>
#include "../num/matlab.hpp"
#include "paraview.hpp"





namespace example2D{
  // RHS for variable in Omega 1
  double fun_rhs1(const R2 P) {
    return -exp(- pow(P.x,2) - pow(P.y,2) + 1)*(12*pow(P.x,4)*P.y + 3*pow(P.x,3) + 8*pow(P.x,2)*pow(P.y,3) - 48*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 4*pow(P.y,5) + 16*pow(P.y,3));
  }

  // RHS for variable in Omega 2
  double fun_rhs2(const R2 P) {
    return -2*exp(- pow(P.x,2) - pow(P.y,2) + 1)*(6*pow(P.x,4)*P.y + 3*pow(P.x,3) + 4*pow(P.x,2)*pow(P.y,3) - 24*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 2*pow(P.y,5) + 8*pow(P.y,3));
  }

  // RHS for bulk variables
  double fun_rhsBulk(const R2 P, int elementComp, int dom) {
    if(dom == 0) return fun_rhs1(P);
    else return fun_rhs2(P);
  }

  // RHS for surface variable
  double fun_rhs0(const R2 P, int elementComp) {
    return 3*(- pow(P.x,5) + 2*pow(P.x,3)*pow(P.y,2) + 9*pow(P.x,2)*P.y + 3*P.x*pow(P.y,4) - 3*pow(P.y,3))/(pow(P.x,2) + pow(P.y,2));
  }

  // Exact solution in the bulk
  double fun_uBulk(const R2 P,  int elementComp, int domain) {
    if(domain == 0) return     exp(1.0-P.x*P.x-P.y*P.y)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
    else            return 2.0*exp(1.0-P.x*P.x-P.y*P.y)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
  }

  // Exact solution on surface
  double fun_uSurface(const R2 P,  int elementComp) {
    return 3.0*P.x*P.x*P.y - pow(P.y,3);
  }

  // // Exact solution on surface for time-step t
  // double fun_uSurfaceT(const R2 P,  int elementComp, double t) {
  //   return 3.0*P.x*P.x*P.y - pow(P.y,3);
  // }

  // Level-set function
  double fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x)*(P.x) + (P.y)*(P.y)) - 1.0;
  }

  // Velocity Field
  R fun_velocity(const R2 P, const int i){
    if(i == 0) return P.y;
    else       return -P.x;
  }
}
namespace example3D{
  // RHS for variable in Omega 1
  double fun_rhs1(const R3 P) {
    return -exp(- pow(P.x,2) - pow(P.y,2) + 1)*(12*pow(P.x,4)*P.y + 3*pow(P.x,3) + 8*pow(P.x,2)*pow(P.y,3) - 48*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 4*pow(P.y,5) + 16*pow(P.y,3));
  }

  // RHS for variable in Omega 2
  double fun_rhs2(const R3 P) {
    return -2*exp(- pow(P.x,2) - pow(P.y,2) + 1)*(6*pow(P.x,4)*P.y + 3*pow(P.x,3) + 4*pow(P.x,2)*pow(P.y,3) - 24*pow(P.x,2)*P.y - 9*P.x*pow(P.y,2) - 2*pow(P.y,5) + 8*pow(P.y,3));
  }

  // RHS for bulk variables
  double fun_rhsBulk(const R3 P, int elementComp, int dom) {
    if(dom == 0) return fun_rhs1(P);
    else return fun_rhs2(P);
  }

  // RHS for surface variable
  double fun_rhs0(const R3 P, int elementComp) {
    return 3*(- pow(P.x,5) + 2*pow(P.x,3)*pow(P.y,2) + 9*pow(P.x,2)*P.y + 3*P.x*pow(P.y,4) - 3*pow(P.y,3))/(pow(P.x,2) + pow(P.y,2));
  }

  // Exact solution in the bulk
  double fun_uBulk(const R3 P,  int elementComp, int domain) {
    if(domain == 0) return     exp(1.0-P.x*P.x-P.y*P.y-P.z*P.z)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
    else            return 2.0*exp(1.0-P.x*P.x-P.y*P.y-P.z*P.z)*( 3.0*P.x*P.x*P.y - pow(P.y,3));
  }

  // Exact solution on surface
  double fun_uSurface(const R3 P,  int elementComp) {
    return 3.0*P.x*P.x*P.y - pow(P.y,3);
  }

  // Exact solution on surface for time-step t
  double fun_uSurfaceT(const R3 P,  int elementComp, double t) {
    return 3.0*P.x*P.x*P.y - pow(P.y,3);
  }

  // Level-set function
  double fun_levelSet(const R3 P, const int i) {
    return sqrt((P.x)*(P.x) + (P.y)*(P.y)+ (P.z)*(P.z)) - 1.0;
  }

  // Velocity Field
  R fun_velocity(const R3 P, const int i){
    if(i == 0) return P.z/10;
    else if(i == 1) return 0;
    else       return -P.x/10;
  }
}

struct kappa_E1 : public Virtual_Parameter {
  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    return meas_Cut/measK;
  }
};
struct kappa_E2 : public Virtual_Parameter {
  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    return 1-meas_Cut/measK;
  }
};


// 2 DIMENSIONAL PROBLEM
// =====================================================


void solve(int argc, char** argv, int nn, int i) {
  using namespace example2D;
  typedef Mesh2 Mesh;
  typedef FESpace2 Space;
  typedef CutFESpace<Mesh> CutSpace;

  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

    // 36
    int nx = nn, ny = nn;
    double lx = 3., ly = 3.;
    double h = lx/(nx-1);

    // Problem parameters
    double A0 = 1.0, A1 = 1.0, A2 = 0.5;
    double kappa1 = 2.0, kappa2 = 0.5, kappa01 = 1, kappa02 = 2;
    double kappaTilde1 = kappa1/kappa01, kappaTilde2 = kappa2/kappa02;
    ParameterCutFEM kappaTilde(kappaTilde1    , kappaTilde2);
    ParameterCutFEM kappaTildeA(kappaTilde1*A1, kappaTilde2*A2);

    // FULLSTAB PARAMETERS
    double tau_a0 = 1e0, tau_b0 = 0;
    double tau_a1 = 1e1, tau_b1 = 0;
    double tau_a2 = 1e1, tau_b2 = 0;
    double lambdaB = A1*1e1/h;
    double lambdaA0 = A0*5e3/h;

    // WANT
    //    double lambdaB = A1*1e1/h;    // coefficient on the outer boundary
    //    double lambdaA0 = A0*1e1/h;   // for the surface problem

    // Local bulk penalty parameters
    // ParameterCutFEM A(A1, A2);
    // ParameterCutFEM penalty1(Parameter::lambdaA);

    // Global bulk penalty parameter (if this is commented, the two lines above should be out-commented
    ParameterCutFEM penalty1(kappaTilde1*A1*tau_a1/h, kappaTilde2*A2*tau_a2/h);
    ParameterCutFEM penalty2(kappaTilde1*tau_b1, kappaTilde2*tau_b2);
    double penalty0 = tau_b0*sqrt(2), penalty0div = tau_b0/sqrt(2);
    double lambda1 = 1.0;

    // Stabilization parameters
    double tau11 = 0.1*A1, tau10 = A1, tau20 = A2, tau21 = 0.1*A2;
    double tau00 = A0, tau01 = A0, tau02 = 0.1*A0;
    ParameterCutFEM tau_i0(tau10, tau20);
    ParameterCutFEM tau_i1(tau11, tau21);

    // CONSTRUCTION OF THE MESH
    // =====================================================
    Mesh Kh(nx, ny, -1.5, -1.5, lx, ly);
    Space Vh(Kh, DataFE<Mesh>::P1dc);
    Vh.info();

    // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
    // =====================================================
    Space Lh(Kh, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);

    // Space for the velocity field ! THIS IS FROM BULK SURFACE,
    Lagrange2 FEvelocity(1);
    Space VelVh(Kh, FEvelocity);
    Fun_h vel(VelVh, fun_velocity);

    // CONSTRUCTION OF THE ACTIVE MESH
    ActiveMesh<Mesh> Khi(Kh, interface);
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.create_surface_mesh(interface);
    // Khi.info();
    // Kh0.info();

    CutSpace Wh(Khi, Vh);
    CutSpace Sh(Kh0, Vh);
    Wh.info();
    Sh.info();


    // OBJECTS NEEDED FOR THE PROBLEM
    // =====================================================
    CutFEM<Mesh> convdiff(Wh);
    convdiff.add(Sh);
    Normal n;
    Tangent t;
    Conormal conormal;

    // variable are [v1,v2,v0]
    int idx_s0 = Wh.NbDoF(); // where v0 start in the solution array
    int idx_intf = Sh.NbDoF();

    KN_<double> data_uh(convdiff.rhs_(SubArray(Wh.NbDoF(),0)));
    KN_<double> data_sh(convdiff.rhs_(SubArray(Sh.NbDoF(),idx_s0)));
    Fun_h uh(Wh, data_uh);  // fun representing bulk variables  [v1, v2]
    Fun_h us(Sh, data_sh);  // fun representing surface variable v0

    Fun_h fh (Wh, fun_rhsBulk); // create a FE-function on the cutSpace [v1,v2]
    Fun_h fh0(Sh, fun_rhs0);    // create a function on the cutMesh [v0]
    Fun_h g  (Wh, fun_uBulk);    // create an FE-function of the exact bulk solution

    FunTest u(Wh,1), v(Wh,1); // Omega (both subdomains)
    FunTest u0(Sh,1), v0(Sh,1);   // Omega 0, i.e on the interface.
    FunTest u1(Wh,1,0,0), v1(Wh,1,0,0); // Omega1
    FunTest u2(Wh,1,0,1), v2(Wh,1,0,1); // Omega2

    // Assembly of the linear system
    // Integral on element for bulk variables (K_{h,1} and K_{h,2})
    convdiff.addBilinear(
            innerProduct(kappaTildeA*grad(u),grad(v))             // (3.12)
            + innerProduct(kappaTilde*(vel.expression(),grad(u)), v)*0.5 // (3.14)
            - innerProduct(kappaTilde*u,(vel.expression(),grad(v)))*0.5  // (3.14)
            , Khi
    );
    matlab::Export(convdiff.mat_, "mat1_new.dat");
    convdiff.cleanMatrix();
    // integral on Edges for bulk variables (i.e. E_{h,1} and E_{h,2})
    // GLOBAL WEIGHTS
    convdiff.addBilinear(
            - innerProduct(kappaTildeA*average(grad(u)*n),jump(v))      // (3.12)
            - innerProduct(kappaTildeA*jump(u), average(grad(v)*n))         // (3.13)
            + innerProduct(average((vel*n)*u), kappaTilde*jump(v))*0.5         // (3.15)
            - innerProduct(jump((vel*n)*u),    kappaTilde*average(v))*0.5      // (3.15)
            , Khi
            , innerFacet
    );
    matlab::Export(convdiff.mat_, "mat2_new.dat");
    convdiff.cleanMatrix();
//
//     // LOCAL WEIGHTS
// //    convdiff.addBilinear(
// //            - innerProduct(kappaTildeA*average(grad(u)*n, kappa_E1, kappa_E2),jump(v))      // (3.12)
// //            - innerProduct(kappaTildeA*jump(u), average(grad(v)*n, kappa_E1, kappa_E2))         // (3.13)
// //            + innerProduct(average((vel*n)*u, kappa_E1, kappa_E2), kappaTilde*jump(v))*0.5         // (3.15)
// //            - innerProduct(jump((vel*n)*u),    kappaTilde*average(v, kappa_E1, kappa_E2))*0.5      // (3.15)
// //            , innerEdge
// //    );
//     //getchar();
//     //std::cout << "AFTER" << std::endl;

    // Penalty terms on E_{h,1} and E_{h,2}
    convdiff.addBilinear(
            innerProduct(penalty1*jump(u),jump(v))                  // (3.13)
            + innerProduct(penalty2*jump(u*fabs(vel*n)),jump(v))    // (3.16)
            , Khi
            , innerFacet
    );
    matlab::Export(convdiff.mat_, "mat3_new.dat");
    convdiff.cleanMatrix();

    // Integral on interface for surface variable (K_{h,0})
    convdiff.addBilinear(
            innerProduct(A0*gradS(u0),gradS(v0))             // (3.12)
            + innerProduct((vel.expression(),gradS(u0)),v0)*0.5     // (3.14)
            - innerProduct(u0,(vel.expression(),gradS(v0)))*0.5     // (3.14)
            , interface
    );
    matlab::Export(convdiff.mat_, "mat4_new.dat");
    convdiff.cleanMatrix();
    //// -------- Point evaluation on E_{h,0} ------- //
    // VARIANT 2
    convdiff.addBilinear(
      + innerProduct(average((*(vel.expression(1).begin())*u0),0.5,-0.5), jump((*(vel.expression(1).begin())*v0)))*0.5       // (3.15)/(3.84)
      // - innerProduct(average(A0*gradS(u0)*conormal,0.5,-0.5),jump(v0)) // (3.12)/(3.81)
      // - innerProduct(jump(u0), average(A0*gradS(v0)*conormal,0.5,-0.5))    // (3.13)/(3.82) added for symmetry
      // + innerProduct(lambdaA0*jump(u0),jump(v0))                             // (3.13)/(3.82)
      // + innerProduct(average((vel*conormal)*u0,0.5,-0.5), jump(v0))*0.5       // (3.15)/(3.84)
      // - innerProduct(jump(u0),    average((vel*conormal)*v0,0.5,-0.5))*0.5    // (3.15)/(3.84)
      // + innerProduct(penalty0*jump(u0), jump(v0))                  // (3.16)
      , interface
      , innerRidge
    );
    matlab::Export(convdiff.mat_, "mat5_new.dat");
    convdiff.cleanMatrix();

    // Mixed terms
    convdiff.addBilinear(
            + innerProduct(jump(kappa1*u1,kappa01*u0), jump(kappa1*v1, kappa01*v0))*(1.0/kappa01)
            + innerProduct(jump(kappa2*u2,kappa02*u0), jump(kappa2*v2, kappa02*v0))*(1.0/kappa02)
            , interface
    );
    matlab::Export(convdiff.mat_, "mat6_new.dat");
    convdiff.cleanMatrix();

    //// Stabilization of the bulk
    //MacroElement macro(Wh, 0.125);
    convdiff.addFaceStabilization(
      innerProduct(1./h*tau_i0*jump(u)      , kappaTilde*jump(v))
      + innerProduct(h*tau_i1*jump(grad(u)), kappaTilde*jump(grad(v)))
      , Khi
      //            , macro
    );

    //// Stabilization of the SURFACE
    //MacroElementSurface macroInterface(interface, 0.25);
    convdiff.addBilinear(
      innerProduct(tau00*1./h/h*jump(u0), jump(v0))
      + innerProduct(tau01*jump(grad(u0)), jump(grad(v0)))
      , Kh0
      , innerFacet
      //      , macroInterface
    );

    convdiff.addBilinear(
      innerProduct(tau02*h*h*grad(u0)*n, grad(v0)*n)
      , interface
    );
    matlab::Export(convdiff.mat_, "mat7_new.dat");
    convdiff.cleanMatrix();
//    gnuplot::save(macro,"../../outputFiles/statConvectionDiffusion/gnuplot/");
//    gnuplot::save(macroInterface,"../../outputFiles/statConvectionDiffusion/gnuplot/");

    // ----- Add Dirichlet boundary conditions ------
    convdiff.addBilinear(
      - innerProduct(kappaTilde1*A1*grad(u1)*n, v1)
      - innerProduct(kappaTilde1*A1*u1, grad(v1)*n)
      + innerProduct(kappaTilde1*lambdaB*u1, v1)
      , Khi
      , boundary
    );

    // RHS on outer boundary
    convdiff.addLinear(
      - innerProduct(g.expression(), (kappaTilde1*A1)*grad(v1)*n)
      + innerProduct(g.expression(), (kappaTilde1*lambdaB)*v1)
      - innerProduct(g.expression(), (vel*n)*kappaTilde1*v1)*0.5
      , Khi
      , boundary
    );
    matlab::Export(convdiff.mat_, "mat8_new.dat");
    convdiff.cleanMatrix();
   // Add linear part on faces
   convdiff.addLinear(
     innerProduct(fh0.expression(),v0)
     , interface
   );
   // Add RHS on bulk
   convdiff.addLinear(
     innerProduct(fh.expression(),kappaTilde*v)
     , Khi
   );
   matlab::Export(convdiff.rhs_, "rhs_new.dat");
   convdiff.cleanMatrix();
   // matlab::Export(convdiff.mat_, "matCONV.dat");


   // Solve linear system
   convdiff.solve();

   // Compute L2 error
   R errU  = L2normCut(uh, fun_uBulk   , 0,1);
   R errU0 = L2normSurf(us, fun_uSurface, interface, 0,1);

   // std::cout << convdiff.rhs_ << std::endl;

   std::cout << " || u  - u_ex ||_2  = \t" << errU  << std::endl;
   std::cout << " || u0 - u0_ex||_2  = \t" << errU0 << std::endl;

    // PRINT THE SOLUTION TO PARAVIEW
    // =====================================================
    if(MPIcf::IamMaster()){
        Fun_h uex (Wh, fun_uBulk);
        Paraview<Mesh> writer(Khi, "convDiff_Bulk"+to_string(i)+".vtk");
        writer.add(uh, "convDiff", 0, 1);
        writer.add(uex,"convDiff_ex", 0, 1);

        Fun_h uexsurf(Sh, fun_uSurface);
        Paraview<Mesh> writerU0(Kh0, "convDiff_Surface"+to_string(i)+".vtk");
        writerU0.add(us, "convDiff", 0, 1);
        writerU0.add(levelSet, "levelSet", 0, 1);
        writerU0.add(uexsurf, "convDiff_ex_surface", 0, 1);
    }

}

//
// void solve3D(int argc, char** argv) {
//   using namespace example3D;
//   typedef Mesh3 Mesh;
//   typedef FESpace3 FESpace;
//   typedef TestFunction<3> FunTest;
//   typedef FunFEM<Mesh3> Fun_h;
//     // INITIALIZE MPI
//     // =====================================================
//     MPIcf cfMPI(argc,argv);
//
//     // 36
//     int nx = 5, ny = 5, nz = 5;
//     double lx = 3., ly = 3.,lz = 3.;
//     double h = lx/(nx-1);
//     //double h = sqrt(2*lx/(nx-1)*lx/(nx-1));     // hypotenuse of triangle (all triangles are the same)
//     std::cout << "h = " << h << std::endl;
//
//     // Problem parameters
//     double A0 = 1.0, A1 = 1.0, A2 = 0.5;
//     double kappa1 = 2.0, kappa2 = 0.5, kappa01 = 1, kappa02 = 2;
//     double kappaTilde1 = kappa1/kappa01, kappaTilde2 = kappa2/kappa02;
//     ParameterCutFEM kappaTilde("kappaTilde", kappaTilde1, kappaTilde2);
//     ParameterCutFEM kappaTildeA("kappaTildeA", kappaTilde1*A1, kappaTilde2*A2);
//
//     // Local weights for average function across bulk faces
//     ParameterCutFEM kappa_E1("kappa_E1", fun_kappa_E1);
//     ParameterCutFEM kappa_E2("kappa_E2", fun_kappa_E2);
//
//     // Penalty parameters, arbitrarily chosen at this time.
//
//     // FULLSTAB PARAMETERS
//     double tau_a0 = 1e0, tau_b0 = 0;
//     double tau_a1 = 1e1, tau_b1 = 0;
//     double tau_a2 = 1e1, tau_b2 = 0;
// //    double lambdaB = A1*2e2/h;    // coefficient on the outer boundary
// //    double lambdaA0 = A0*7e3/h;   // for the surface problem
//     double lambdaB = A1*1e1/h;    // coefficient on the outer boundary
//     //double lambdaA0 = A0*5e1/h;   // for the surface problem
//     double lambdaA0 = A0*5e3/h;
//
//     // WANT
// //    double lambdaB = A1*1e1/h;    // coefficient on the outer boundary
// //    double lambdaA0 = A0*1e1/h;   // for the surface problem
//
//     // Local bulk penalty parameters
//     //ParameterCutFEM A("A", A1, A2);
//     //ParameterCutFEM penalty1(Parameter::lambdaA);
//
//     // Global bulk penalty parameter (if this is commented, the two lines above should be out-commented
//     ParameterCutFEM penalty1("penalty1", kappaTilde1*A1*tau_a1/h, kappaTilde2*A2*tau_a2/h);
//
//     // Other global penalty parameter (do not change)
//     ParameterCutFEM penalty2("penalty2", kappaTilde1*tau_b1, kappaTilde2*tau_b2);
//     double penalty0 = tau_b0*sqrt(2), penalty0div = tau_b0/sqrt(2);
//     double lambda1 = 1.0;
//
//     // Stabilization parameters
//     double tau11 = 0.1*A1, tau10 = A1, tau20 = A2, tau21 = 0.1*A2;
//     double tau00 = A0, tau01 = A0, tau02 = 0.1*A0;
//     //double tau00 = 500000*A0, tau01 = 1e-2*A0, tau02 = 0.01*A0;
//     ParameterCutFEM tau_i0("taui0", tau10, tau20);
//     ParameterCutFEM tau_i1("taui1", tau11, tau21);
//
//     // CONSTRUCTION OF THE MESH
//     // =====================================================
//     Mesh Kh(nx, ny, nz, -1.5, -1.5, -1.5, lx, ly, lz);
//     Kh.info();
//     // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
//     // =====================================================
//     FESpace Lh(Kh, DataFE<Mesh>::P1); //dc are the discontinuous elements
//     Lh.info();
//
//     Fun_h levelSet(Lh, fun_levelSet);
//     Interface3 interface(Kh, levelSet.v);
//
//     // Space for the velocity field ! THIS IS FROM BULK SURFACE,
//     // ASK HOW VEL.EXPRESSION WORKS, AND WHAT SPACE TO USE
//     Lagrange3 FEvelocity(1);
//     FESpace VelVh(Kh, FEvelocity);
//     Fun_h vel(VelVh, fun_velocity);
//
//     // CONSTRUCTION OF THE FE SPACE AND THE CUT SPACE
//     FESpace Vh(Kh, DataFE<Mesh>::P1dc);
//     Vh.info();
//
//     CutFESpace3 Wh(Vh, interface, {1,-1});
//
//     Mesh cutTh(interface);
//     FESpace Sh(cutTh, interface, DataFE<Mesh>::P1dc);
//     Sh.info();
//     Sh.backSpace = &Vh;
//
//     // OBJECTS NEEDED FOR THE PROBLEM
//     // =====================================================
//     CutFEM<Mesh> convdiff(Wh);
//     convdiff.add(Sh);
//     Normal n;
//
//     // variable are [v1,v2,v0]
//     int idx_s0 = Wh.NbDoF(); // where v0 start in the solution array
//     int idx_intf = Sh.NbDoF();
//
//
//     KN_<double> data_uh(convdiff.rhs(SubArray(Wh.NbDoF(),0)));
//     KN_<double> data_sh(convdiff.rhs(SubArray(Sh.NbDoF(),idx_s0)));
//     Fun_h uh(Wh, data_uh);  // fun representing bulk variables  [v1, v2]
//     Fun_h us(Sh, data_sh);  // fun representing surface variable v0
//
//     Fun_h fh (Wh, fun_rhsBulk); // create a FE-function on the cutSpace [v1,v2]
//     Fun_h fh0(Sh, fun_rhs0);  // create a function on the cutMesh [v0]
//     Fun_h g (Wh, fun_uBulk);    // create an FE-function of the exact bulk solution
//
//     FunTest u(Wh,1), v(Wh,1); // Omega (both subdomains)
//     FunTest u0(Sh,1), v0(Sh,1);   // Omega 0, i.e on the interface.
//     FunTest u1(Wh,1,0,0), v1(Wh,1,0,0); // Omega1
//     FunTest u2(Wh,1,0,1), v2(Wh,1,0,1); // Omega2
//
//     // Assembly of the linear system
//
//     // Integral on element for bulk variables (K_{h,1} and K_{h,2})
//     convdiff.addBilinear(
//             innerProduct(kappaTildeA*grad(u),grad(v))             // (3.12)
//             + innerProduct(kappaTilde*(vel.expression(),grad(u)), v)*0.5 // (3.14)
//             - innerProduct(kappaTilde*u,(vel.expression(),grad(v)))*0.5  // (3.14)
//     );
//
//     // integral on Edges for bulk variables (i.e. E_{h,1} and E_{h,2})
//      // GLOBAL WEIGHTS
//      convdiff.addBilinear(
//              - innerProduct(kappaTildeA*average(grad(u)*n),jump(v))      // (3.12)
//              - innerProduct(kappaTildeA*jump(u), average(grad(v)*n))         // (3.13)
//              + innerProduct(average((vel*n)*u), kappaTilde*jump(v))*0.5         // (3.15)
//              - innerProduct(jump((vel*n)*u),    kappaTilde*average(v))*0.5      // (3.15)
//              , innerFace
//      );
// //
// //     // LOCAL WEIGHTS
// // //    convdiff.addBilinear(
// // //            - innerProduct(kappaTildeA*average(grad(u)*n, kappa_E1, kappa_E2),jump(v))      // (3.12)
// // //            - innerProduct(kappaTildeA*jump(u), average(grad(v)*n, kappa_E1, kappa_E2))         // (3.13)
// // //            + innerProduct(average((vel*n)*u, kappa_E1, kappa_E2), kappaTilde*jump(v))*0.5         // (3.15)
// // //            - innerProduct(jump((vel*n)*u),    kappaTilde*average(v, kappa_E1, kappa_E2))*0.5      // (3.15)
// // //            , innerEdge
// // //    );
// //     //getchar();
// //     //std::cout << "AFTER" << std::endl;
// //
//     // Penalty terms on E_{h,1} and E_{h,2}
//     // convdiff.addBilinear(
//     //     innerProduct(penalty1*jump(u),jump(v))                  // (3.13)
//     //   + innerProduct(penalty2*jump(u*fabs(vel*n)),jump(v))    // (3.16)
//     //   , innerFace
//     // );
//     //
//     // // Integral on interface for surface variable (K_{h,0})
//     // convdiff.addBilinear(
//     //         innerProduct(A0*gradS(u0),gradS(v0))             // (3.12)
//     //         + innerProduct((vel.expression(),gradS(u0)),v0)*0.5     // (3.14)
//     //         - innerProduct(u0,(vel.expression(),gradS(v0)))*0.5     // (3.14)
//     //         , interface
//     // );
// //
// //     //// -------- Point evaluation on E_{h,0} ------- //
// //     // VARIANT 2
// //     convdiff.addBilinear(
// //             - innerProduct(average(A0*gradS(u0)*t,0.5,-0.5),jump(v0)) // (3.12)/(3.81)
// //             - innerProduct(jump(u0), average(A0*gradS(v0)*t,0.5,-0.5))    // (3.13)/(3.82) added for symmetry
// //             //+ innerProduct(tau_a0*A0/h*jump(u0),jump(v0))                             // (3.13)/(3.82)
// //             + innerProduct(lambdaA0*jump(u0),jump(v0))                             // (3.13)/(3.82)
// //             + innerProduct(average((vel*t)*u0,0.5,-0.5), jump(v0))*0.5       // (3.15)/(3.84)
// //             - innerProduct(jump(u0),    average((vel*t)*v0,0.5,-0.5))*0.5    // (3.15)/(3.84)
// //             + innerProduct(penalty0*jump(u0),
// //                            jump(v0))                  // (3.16)
// //             , interface
// //            , innerEdge
// //     );
// //
// //
//     // Mixed terms
//     convdiff.addBilinear(
//             + innerProduct(jump(kappa1*u1,kappa01*u0), jump(kappa1*v1, kappa01*v0))*(1.0/kappa01)
//             + innerProduct(jump(kappa2*u2,kappa02*u0), jump(kappa2*v2, kappa02*v0))*(1.0/kappa02)
//             , interface
//     );
// //
// //
// //     //// Stabilization of the bulk
// //     //MacroElement macro(Wh, 0.125);
// //     convdiff.addFaceStabilization(
// //             innerProduct(1./h*tau_i0*jump(u)      , kappaTilde*jump(v))
// //             + innerProduct(h*tau_i1*jump(grad(u)), kappaTilde*jump(grad(v)))
// //       //            , macro
// //     );
// //
// //     //// Stabilization of the SURFACE
// //     //MacroElementSurface macroInterface(interface, 0.25);
// //     convdiff.addBilinear(
// //             innerProduct(tau00*1./h/h*jump(u0), jump(v0))
// //             + innerProduct(tau01*jump(grad(u0)), jump(grad(v0)))
// //             , innerEdge
// //       //      , macroInterface
// //     );
// //
// //     convdiff.addBilinear(
// //             innerProduct(tau02*h*h*grad(u0)*n, grad(v0)*n)
// //             , interface
// //     );
// //
// // //    gnuplot::save(macro,"../../outputFiles/statConvectionDiffusion/gnuplot/");
// // //    gnuplot::save(macroInterface,"../../outputFiles/statConvectionDiffusion/gnuplot/");
// //
// //     // ----- Add Dirichlet boundary conditions ------
// //     convdiff.addBilinear(
// //             - innerProduct(kappaTilde1*A1*grad(u1)*n, v1)
// //             - innerProduct(kappaTilde1*A1*u1, grad(v1)*n)
// //             + innerProduct(kappaTilde1*lambdaB*u1, v1)
// //             , boundary
// //     );
// //
// //     // RHS on outer boundary
// //     convdiff.addLinear(
// //             - innerProduct(g.expression(), (kappaTilde1*A1)*grad(v1)*n)
// //             + innerProduct(g.expression(), (kappaTilde1*lambdaB)*v1)
// //             - innerProduct(g.expression(), (vel*n)*kappaTilde1*v1)*0.5,
// //             boundary
// //    );
// //
// //     // Add linear part on faces
// //     convdiff.addLinear(
// //             innerProduct(fh0.expression(),v0),
// //             interface
// //     );
// //     // Add RHS on bulk
// //     convdiff.addLinear(
// //             innerProduct(fh.expression(),kappaTilde*v)
// //     );
// //     std::cout << "h = " << h << std::endl;
// //
// //     // Solve linear system
// //     convdiff.solve();
//
//     // Compute L2 error
//     // R errU  = L2normCut(uh, fun_uBulk   , 0,1);
//     // R errU0 = L2normSurf(us, fun_uSurfaceT, 0,0,1);
//     // //
//     // std::cout << " || u  - u_ex ||_2  = \t" << errU  << std::endl;
//     // std::cout << " || u0 - u0_ex||_2  = \t" << errU0 << std::endl;
//
//     // PRINT THE SOLUTION TO PARAVIEW
//     // =====================================================
//     if(MPIcf::IamMaster()){
//         Fun_h uex (Wh, fun_uBulk);
//         Paraview<Mesh3> writer(Wh, levelSet, "convDiff_Bulk3.vtk");
//         // writer.add(uh, "convDiff", 0, 1);
//         // writer.add(uex,"convDiff_ex", 0, 1);
//
//         // Fun_h uexsurf(Sh, fun_uSurface);
//         Paraview<Mesh3> writerU0(Sh, "convDiff_Surface3.vtk");
//         // // writerU0.add(us, "convDiff", 0, 1);
//         writerU0.add(levelSet, "levelSet", 0, 1);
//         // // writerU0.add(uexsurf, "convDiff_ex_surface", 0, 1);
//     }
//
// }

int main(int argc, char** argv ) {

  // INITIALIZE MPI
  // =====================================================
  MPIcf cfMPI(argc,argv);
  int nx = 11;
  for(int i=0; i<1; ++i){

    solve(argc, argv, nx, i);

    nx = 2*nx-1;
  }
}
