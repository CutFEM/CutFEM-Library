#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <omp.h>
#include "../util/cputime.h"
#include "cfmpi.hpp"


#include "baseCutProblem.hpp"
#include "levelSet.hpp"
#include "../num/gnuplot.hpp"

typedef std::map<std::pair<int,int>,R> MatMap;

int x_0 = 0;
int x_L = -1;
double c0 = 0.5;
// R fun_solution(const R2 P, int elementComp, int domain) {
//   double t = 0.;
//   if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
//   else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
// }
// R fun_init(const R2 P, int elementComp, int domain, double t) {
//   if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
//   else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
// }



R fun_levelSet1(const R2 P, const int i) {
  double r = P.norm();
  double dy = 1.;
  return r - dy;
}
R fun_levelSet2(const R2 P, const int i) {
  double r = P.norm();
  double dy = 1.384;
  return - r + dy;
}
R fun_levelSet(const R2 P, const int i) {
  return fun_levelSet1(P,i)*fun_levelSet2(P,i);
}
R fun_solution(const R2 P, int elementComp, int domain) {
  double gamma = 1.4;
  double rho_i = 1.;
  double M_i = 2.25;
  double r_i = 1.;
  double r = P.norm();
  double r_ir = r_i/r;
  double rho = rho_i*pow( (1 + (gamma - 1)/2*M_i*M_i*(1-(r_ir*r_ir))), (1./(gamma-1)));
  if(elementComp == 0){
    return rho;
  }
  double cosTheta = (P.x/r);
  double sinTheta = (P.y/r);
  double a_i = 1.;
  double u = a_i*M_i*r_ir*sinTheta;
  double v = -a_i*M_i*r_ir*cosTheta;
  if(elementComp == 1){
    return rho*u;
  }
  if(elementComp == 2){
    return rho*v;
  }
  if(elementComp == 3){
    return pow(rho, gamma)/(gamma*(gamma-1))+0.5*rho*(u*u+v*v);
  }
}
R fluxExact_X(const R2 P, int elementComp, int domain) {
  double rho = fun_solution(P, 0, domain);
  double m   = fun_solution(P, 1, domain);
  double n   = fun_solution(P, 2, domain);
  double E   = fun_solution(P, 3, domain);
  double gamma = 1.4;
  double p = (gamma - 1)*(E-0.5*(m*m/rho+n*n/rho));
  if(elementComp == 0){
    // return 3*rho;
    return m;
  }
  if(elementComp == 1){
    return m*m/rho + p;
  }
  if(elementComp == 2){
    return m*n/rho;
  }
  else{
    return m/rho*(E + p);
  }
}
R fluxExact_Y(const R2 P, int elementComp, int domain) {
  double rho = fun_solution(P, 0, domain);
  double m   = fun_solution(P, 1, domain);
  double n   = fun_solution(P, 2, domain);
  double E   = fun_solution(P, 3, domain);
  double gamma = 1.4;
  double p = (gamma - 1)*(E-0.5*(m*m/rho+n*n/rho));
  if(elementComp == 0){
    // return 1.*rho;
    return n;
  }
  if(elementComp == 1){
    return m*n/rho;
  }
  if(elementComp == 2){
    return n*n/rho+p ;
  }
  else{
    return n/rho*(E + p);
  }
}
R fluxExact_P(const R2 P, int elementComp, int domain) {
  double rho = fun_solution(P, 0, domain);
  double m   = fun_solution(P, 1, domain);
  double n   = fun_solution(P, 2, domain);
  double E   = fun_solution(P, 3, domain);
  double gamma = 1.4;
  double p = (gamma - 1)*(E-0.5*(m*m/rho+n*n/rho));
  return p;
}


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE BACKGROUND MESH
  // =====================================================
  int nx = 101;
  int ny = 101;
  Mesh2 Th(nx, ny, 0., 0., 1.5, 1.5);
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;
  int niteration = 1;
  int ifig = 0;
  double gamma = 1.4;
  double Cstab = 1e-2;
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",12, 2);
  CutFEM_Parameter lambda("lambda",0, 1.);
  double lambdaB_rho = 3;
  double lambdaB_E = 12;
  double lambdaB_m = 7;
  double lambdaB_n = 7;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  Fun2_h levelSetDown(Lh, fun_levelSet1);
  Interface2 interface(Th, levelSetDown.v, 1);
  Fun2_h levelSetUp(Lh, fun_levelSet2);
  interface.add(levelSetUp.v, 2);

  Rn levelSet_data(levelSetDown.v.size(), 0.);
  for(int i=0;i<levelSet_data.size();++i){
    levelSet_data(i) = min(levelSetDown.v(i), levelSetUp.v(i));
  }
  // levelSet_data = levelSetDown.v * levelSetUp.v;
  Fun2_h levelSet(Lh, levelSet_data);

  // Interface2 interface(Th, levelSet.v);

  // if(MPIcf::IamMaster()) {
  //   Paraview2 writer(Lh, levelSet, "levelSet_supersonic_"+to_string(ifig++)+".vtk");
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.add(levelSetUp, "levelSetUp", 0, 1);
  //   writer.add(levelSetDown, "levelSetDown", 0, 1);
  // }

  // CONSTRUCTION OF THE FE SPACES AND THE CUT SPACE
  KN<const GTypeOfFE<Mesh2>* > arrayFE(4);
  arrayFE(0) = &DataFE<Mesh2>::P1dc;
  arrayFE(1) = &DataFE<Mesh2>::P1dc;
  arrayFE(2) = &DataFE<Mesh2>::P1dc;
  arrayFE(3) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> eulerFE(arrayFE);
  FESpace2 Vh(Th, eulerFE);
  CutFESpace2 Wh(Vh, interface, {1});

  // FOR THE PRESSURE
  FESpace2 Ph(Th, DataFE<Mesh2>::P1dc);
  CutFESpace2 W0h(Ph, interface, {1});

  // FOR THE FLUX
  KN<const GTypeOfFE<Mesh2>* > arrayFE_Flux(2);
  arrayFE_Flux(0) = &DataFE<Mesh2>::P1dc;
  arrayFE_Flux(1) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> eulerFE_Flux(arrayFE_Flux);
  FESpace2 Flh(Th, eulerFE_Flux);
  CutFESpace2 WFlh(Flh, interface, {1});

  std::cout << " dof " << Wh.nbDoF << std::endl;

  // gnuplot::save(Th);
  // gnuplot::save(interface, "gamma.dat");
  // gnuplot::save(Wh,0,"Vh1.dat");

  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh2> euler(Wh);
  Normal N;
  MatMap Mh, Ah;

  // variable are [rho, rho*u, rho*v, E]
  Rn data_Un(Wh.NbDoF(), 0.);
  interpolate(Wh, data_Un, fun_solution);
  Fun2_h Un(Wh, data_Un);
  Fun2_h gh(Wh, fun_solution);
  Expression2 gx_rho(gh, 0, op_id);
  Expression2 gx_m(gh, 1, op_id);
  Expression2 gx_n(gh, 2, op_id);
  Expression2 gx_E(gh, 3, op_id);

  Fun2_h F(Wh, fluxExact_X);
  Fun2_h G(Wh, fluxExact_Y);
  Expr Fex_rho(F, 0, op_id);
  Expr Fex_m(F, 1, op_id);
  Expr Fex_n(F, 2, op_id);
  Expr Fex_E(F, 3, op_id);
  Expr Gex_rho(G, 0, op_id);
  Expr Gex_m(G, 1, op_id);
  Expr Gex_n(G, 2, op_id);
  Expr Gex_E(G, 3, op_id);
  FunTest FluxEx_rho(W0h, Fex_rho, Gex_rho);
  FunTest FluxEx_m  (W0h,Fex_m, Gex_m);
  FunTest FluxEx_n  (W0h,Fex_n, Gex_n);
  FunTest FluxEx_E  (W0h,Fex_E, Gex_E);

  if(MPIcf::IamMaster()) {
    Paraview2 writer(Wh, levelSet, "euler_exact.vtk");
    writer.add(gh, "rho", 0, 1);
    writer.add(gh, "u", 1, 2);
    writer.add(gh, "E", 3, 1);
    writer.add(F, "F_rho", 0, 1);
    writer.add(F, "F_u", 1, 2);
    writer.add(F, "F_E", 3, 1);
    writer.add(G, "G_rho", 0, 1);
    writer.add(G, "G_u", 1, 2);
    writer.add(G, "G_E", 3, 1);
    Paraview2 writer2(Lh, levelSet, "eulerMesh.vtk");
  }

  // DEFINE FUNCTION ON THE DATA
  Expr f_rho(Un, 0, op_id);
  Expr f_m  (Un, 1, op_id);
  Expr f_n  (Un, 2, op_id);
  Expr f_E  (Un, 3, op_id);
  FunTest f_Un(Wh, f_rho, f_m, f_n, f_E);
  FunTest rho_n(W0h, f_rho);
  FunTest m_n(W0h, f_m);
  FunTest n_n(W0h, f_n);
  FunTest E_n(W0h, f_E);

  // U = [rho, rho*u, rho*v, E]
  FunTest U(Wh,4), V(Wh,4);
  FunTest rho(Wh,1,0), m(Wh, 1,1), n(Wh,1,2), E(Wh,1,3);

  //COMPUTE THE MASS MATRIX
  euler.setMatrixTo(Ah);
  euler.addFaceStabilization(
      - innerProduct(jump(U), Cstab*jump(V))
      - innerProduct((h^2)*jump(grad(U)), Cstab*jump(grad(V)))
  );
  euler.setMatrixTo(Mh);   // to fill Mh
  euler.addBilinear(
           innerProduct(U, V)
         );
  euler.addFaceStabilization(
         innerProduct(h*jump(U), Cstab*jump(V))
       + innerProduct((h^3)*jump(grad(U)), Cstab*jump(grad(V)))
         );
  euler.reset();   // unset point Mh

  // Fun2_h ftest(W0h, fun_test); Expr ff(ftest, 0, op_id);
  for(int i=0; i<niteration; ++i) {

    double t0 = MPIcf::Wtime();
    // COMPUTE P and ALPHA
    Fun2_h p(W0h, (gamma - 1)*(f_E - 0.5*(f_m*f_m/f_rho + f_n*f_n/f_rho)));
    Expr f_p(p, 0, op_id);
    Fun2_h alpha(W0h, fabs(f_m/f_rho + f_n/f_rho)+sqrt(gamma*f_p/f_rho));
    Expr f_alpha(alpha, 0, op_id);

    // DEFINE FUNCTION AND EXPRESSION FOR THE FLUX
    Fun2_h flux_rho1(W0h,f_m);                 Expr f_1(flux_rho1, 0, op_id);
    Fun2_h flux_m1(W0h,f_m*f_m/f_rho + f_p);   Expr f_2(flux_m1, 0, op_id);
    Fun2_h flux_n1(W0h,f_n*f_m/f_rho);         Expr f_3(flux_n1, 0, op_id);
    Fun2_h flux_E1(W0h,f_m/f_rho*(f_E + f_p)); Expr f_4(flux_E1, 0, op_id);
    // g
    Fun2_h flux_rho2(W0h,f_n);                 Expr g_1(flux_rho2, 0, op_id);
    Fun2_h flux_m2(W0h,f_n*f_m/f_rho);         Expr g_2(flux_m2, 0, op_id);
    Fun2_h flux_n2(W0h,f_n*f_n/f_rho + f_p);   Expr g_3(flux_n2, 0, op_id);
    Fun2_h flux_E2(W0h,f_n/f_rho*(f_E + f_p)); Expr g_4(flux_E2, 0, op_id);


    FunTest Flux_rho(W0h,f_1, g_1);
    FunTest Flux_m  (W0h,f_2, g_2);
    FunTest Flux_n  (W0h,f_3, g_3);
    FunTest Flux_E  (W0h,f_4, g_4);



    // if( MPIcf::IamMaster() ) {
    //   Fun2_h Err_Frho(W0h,fabs(Gex_rho-g_1));
    //   Fun2_h flux_rho(WFlh, f_1, g_1);
    //   Fun2_h fluxEx_rho(WFlh, Fex_rho, Gex_rho);
    //
    //   Paraview2 writerFrho(W0h, levelSet, "variable_rho.vtk");
    //   writerFrho.add(Err_Frho, "err0r_Frho", 0, 1);
    //   writerFrho.add(flux_rho, "Frho_h", 0, 2);
    //   writerFrho.add(fluxEx_rho, "Frho_Ex", 0, 2);
    //
    //
    //   Fun2_h Err_Fm(W0h,fabs(Gex_m-g_2));
    //   Fun2_h flux_m(WFlh, f_2, g_2);
    //   Fun2_h fluxEx_m(WFlh, Fex_m, Gex_m);
    //   Paraview2 writerFm(W0h, levelSet, "variable_m.vtk");
    //   writerFm.add(Err_Fm, "err0r_Fm", 0, 1);
    //   writerFm.add(flux_m, "Fm_h", 0, 2);
    //   writerFm.add(fluxEx_m, "Fm_Ex", 0, 2);
    //
    //   Fun2_h Err_Fn(W0h,fabs(Gex_n-g_3));
    //   Fun2_h flux_n(WFlh, f_3, g_3);
    //   Fun2_h fluxEx_n(WFlh, Fex_n, Gex_n);
    //   Paraview2 writerFn(W0h, levelSet, "variable_n.vtk");
    //   writerFn.add(Err_Fn, "err0r_Fn", 0, 1);
    //   writerFn.add(flux_n, "Fn_h", 0, 2);
    //   writerFn.add(fluxEx_n, "Fn_Ex", 0, 2);
    //
    //   Fun2_h Err_FE(W0h,fabs(Gex_E-g_4));
    //   Fun2_h flux_E(WFlh, f_4, g_4);
    //   Fun2_h fluxEx_E(WFlh, Fex_E, Gex_E);
    //   Paraview2 writerFE(W0h, levelSet, "variable_E.vtk");
    //   writerFE.add(Err_FE, "err0r_FE", 0, 1);
    //   writerFE.add(flux_E, "FE_h", 0, 2);
    //   writerFE.add(fluxEx_E, "FE_Ex", 0, 2);
    // }
    std::cout << " Time create function Flux => " <<  MPIcf::Wtime()-t0 << std::endl;
    double t1 = MPIcf::Wtime();

    // CONSTRUCT THE RHS
    euler.addLinear(
        innerProduct(Flux_rho  , grad(rho))
      + innerProduct(Flux_m  , grad(m))
      + innerProduct(Flux_n  , grad(n))
      + innerProduct(Flux_E  , grad(E))
    );

    euler.addLinear(
      - innerProduct(average(Flux_rho*N), jump(rho))
      - innerProduct(average(Flux_m*N), jump(m))
      - innerProduct(average(Flux_n*N), jump(n))
      - innerProduct(average(Flux_E*N), jump(E))
      - innerProduct(0.5*lambdaE*jump(f_Un)  , jump(V))
      , innerEdge
    );


    // RHO && E BOUNDARY CONDITION IN RHS
    euler.addLinear(
      - innerProduct(Flux_rho, (0.5*rho)*N)
      - innerProduct(Flux_E    , (0.5*E)*N)
      - innerProduct(rho_n   , lambdaB_rho*0.5*rho)
      - innerProduct(E_n     , lambdaB_E*0.5*E)
      , boundary
      , {4}          // boundary in
    );
    euler.addLinear(
      - innerProduct(Flux_rho*N  , rho )
      - innerProduct(Flux_E*N    , E )
      , boundary
      , {1}        // boundary out
    );


    // M && N BOUNDARY CONDITION IN RHS
    euler.addLinear(
      - innerProduct(Flux_m , (0.5*m)*N)
      - innerProduct(Flux_n , (0.5*n)*N)
      - innerProduct(m_n   , lambdaB_m*0.5*m)
      - innerProduct(n_n   , lambdaB_n*0.5*n)
      , interface
      , {1}          // boundary in
    );
    euler.addLinear(
      - innerProduct(Flux_m , (0.5*m)*N)
      - innerProduct(m_n    , lambdaB_m*0.5*m)
      , boundary
      , {4}          // boundary in
    );
    euler.addLinear(
      - innerProduct(Flux_n , (0.5*n)*N)
      - innerProduct(n_n    , lambdaB_n*0.5*n)
      , boundary
      , {1}          // boundary in
    );
    euler.addLinear(
      - innerProduct(Flux_m*N  , m )
      - innerProduct(Flux_n*N  , n )
      , interface
      , {2}        // boundary out
    );

    MatriceMap<double> mAh(euler.nDoF, euler.nDoF, Ah);
    mAh.addMatMul(data_Un, euler.rhs);

    // RHO && E BOUNDARY CONDITION  IN LHS
    euler.addLinear(
      - innerProduct(FluxEx_rho, (0.5*rho)*N)
      - innerProduct(FluxEx_E  , (0.5*E)*N)
      + innerProduct(gx_rho       , lambdaB_rho*0.5*rho)
      + innerProduct(gx_E         , lambdaB_E*0.5*E)
      , boundary
      , {4}          // boundary in
    );
    // M && N BOUNDARY CONDITION IN RHS
    euler.addLinear(
      - innerProduct(FluxEx_m , (0.5*m)*N)
      - innerProduct(FluxEx_n , (0.5*n)*N)
      + innerProduct(FunTest(W0h,gx_m)   , lambdaB_m*0.5*m)
      + innerProduct(FunTest(W0h,gx_n)   , lambdaB_n*0.5*n)
      , interface
      , {1}          // boundary in
    );
    euler.addLinear(
      - innerProduct(FluxEx_m , (0.5*m)*N)
      + innerProduct(gx_m    , lambdaB_m*0.5*m)
      , boundary
      , {4}          // boundary in
    );
    euler.addLinear(
      - innerProduct(FluxEx_n , (0.5*n)*N)
      + innerProduct(gx_n    , lambdaB_n*0.5*n)
      , boundary
      , {1}          // boundary in
    );


    std::cout << " Time assembly RHS => " << MPIcf::Wtime()-t1 << std::endl;

    euler.clearMatrix = false;
    euler.solve(Mh, euler.rhs);

    data_Un += dt * euler.rhs;
    // data_Un = euler.rhs;
    euler.rhs = 0.;
    if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
      Paraview2 writer(Wh, levelSet, "euler_"+to_string(ifig++)+".vtk");
      writer.add(Un, "rho", 0, 1);
      writer.add(Un, "m", 1, 1);
      writer.add(Un, "n", 2, 1);
      writer.add(Un, "E", 3, 1);
    }

  }


}
