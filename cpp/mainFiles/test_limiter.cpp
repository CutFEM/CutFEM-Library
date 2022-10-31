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
#include "../FESpace/limiter.hpp"

#define LINEARIZED_EULER

#ifdef CG_TEST
R fun_solution(const R2 P, int elementComp) {
  if(-P.x > 0) return 1;
  else return 20;
}

int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE MESH
  // =====================================================
  int nx = 10;
  int ny = 10;
  Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  FESpace2 Lh(Th, DataFE<Mesh2>::P0);
  FESpace2 Vh(Th, DataFE<Mesh2>::P0);
  LagrangeDC2 peiFE_Flux(0);
  FESpace2 Fh(Th, peiFE_Flux);

  Fun2_h Un(Vh, fun_solution);
  Expression2   expr_Un(Un, 0, op_id);
  Fun2_h flux_Un(Fh, expr_Un, expr_Un);

  Limiter limiter;
  limiter.KXRCF_indicator(Un, flux_Un);
  Fun2_h indicator(Lh, limiter.indicator);

  Rn alphaK(Un.size());
  limiter.limiting(Un, alphaK);
  Fun2_h alpha(Vh, alphaK);

  if(MPIcf::IamMaster() ) {
    Paraview2 writer(Vh,  "limiter.vtk");
    writer.add(Un       , "uh", 0, 1);
    writer.add(indicator, "Ih", 0, 1);
    writer.add(alpha, "alpha_i", 0, 1);
  }

}
#endif


#ifdef LINEARIZED_EULER
double frequency = 1500;
R fun_boundary_left(const R2 P, int elementComp, double t) {
  return 1. + 0.5*sin(2.*pi*frequency*t);
}
// R fun_solution(const R2 P, int elementComp) {
//   if(-P.x > 0) return 1;
//   else return 20;
// }

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;
  typedef std::map<std::pair<int,int>,R> MatMap;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE MESH
  // =====================================================
  int nx = 100;
  int ny = 10;
  Mesh2 Th(nx, ny, 0., 0.05, 1., 0.1);

  // DEFINITION OF THE PARAMETERS
  double dt = 0.000001;
  int ifig = 0;
  double tid = 0.;
  double Tend = 0.0015;
  int niteration = 10;//Tend / dt;
  dt = Tend / niteration;
  double c = 340;
  double rho_0 = 1.224;
  double lambdaE = 1.;


  // DEFINITION OF THE SPACES
  // 3 variables [P, u, v]
  KN<const GTypeOfFE<Mesh2>* > arrayFE(3);
  arrayFE(0) = &DataFE<Mesh2>::P1dc;
  arrayFE(1) = &DataFE<Mesh2>::P1dc;
  arrayFE(2) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> eulerFE(arrayFE);
  FESpace2 Vh(Th, eulerFE);

  // for the flux
  LagrangeDC2 peiFE_Flux(1);
  FESpace2 Fh(Th, peiFE_Flux);

  // For the indicator
  FESpace2 Lh(Th, DataFE<Mesh2>::P0);


  FEM<Mesh2> euler(Vh);
  Normal N;
  MatMap Mh, Ah;


  Rn data_Un(Vh.NbDoF(), 0.);
  Fun2_h fun_Un(Vh, data_Un);

  // U = [rho, rho*u, rho*v, E]
  FunTest U(Vh,3), V(Vh,3);
  FunTest P(Vh,1,0), u(Vh,1,1), v(Vh,1,2);
  FunTest uv(Vh,2,1);


  //COMPUTE THE MASS MATRIX
  euler.setMatrixTo(Ah);
  euler.addBilinear(       // (F(U), grad(V))
            innerProduct(rho_0*c*c*uv, grad(P))
          + innerProduct(1./rho_0*P  , dx(u))
          + innerProduct(1./rho_0*P  , dy(v))
         );

  euler.addBilinear(
    - innerProduct(average(rho_0*c*c*uv*N) , jump(P))
    - innerProduct(average(1./rho_0*P*N.x) , jump(u))
    - innerProduct(average(1./rho_0*P*N.y) , jump(v))
    - innerProduct(0.5*lambdaE*c*jump(U)  , jump(V))

    // - innerProduct(0.5*lambdaE*c*jump(P)           , jump(P))
    // - innerProduct(0.5*lambdaE*c*jump(uv*N)*N.x    , jump(u))
    // - innerProduct(0.5*lambdaE*c*jump(uv*N)*N.y    , jump(v))
    , innerEdge
  );

   // [(rho_0*c*c*u, rho_0*c*c*v),
  // (P/rh0_0, 0 )
  // ( 0, P/rho_0)]
  euler.addBilinear(
    - innerProduct(rho_0*c*c*uv*N, 0.5*P)
    // - innerProduct(1./rho_0*P*N.x, 0.5*u)
    // - innerProduct(1./rho_0*P*N.y, 0.5*v)
    - innerProduct(u, c*0.5*u)
    - innerProduct(v, c*0.5*v)
    , boundary
    , {4}          // label left boundary
  );
  // euler.addBilinear(
  //   - innerProduct(rho_0*c*c*uv*N, 0.5*P)
  //   // - innerProduct(1./rho_0*P*N.x, 0.5*u)
  //   // - innerProduct(1./rho_0*P*N.y, 0.5*v)
  //   // - innerProduct(rho_0*u - P, c*u)
  //   , boundary
  //   , {2}          // label left boundary
  // );
  // euler.addBilinear(
  //   - innerProduct(rho_0*c*c*uv*N, P)
  //   - innerProduct(1./rho_0*P*N.x, u)
  //   - innerProduct(1./rho_0*P*N.y, v)
  //   , boundary
  //   , {2}          // label left boundary
  // );

  euler.setMatrixTo(Mh);   // to fill Mh
  euler.addBilinear(       // time discretization
           innerProduct(U, V)
         );
  euler.reset();


  for(int iter=0;iter<niteration;++iter) {

    std::cout << " Iteration " << iter << " / " << niteration << std::endl;
    std::cout << " Time      " << tid << " / " << Tend << std::endl;

    Fun2_h fun_b4(Vh, fun_boundary_left, tid);
    Expr expr_Uex(fun_b4, 0,op_id);
    FunTest U_b4(fun_b4,0,1);


    MatriceMap<double> mAh(euler.nDoF, euler.nDoF, Ah);
    mAh.addMatMul(data_Un, euler.rhs);

    euler.addLinear(
      - innerProduct(rho_0*c*c*U_b4, (0.5*P)*N.x)
      + innerProduct(U_b4          , c*0.5*u)
      , boundary
      , {4}          // boundary in
    );

    euler.clearMatrix = false;
    euler.solve(Mh, euler.rhs);
    data_Un += dt * euler.rhs;
    euler.rhs = 0.;



    if(MPIcf::IamMaster() ) {
      Paraview2 writer(Vh, "linearizedEuler"+to_string(ifig++)+".vtk");
      writer.add(fun_Un, "P", 0, 1);
      writer.add(fun_Un, "uh", 1, 1);
      writer.add(fun_Un, "vh", 2, 1);
      writer.add(fun_b4, "b4", 0, 1);

    }
    tid += dt;
  }




  // Fun2_h Un(Vh, fun_solution);
  // Expression2   expr_Un(Un, 0, op_id);
  // Fun2_h flux_Un(Fh, expr_Un, expr_Un);
  //
  //
  //














  // Limiter limiter;
  // limiter.KXRCF_indicator(Un, flux_Un);
  // Fun2_h indicator(Lh, limiter.indicator);
  //
  // Rn alphaK(Un.size());
  // limiter.limiting(Un, alphaK);
  // Fun2_h alpha(Vh, alphaK);
  //
  // if(MPIcf::IamMaster() ) {
  //   Paraview2 writer(Vh,  "limiter.vtk");
  //   writer.add(Un       , "uh", 0, 1);
  //   writer.add(indicator, "Ih", 0, 1);
  //   writer.add(alpha, "alpha_i", 0, 1);
  // }

}




#endif
