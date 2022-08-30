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
#include "../FESpace/limiter.hpp"
#include "../problem/projection.hpp"
#include "../num/gnuplot.hpp"


#define EXAMPLE_FEM
typedef std::map<std::pair<int,int>,R> MatMap;

#ifdef EXAMPLE_FEM

R fun_init_solution(const R2 P, int elementComp) {
  return (P.x <= 0);
}
R fun_solution(const R2 P, int elementComp, double t) {
  double x_shock = 0.5*t;
  return (P.x <= x_shock);
}


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  int nx = 51;
  int ny = 51;
  Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;
  int ifig = 0;
  double tid = 0.;
  double Tend = 0.1;
  int niteration = Tend / dt;
  dt = Tend / niteration;


  // DEFINITION OF PARAMETERS
  // =====================================================
  double Cstab = 1e-2;
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",1, 1);
  double lambdaB = 1.;


  // CONSTRUCTION OF THE FE SPACES AND THE CUT SPACE FOR Uh AND FLUX
  // =====================================================
  // For the exact solution
  FESpace2 Eh(Th, DataFE<Mesh2>::P1dc);

  // For the computation
  FESpace2 Vh(Th, DataFE<Mesh2>::P1dcTaylor);

  // For the flux
  KN<const GTypeOfFE<Mesh2>* > arrayFE_Flux(2);
  arrayFE_Flux(0) = &DataFE<Mesh2>::P1dcTaylor;
  arrayFE_Flux(1) = &DataFE<Mesh2>::P1dcTaylor;
  GTypeOfFESum<Mesh2> burgerFE_Flux(arrayFE_Flux);
  FESpace2 Fh(Th, burgerFE_Flux);

  // For the indicator
  FESpace2 Ih(Th, DataFE<Mesh2>::P0);


  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  FEM<Mesh2> burger(Vh);
  Normal N;
  MatMap Mh, Ah;

  // FUNCTION Un AND TESTFUNCTION
  // =====================================================
  Rn data_Un(Vh.NbDoF(), 0.);
  interpolate(Vh, data_Un, fun_init_solution);
  Fun2_h fun_Un(Vh, data_Un);
  // projection(fun_Un, fun_init_solution);

  FunTest u(Vh,1), v(Vh,1), Un(fun_Un, 0, 1);

  // FOR THE LIMITING PROCEDURE
  // =====================================================
  Limiter limiter;
  Rn data_limiting_Un(Vh.NbDoF(), 0.);


  //COMPUTE THE MASS MATRIX
  // =====================================================
  burger.setMatrixTo(Ah);
  // burger.addFaceStabilization(
  //     - innerProduct(jump(u), Cstab*jump(v))
  //     - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
  // );
  burger.setMatrixTo(Mh);   // to fill Mh
  burger.addBilinear(
           innerProduct(u, v)
         );
  // burger.addFaceStabilization(
  //        innerProduct(h*jump(u), Cstab*jump(v))
  //      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
  //        );
  burger.reset();   // unset point Mh

  for(int iter=0;iter<niteration;++iter) {

    std::cout << " Iteration " << iter << " / " << niteration << std::endl;
    std::cout << " Time      " << tid << " / " << Tend << std::endl;

    // DEFINE FUNCTION AND EXPRESSION FOR THE FLUX
    Expr   expr_Un(fun_Un, 0, op_id);
    Fun2_h fun_flux_Un0(Fh, 0.5*expr_Un*expr_Un, 0.5*expr_Un*expr_Un);

    Fun2_h fun_Ex(Eh, fun_solution, tid);
    Expr   expr_Ex(fun_Ex, 0, op_id);
    Fun2_h fun_flux_Ex(Fh, 0.5*expr_Ex*expr_Ex, 0.5*expr_Ex*expr_Ex);

    FunTest Uex(fun_Ex,0,1), flux_Ex(fun_flux_Ex, 0, 2);
    //
    limiter.KXRCF_indicator(fun_Un, fun_flux_Un0);
    Fun2_h indicator(Ih, limiter.indicator);
    limiter.limiting(fun_Un, data_limiting_Un);
    //
    // Fun2_h limiting_Un(Vh, data_limiting_Un);

    data_Un = data_limiting_Un;
    Fun2_h fun_flux_Un(Fh, 0.5*expr_Un*expr_Un, 0.5*expr_Un*expr_Un);
    FunTest flux_Un(fun_flux_Un, 0, 2);

    if( MPIcf::IamMaster() && iter%5==0 ) {
      Paraview2 writer(Vh, "burgerFEM_"+to_string(ifig++)+".vtk");
      writer.add(fun_Un, "uh"  , 0, 1);
      // writer.add(fun_Ex, "u_ex", 0, 1);
      // writer.add(indicator, "indicator", 0, 1);
      // writer.add(limiting_Un, "u_limiting" ,0 ,1);
      // writer.add(fun_flux_Un, "flux", 0, 2);
    }

    double t1 = MPIcf::Wtime();

    // CONSTRUCT THE RHS
    burger.addLinear(
        innerProduct(flux_Un  , grad(v))
    );
    burger.addLinear(
      - innerProduct(average(flux_Un*N), jump(v))
      - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
      , innerEdge
    );
    // burger.addLinear(
    //   - innerProduct(average(flux_Un*N), jump(v))
    //   - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
    //   // - jump(innerProduct(flux_Un*N,v))
    //   // - innerProduct(jump(flux_Un*N), jump(lambda*v))
    //   , interface
    // );
    //
    // BOUNDARY CONDITION IN RHS
    burger.addLinear(
      - innerProduct(flux_Un, (0.5*v)*N)
      - innerProduct(Un    , lambdaB*0.5*v)
      , boundary
      , {4}          // boundary in
    );
    burger.addLinear(
      - innerProduct(flux_Un*N  , v )
      , boundary
      , {1,2,3}        // boundary out
    );

    MatriceMap<double> mAh(burger.nDoF, burger.nDoF, Ah);
    mAh.addMatMul(data_Un, burger.rhs);

    //  BOUNDARY CONDITION  IN LHS
    burger.addLinear(
      - innerProduct(flux_Ex, (0.5*v)*N)
      + innerProduct(Uex    , lambdaB*0.5*v)
      , boundary
      , {4}          // boundary in
    );
    std::cout << " Time assembly RHS => " << MPIcf::Wtime()-t1 << std::endl;

    burger.clearMatrix = false;
    burger.solve(Mh, burger.rhs);
    data_Un += dt * burger.rhs;
    burger.rhs = 0.;




    // if(MPIcf::IamMaster() &&  iter+1 == niteration) {
    //   Paraview2 writer(Wh, levelSet, "burger2_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_Un, "uh", 0, 1);
    //   // writer.add(flux_Un, "flux", 0, 2);
    //   writer.add(fun_Ex, "u_ex", 0, 1);
    //   // writer.add(indicator, "Ih", 0, 1);
    // }

    tid += dt;
    // return 0;
  }


  };


#endif

#ifdef EXAMPLE1
R fun_levelSet(const R2 P, const int i) {
  return - P.x - P.y + 1.66;
}
R fun_init_solution(const R2 P, int elementComp, int domain) {
  return sin(pi*(P.x+P.y));
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  double C = P.x+P.y;
  double x1 = 0;
  double x2 = 4;
  double x0 = 0.5*(x1+x2);
  double f1 = x1 + 2.*sin(pi*x1)*t - C;
  if(fabs(f1) < 1e-10) {
    return sin(pi*x1);
  }
  double f2 = x2 + 2.*sin(pi*x2)*t - C;
  if(fabs(f2) < 1e-10) {
    return sin(pi*x2);
  }

  for(int i=0;i<=50;++i) {
    double fa = x1 + 2.*sin(pi*x1)*t - C;
    double fm = x0 + 2.*sin(pi*x0)*t - C;

    if(fa*fm < 0) x2 = x0;
    else          x1 = x0;

    if(fabs(fm) < 1e-8) break;
    x0 = 0.5*(x1+x2);

    if(i == 50) {
      std::cout << " f(x0) = " << fm << std::endl;
      assert(0);
      std::cout << " Burger Newton not converged" << std::endl;
    }
  }
  return sin(pi*x0);
}


// R fun_solution(const R2 P, double un, double t) {
//   double C = P.x+P.y;
//   // if(un > 1) un = 1;
//   // else if (un < -1) un = -1;
//   // double x0 = C - 2*un*t;
//   double x0 = C;
//   double xx0 = x0;
//
//   if(P.x == 0.875 && P.y == 0.025) std::cout << "x0 " << x0 << std::endl;
//   for(int i=0;i<=20;++i) {
//     double Fx  = x0 + 2.*sin(pi*x0)*t - C;
//     double dFx = 1. + 2.*pi*cos(pi*x0)*t;
//
//     double dx = Fx / dFx;
//     x0 -= dx;
//
//     if(P.x == 0.875 && P.y == 0.025) std::cout << x0 << "\t" << dx << std::endl;
//
//     if(fabs(dx) < 1e-8) break;
//     if(i == 10) {
//       std::cout << " dx = " << dx << std::endl;
//       std::cout << " xx0 = " << xx0 << std::endl;
//       std::cout << un << "\t" << asin(un) << std::endl;
//       std::cout << " x+y = " << C << "\n"
//                 << " x0  = " << x0 << std::endl;
//       std::cout << P << std::endl;
//       std::cout << " Burger Newton not converged" << std::endl;
//       assert(0);
//     }
//   }
//   return sin(pi*x0);
// }
// void fun_solution(Fun2_h& f,  const Rn& un, double t) {
//   const FESpace2& Vh(*f.Vh);
//   assert(un.size() == Vh.nbDoF);
//   double dataSend[Vh.nbDoF];
//   Rn_ fhSend(dataSend, Vh.nbDoF); fhSend = 1e+50;
//   for (int k=Vh.first_element();k<Vh.last_element(); k+= Vh.next_element()) {
//     const typename FESpace2::FElement& FK((Vh)[k]);
//     const int domain = FK.whichDomain();
//     const int kb = Vh.idxElementInBackMesh(k);
//     for(int j=FK.dfcbegin(0);j<FK.dfcend(0);++j) {
//       R2 mip = FK.Pt(j);
//       fhSend(FK.loc2glb(j)) = fun_solution(mip, un(FK.loc2glb(j)), t);
//     }
//   }
//   MPIcf::AllReduce(dataSend, f.data, fhSend.size(),MPI_MIN);
// }


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  int nx = 51;
  int ny = 51;
  Mesh2 Th(nx, ny, 0., 0., 2., 2.);
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;
  int ifig = 0;
  double tid = 0.;
  double Tend = 0.7/pi;
  int niteration = Tend / dt;
  dt = Tend / niteration;


  // DEFINITION OF PARAMETERS
  // =====================================================
  double Cstab = 1e-2;
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",1, 1);
  CutFEM_Parameter lambda("lambda",0, 1.);
  double lambdaB = 1.;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  Fun2_h levelSet(Lh, fun_levelSet);
  Interface2 interface(Th, levelSet.v, 1);

  // CONSTRUCTION OF THE FE SPACES AND THE CUT SPACE FOR Uh AND FLUX
  // =====================================================
  FESpace2 Vh(Th, DataFE<Mesh2>::P1dc);
  CutFESpace2 Wh(Vh, interface, {1, -1});

  KN<const GTypeOfFE<Mesh2>* > arrayFE_Flux(2);
  arrayFE_Flux(0) = &DataFE<Mesh2>::P1dc;
  arrayFE_Flux(1) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> burgerFE_Flux(arrayFE_Flux);
  FESpace2 Fh(Th, burgerFE_Flux);
  CutFESpace2 WFh(Fh, interface, {1, -1});

  // For the indicator
  FESpace2 Ih(Th, DataFE<Mesh2>::P0);
  CutFESpace2 WIh(Ih, interface, {1, -1});


  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh2> burger(Wh);
  Normal N;
  MatMap Mh, Ah;

  // FUNCTION Un AND TESTFUNCTION
  // =====================================================
  Rn data_Un(Wh.NbDoF(), 0.);
  interpolate(Wh, data_Un, fun_init_solution);
  Fun2_h fun_Un(Wh, data_Un);
  FunTest u(Wh,1), v(Wh,1), Un(fun_Un, 0, 1);


  // FOR THE LIMITING PROCEDURE
  // =====================================================
  Limiter limiter;
  Rn data_limiting_Un(Wh.NbDoF(), 0.);


  //COMPUTE THE MASS MATRIX
  // =====================================================
  burger.setMatrixTo(Ah);
  burger.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
  );
  burger.setMatrixTo(Mh);   // to fill Mh
  burger.addBilinear(
           innerProduct(u, v)
         );
  burger.addFaceStabilization(
         innerProduct(h*jump(u), Cstab*jump(v))
       + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
         );
  burger.reset();   // unset point Mh

  for(int iter=0;iter<niteration;++iter) {

    std::cout << " Iteration " << iter << " / " << niteration << std::endl;
    std::cout << " Time      " << tid << " / " << Tend << std::endl;

    // DEFINE FUNCTION AND EXPRESSION FOR THE FLUX
    Expr   expr_Un(fun_Un, 0, op_id);
    Fun2_h fun_flux_Un(WFh, 0.5*expr_Un*expr_Un, 0.5*expr_Un*expr_Un);
    FunTest flux_Un(fun_flux_Un, 0, 2);

    Fun2_h fun_Ex(Wh, fun_solution, tid);
    // fun_solution(fun_Ex, data_Un, tid);
    Expr   expr_Ex(fun_Ex, 0, op_id);
    Fun2_h fun_flux_Ex(WFh, 0.5*expr_Ex*expr_Ex, 0.5*expr_Ex*expr_Ex);
    FunTest Uex(fun_Ex,0,1), flux_Ex(fun_flux_Ex, 0, 2);

    limiter.KXRCF_indicator(fun_Un, fun_flux_Un);
    Fun2_h indicator(WIh, limiter.indicator);
    // if(iter > 100){
      limiter.limiting(fun_Un, data_limiting_Un);
    // }
    Fun2_h limiting_Un(Wh, data_limiting_Un);


    // std::cout << data_Un << std::endl;
    // std::cout << fun_Un.v << std::endl;

    // std::cout << data_limiting_Un << std::endl;


    if( MPIcf::IamMaster() && iter%5==0 ) {
      Paraview2 writer(Wh, levelSet, "burger_"+to_string(ifig++)+".vtk");
      writer.add(fun_Un, "uh"  , 0, 1);
      writer.add(fun_Ex, "u_ex", 0, 1);
      writer.add(indicator, "indicator", 0, 1);
      writer.add(limiting_Un, "u_limiting" ,0 ,1);
      // writer.add(fun_flux_Un, "flux", 0, 2);
    }

    double t1 = MPIcf::Wtime();

    // CONSTRUCT THE RHS
    burger.addLinear(
        innerProduct(flux_Un  , grad(v))
    );
    burger.addLinear(
      - innerProduct(average(flux_Un*N), jump(v))
      - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
      , innerEdge
    );
    burger.addLinear(
      - innerProduct(average(flux_Un*N), jump(v))
      - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
      // - jump(innerProduct(flux_Un*N,v))
      // - innerProduct(jump(flux_Un*N), jump(lambda*v))
      , interface
    );

    // BOUNDARY CONDITION IN RHS
    burger.addLinear(
      - innerProduct(flux_Un, (0.5*v)*N)
      - innerProduct(Un    , lambdaB*0.5*v)
      , boundary
      , {1,4}          // boundary in
    );
    burger.addLinear(
      - innerProduct(flux_Un*N  , v )
      , boundary
      , {2,3}        // boundary out
    );

    MatriceMap<double> mAh(burger.nDoF, burger.nDoF, Ah);
    mAh.addMatMul(data_Un, burger.rhs);

    //  BOUNDARY CONDITION  IN LHS
    burger.addLinear(
      - innerProduct(flux_Ex, (0.5*v)*N)
      + innerProduct(Uex    , lambdaB*0.5*v)
      , boundary
      , {1, 4}          // boundary in
    );
    std::cout << " Time assembly RHS => " << MPIcf::Wtime()-t1 << std::endl;

    burger.clearMatrix = false;
    burger.solve(Mh, burger.rhs);
    data_Un += dt * burger.rhs;
    burger.rhs = 0.;




    if(MPIcf::IamMaster() &&  iter+1 == niteration) {
      Paraview2 writer(Wh, levelSet, "burger_"+to_string(ifig++)+".vtk");
      writer.add(fun_Un, "uh", 0, 1);
      // writer.add(flux_Un, "flux", 0, 2);
      writer.add(fun_Ex, "u_ex", 0, 1);
      writer.add(indicator, "Ih", 0, 1);

    }

    tid += dt;
  }


  };


#endif
#ifdef EXAMPLE2

R fun_levelSet(const R2 P, const int i) {
  return -P.x -P.y + 0.5;
}
R fun_init_solution(const R2 P, int elementComp, int domain) {
  return (P.x <= 0);
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  double x_shock = 0.5*t;
  return (P.x <= x_shock);
}



int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  int nx = 51;
  int ny = 51;
  Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;
  int ifig = 0;
  double tid = 0.;
  double Tend = 0.1;
  int niteration = Tend / dt;
  dt = Tend / niteration;


  // DEFINITION OF PARAMETERS
  // =====================================================
  double Cstab = 1e-2;
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",1, 1);
  CutFEM_Parameter lambda("lambda",0, 1.);
  double lambdaB = 1.;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  Fun2_h levelSet(Lh, fun_levelSet);
  Interface2 interface(Th, levelSet.v, 1);

  // CONSTRUCTION OF THE FE SPACES AND THE CUT SPACE FOR Uh AND FLUX
  // =====================================================
  FESpace2 Eh(Th, DataFE<Mesh2>::P1dc);
  CutFESpace2 WEh(Eh, interface, {1, -1});


  FESpace2 Vh(Th, DataFE<Mesh2>::P1dcOrt);
  CutFESpace2 Wh(Vh, interface, {1, -1});

  KN<const GTypeOfFE<Mesh2>* > arrayFE_Flux(2);
  arrayFE_Flux(0) = &DataFE<Mesh2>::P1dc;
  arrayFE_Flux(1) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> burgerFE_Flux(arrayFE_Flux);
  FESpace2 Fh(Th, burgerFE_Flux);
  CutFESpace2 WFh(Fh, interface, {1, -1});

  // For the indicator
  FESpace2 Ih(Th, DataFE<Mesh2>::P0);
  CutFESpace2 WIh(Ih, interface, {1, -1});


  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh2> burger(Wh);
  Normal N;
  MatMap Mh, Ah;

  // FUNCTION Un AND TESTFUNCTION
  // =====================================================
  Rn data_Un(Wh.NbDoF(), 0.);
  // interpolate(Wh, data_Un, fun_init_solution);
  Fun2_h fun_Un(Wh, data_Un);
  projection(fun_Un, fun_init_solution);
  // for(int i=0;i<data_Un.size();i+=3){
  //   if(fabs(data_Un(i))>1e-8) data_Un(i) = 1.;
  // }


  FunTest u(Wh,1), v(Wh,1), Un(fun_Un, 0, 1);


  // FOR THE LIMITING PROCEDURE
  // =====================================================
  Limiter limiter;
  Rn data_limiting_Un(Wh.NbDoF(), 0.);


  //COMPUTE THE MASS MATRIX
  // =====================================================
  burger.setMatrixTo(Ah);
  burger.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
  );
  burger.setMatrixTo(Mh);   // to fill Mh
  burger.addBilinear(
           innerProduct(u, v)
         );
  burger.addFaceStabilization(
         innerProduct(h*jump(u), Cstab*jump(v))
       + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
         );
  burger.reset();   // unset point Mh

  for(int iter=0;iter<niteration;++iter) {

    std::cout << " Iteration " << iter << " / " << niteration << std::endl;
    std::cout << " Time      " << tid << " / " << Tend << std::endl;

    // DEFINE FUNCTION AND EXPRESSION FOR THE FLUX
    Expr   expr_Un(fun_Un, 0, op_id);
    Fun2_h fun_flux_Un0(WFh, 0.5*expr_Un*expr_Un, 0.5*expr_Un*expr_Un);

    Fun2_h fun_Ex(WEh, fun_solution, tid);
    Expr   expr_Ex(fun_Ex, 0, op_id);
    Fun2_h fun_flux_Ex(WFh, 0.5*expr_Ex*expr_Ex, 0.5*expr_Ex*expr_Ex);
    FunTest Uex(fun_Ex,0,1), flux_Ex(fun_flux_Ex, 0, 2);

    limiter.KXRCF_indicator(fun_Un, fun_flux_Un0);
    Fun2_h indicator(WIh, limiter.indicator);
    limiter.limiting(fun_Un, data_limiting_Un);

    Fun2_h limiting_Un(Wh, data_limiting_Un);

    if( MPIcf::IamMaster() && iter%2==0 ) {
      Paraview2 writer(Wh, levelSet, "burger2_"+to_string(ifig++)+".vtk");
      writer.add(fun_Un, "uh"  , 0, 1);
      writer.add(fun_Ex, "u_ex", 0, 1);
      writer.add(indicator, "indicator", 0, 1);
      writer.add(limiting_Un, "u_limiting" ,0 ,1);
      // writer.add(fun_flux_Un, "flux", 0, 2);
    }


    data_Un = data_limiting_Un;
    Fun2_h fun_flux_Un(WFh, 0.5*expr_Un*expr_Un, 0.5*expr_Un*expr_Un);
    FunTest flux_Un(fun_flux_Un, 0, 2);

    double t1 = MPIcf::Wtime();

    // CONSTRUCT THE RHS
    burger.addLinear(
        innerProduct(flux_Un  , grad(v))
    );
    burger.addLinear(
      - innerProduct(average(flux_Un*N), jump(v))
      - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
      , innerEdge
    );
    burger.addLinear(
      - innerProduct(average(flux_Un*N), jump(v))
      - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
      // - jump(innerProduct(flux_Un*N,v))
      // - innerProduct(jump(flux_Un*N), jump(lambda*v))
      , interface
    );

    // BOUNDARY CONDITION IN RHS
    burger.addLinear(
      - innerProduct(flux_Un, (0.5*v)*N)
      - innerProduct(Un    , lambdaB*0.5*v)
      , boundary
      , {4}          // boundary in
    );
    burger.addLinear(
      - innerProduct(flux_Un*N  , v )
      , boundary
      , {1,2,3}        // boundary out
    );

    MatriceMap<double> mAh(burger.nDoF, burger.nDoF, Ah);
    mAh.addMatMul(data_Un, burger.rhs);

    //  BOUNDARY CONDITION  IN LHS
    burger.addLinear(
      - innerProduct(flux_Ex, (0.5*v)*N)
      + innerProduct(Uex    , lambdaB*0.5*v)
      , boundary
      , {4}          // boundary in
    );
    std::cout << " Time assembly RHS => " << MPIcf::Wtime()-t1 << std::endl;

    burger.clearMatrix = false;
    burger.solve(Mh, burger.rhs);
    data_Un += dt * burger.rhs;
    burger.rhs = 0.;




    // if(MPIcf::IamMaster() &&  iter+1 == niteration) {
    //   Paraview2 writer(Wh, levelSet, "burger2_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_Un, "uh", 0, 1);
    //   // writer.add(flux_Un, "flux", 0, 2);
    //   writer.add(fun_Ex, "u_ex", 0, 1);
    //   // writer.add(indicator, "Ih", 0, 1);
    // }

    tid += dt;
  }


  };




#endif
