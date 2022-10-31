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
#define EXAMPLE2

R fun_levelSet1(const R2 P, const int i) {
  double dy = 0.23;
  return P.y - dy;
}
R fun_levelSet2(const R2 P, const int i) {
  double dy = 0.77;
  return - P.y + dy;
}
R fun_levelSet(const R2 P, const int i) {
  return fun_levelSet1(P,i)*fun_levelSet2(P,i);
}

#ifdef EXAMPLE1
int x_0 = 0;
int x_L = -1;
R fun_boundary(const R2 P, int elementComp, double t) {
  return 2*sin(pi*(P.x+P.y-3*t));
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  if(domain == 0) return 2*sin(pi*(P.x+P.y-3*t));
  else return 3*sin(3./2*pi*(P.x-2*t+2./3*P.y-x_0/3));
}
R fluxExact_X(const R2 P, int elementComp, int domain, double t) {
  double u = fun_solution(P, 0, domain, t);
  return 3*u;
}
R fluxExact_Y(const R2 P, int elementComp, int domain, double t) {
  double u = fun_solution(P, 0, domain, t);
  return 0;
}
#endif
#ifdef EXAMPLE2
double c0 = 0.5;
R fun_boundary(const R2 P, int elementComp, double t) {
  return sin(pi*(P.x + P.y- 4*t));
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
  else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
}
R fluxExact_X(const R2 P, int elementComp, int domain, double t) {
  double u = fun_solution(P, 0, domain, t);
  return 3*u;
}
R fluxExact_Y(const R2 P, int elementComp, int domain, double t) {
  double u = fun_solution(P, 0, domain,t);
  return u;
}
#endif

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE BACKGROUND MESH
  // =====================================================
  int nx = 101;
  int ny = 101;
  // Mesh2 Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Mesh2 Th(nx, ny, 0., 0., 1.5, 1.5);
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;
  int niteration = 1;
  int ifig = 0;
  double tid = 0;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  Fun2_h levelSetDown(Lh, fun_levelSet1);
  Interface2 interface(Th, levelSetDown.v, 1);
  Fun2_h levelSetUp(Lh, fun_levelSet2);
  interface.add(levelSetUp.v, 3);

  Rn levelSet_data(levelSetDown.v.size(), 0.);
  levelSet_data = levelSetDown.v * levelSetUp.v;
  Fun2_h levelSet(Lh, levelSet_data);

  // Fun2_h levelSet(Lh, fun_levelSet);
  // Interface2 interface(Th, levelSet.v);


  // CONSTRUCTION OF THE FE SPACES AND THE CUT SPACE
  KN<const GTypeOfFE<Mesh2>* > arrayFE(1);
  arrayFE(0) = &DataFE<Mesh2>::P1dc;
  GTypeOfFESum<Mesh2> eulerFE(arrayFE);
  FESpace2 Vh(Th, eulerFE);
  CutFESpace2 Wh(Vh, interface, {1});


  // DECLARATION OF THE VECTOR CONTAINING THE SOLUTION
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_solution, 0.);
  if(MPIcf::IamMaster()) {
    Fun2_h solex(Wh, fun_solution, tid);
    Fun2_h sol(Wh, u0);
    Paraview2 writer(Wh, levelSet, "hyperbolicDG.vtk");
    writer.add(sol  , "uh" , 0, 1);
    writer.add(levelSet, "levelSet", 0, 1);
  }

  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh2> problem(Wh);
  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);
  MatMap Mh, Ah;
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",3, 2);
  CutFEM_Parameter lambda("lambda",0, 1.);
  #ifdef EXAMPLE1
  CutFEM_R2 beta("beta", R2(3,0), R2(2,0));
  #else
  CutFEM_R2 beta("beta", R2(3,1), R2(2,1));
  #endif
  double lambdaB = 3;
  double Cstab = 1e-2;

  Fun2_h Un(Wh, u0);
  Expr fun_Un(Un, 0, op_id);
  FunTest f_Un(Wh, fun_Un);
  FunTest flux_Un = beta*f_Un;


  // ASSEMBLY CONSTANT PART
  problem.pmat = &Ah;
  problem.addFaceStabilization(
    - innerProduct(jump(u), Cstab*jump(v))
    - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
  );


  problem.pmat = &Mh;
  problem.addBilinear(
    innerProduct(u,v)
  );
  problem.addFaceStabilization(
    innerProduct(h*jump(u), Cstab*jump(v))
    + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
  );

  for(int i=0;i<niteration;++i) {
    Fun2_h gh(Wh, fun_solution, tid);
    Expression2 gx(gh, 0, op_id);
    FunTest ggx(Wh, gx);
    Fun2_h F(Wh, fluxExact_X, tid);
    Fun2_h G(Wh, fluxExact_Y, tid);
    Expression2 Fex_u(F, 0, op_id);
    Expression2 Gex_u(G, 0, op_id);
    FunTest FluxEx_Un(Wh, Fex_u, Gex_u);

    // BUILDING FLUX PART
    // =====================================================
    problem.addLinear(
      innerProduct(flux_Un, grad(v))
    );

    problem.addLinear(
      - innerProduct(average(flux_Un*n), jump(v))
      - innerProduct(0.5*lambdaE*jump(f_Un)  , jump(v))
      , innerEdge
    );
    #ifdef EXAMPLE2
    problem.addLinear(
      -innerProduct(flux_Un*n, v)
      , interface
      , {3}          // label other boundary
    );
    problem.addLinear(
      - innerProduct(flux_Un, (0.5*v)*n)
      - innerProduct(f_Un   , lambdaB*0.5*v)
      , interface
      , {1}          // label left boundary
    );
    #endif
    problem.addLinear(
      -innerProduct(flux_Un*n, v)
      , boundary
      , {2}          // label other boundary
    );

    problem.addLinear(
      - innerProduct(flux_Un, (0.5*v)*n)
      - innerProduct(f_Un   , lambdaB*0.5*v)
      , boundary
      , {4}          // label left boundary
    );


    // Ah*u0 to get stab part
    MatriceMap<double> mAh(problem.nDoF, problem.nDoF, Ah);
    mAh.addMatMul(u0, problem.rhs);
    #ifdef EXAMPLE2
    problem.addLinear(
      - innerProduct(FluxEx_Un, (0.5*v)*n)
      + innerProduct(ggx      , lambdaB*0.5*v)
      , interface
      , {1}          // label left boundary
    );
    #endif
    problem.addLinear(
      - innerProduct(FluxEx_Un, (0.5*v)*n)
      + innerProduct(gx       , lambdaB*0.5*v)
      , boundary
      , {4}          // label left boundary
    );


    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.clearMatrix = false;
    problem.solve(Mh, problem.rhs);

    problem.rhs *= dt;
    problem.rhs += u0;

    u0 = problem.rhs;


    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
      // Fun2_h solex(Wh, fun_solution, tid);
      for(int j=0;j<u0.size();++j) {
        if(fabs(u0(j)) < 1e-16 ) u0(j) = 0.;
      }
      Fun2_h sol(Wh, u0);
      Paraview2 writer(Wh, levelSet, "hyperbolicDG_Example2_"+to_string(ifig++)+".vtk");
      writer.add(sol, "uh", 0, 1);
      // writer.add(solex, "uex", 0, 1);
    }



  }

};
