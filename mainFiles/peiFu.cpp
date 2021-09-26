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

#define EXAMPLE3

#ifdef EXAMPLE1

int x_0 = 0;
int x_L = -1;
R fun_levelSet(const R2 P, const int i) {
  return -P.x;
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return 2*sin(pi*(P.x+P.y-3*t));
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  // return 2*sin(pi*(P.x+P.y-3*t));
  if(domain == 0) return 2*sin(pi*(P.x+P.y-3*t));
  else return 3*sin(3./2*pi*(P.x-2*t+2./3*P.y-x_0/3));
}

void assembly(const FESpace2& Wh, const Interface2& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,0), R2(2,0));
  CutFEM_Parameter lambda("lambda",0., 1.);

  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",3, 2);

  double lambdaB = 3;
  double Cstab = 1e-2;

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);
  TestFunction2 Hessu = grad(grad(du)), Hessv = grad(grad(v));


  // BUILDING A
  // =====================================================
  problem.pmat = &Ah;
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
  );

  problem.addBilinear(
    - jump(innerProduct(beta*u,v*n))
    - innerProduct(jump(beta*u*n), jump(lambda*v))
    , interface
  );
  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      - innerProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
  );

  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addEdgeIntegral(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)   , jump(v))
  );

  problem.addBilinearFormBorder(
     -innerProduct(beta*u*n, v)
    , {2}          // label other boundary
  );
  problem.addBilinearFormBorder(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , {4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.pmat = &Mh;
  problem.addBilinear(
    innerProduct(u,v)
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      + innerProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const FESpace2& Wh, const Interface2& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  // FEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,0), R2(2,0));
  // CutFEM_Parameter lambda("lambda",0., 1.);
  //
  // const CutFEM_Parameter& h(Parameter::h);
  // CutFEM_Parameter lambdaE("lambdaE",3, 2);

  double lambdaB = 3;
  // double Cstab = 1e-2;

  Fun2_h gh(Wh, fun_solution, tn);
  Expression2 gx(gh, 0, op_id);

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  problem.mat = Ah;
  problem.addMatMul(u0);
  problem.addLinearFormBorder(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
    , {4}          // label left boundary
  );
  problem.cleanMatrix();

  problem.mat = Mh;
  std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve();

  uh = problem.rhs;

}

#endif
#ifdef EXAMPLE2
double c0 = 0.5;
R fun_levelSet(const R2 P, const int i) {
  return -P.x - P.y + c0;
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return sin(pi*(P.x + P.y- 4*t));
}
R fun_solution(const R2 P, int elementComp, int domain, double t) {
  // return 2*sin(pi*(P.x+P.y-3*t));
  if(domain == 0) return sin(pi*(P.x + P.y- 4*t));
  else return 4./3*sin(4./3*pi*(P.x + P.y - 3*t - c0/4));
}

void assembly(const FESpace2& Wh, const Interface2& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,1), R2(2,1));
  CutFEM_Parameter lambda("lambda",0., 1.);

  const CutFEM_Parameter& h(Parameter::h);
  // double h = 2./40;
  CutFEM_Parameter lambdaE("lambdaE",3, 2);  // max B = sqrt(10)/sqrt(5)

  double lambdaB = 3;
  double Cstab = 5e-1;
  double Cstabt = 5e-1;

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  TestFunction2 Hessu = grad(grad(u)), Hessv = grad(grad(v));
  // std::cout << Hessu << std::endl;
  // std::cout << jump(Hessu) << std::endl;
  // getchar();
  // BUILDING A
  // =====================================================
  problem.pmat = &Ah;
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
  );

  problem.addBilinear(
    - jump(innerProduct(beta*u,v*n))
    - innerProduct(jump(beta*u*n), jump(lambda*v))
    , interface
  );
  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addEdgeIntegral(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
  );

  problem.addBilinearFormBorder(
     -innerProduct(beta*u*n, v)
    , {2, 3}          // label other boundary
  );
  problem.addBilinearFormBorder(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.pmat = &Mh;
  problem.addBilinear(
    innerProduct(u,v)
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const FESpace2& Wh, const Interface2& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,1), R2(2,1));
  double lambdaB = 3;

  Fun2_h gh(Wh, fun_solution, tn);
  Expression2 gx(gh, 0, op_id);

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  problem.mat = Ah;
  problem.addMatMul(u0);
  problem.addLinearFormBorder(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
    , {1,4}          // label left boundary
  );
  problem.cleanMatrix();

  problem.mat = Mh;
  std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve();

  uh = problem.rhs;

}
#endif
#ifdef EXAMPLE3
double c0 = 0.25;
R fun_levelSet(const R2 P, const int i) {
  return -P.x - P.y + c0;
}
// R fun_boundary(const R2 P, int elementComp, double t) {
//   return ;
// }
R fun_solution(const R2 P, int elementComp, int domain, double t) {

  if(domain == 1) return 0;
  double xs = -0.3, ys = -0.3;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0;
}

void assembly(const FESpace2& Wh, const Interface2& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,1), R2(1,3));
  CutFEM_Parameter lambda("lambda",0.4, 1.);

  const CutFEM_Parameter& h(Parameter::h);
  // double h = 2./40;
  CutFEM_Parameter lambdaE("lambdaE",sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)

  double lambdaB = sqrt(10);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  TestFunction2 Hessu = grad(grad(u)), Hessv = grad(grad(v));
  // std::cout << Hessu << std::endl;
  // std::cout << jump(Hessu) << std::endl;
  // getchar();
  // BUILDING A
  // =====================================================
  problem.pmat = &Ah;
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
  );

  problem.addBilinear(
    - jump(innerProduct(beta*u,v*n))
    - innerProduct(jump(beta*u*n), jump(lambda*v))
    , interface
  );
  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addEdgeIntegral(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
  );

  problem.addBilinearFormBorder(
     -innerProduct(beta*u*n, v)
    , {2, 3}          // label other boundary
  );
  problem.addBilinearFormBorder(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.pmat = &Mh;
  problem.addBilinear(
    innerProduct(u,v)
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const FESpace2& Wh, const Interface2& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,1), R2(1,3));
  double lambdaB = sqrt(10);

  // Fun2_h gh(Wh, fun_solution, tn);
  // Expression2 gx(gh, 0, op_id);

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  problem.mat = Ah;
  problem.addMatMul(u0);
  // problem.addLinearFormBorder(
  //    - innerProduct(gx, beta*  (0.5*v)*n)
  //    + innerProduct(gx, lambdaB*0.5*v)
  //   , {1,4}          // label left boundary
  // );
  problem.cleanMatrix();

  problem.mat = Mh;
  std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve();

  uh = problem.rhs;

}
#endif
int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE MESH
  // =====================================================
  int nx = 100;
  int ny = 100;
  Mesh2 Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 1.;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;
  double qu0 = 0;
  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  Fun2_h levelSet(Lh, fun_levelSet);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  Interface2 interface(Th, levelSet.v);
  FESpace2 Vh(Th, DataFE<Mesh2>::P1dc);
  // FESpace2 Wh(Th, interface, DataFE<Mesh>::P1dc);
  CutFESpace2 Wh(Vh, interface, {1,-1});


  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_solution, 0.);
  Rn uh(u0);

  {
    Fun2_h femSolh(Wh, uh);
    Expression2 femSol(femSolh, 0, op_id);
    qu0 = integral(femSol, Wh);
  }

  if(MPIcf::IamMaster()) {
    Fun2_h solex(Wh, fun_solution, tid);
    Fun2_h sol(Wh, uh);
    Paraview2 writer(Wh, levelSet, "peiFu_ex2_0.vtk");
    writer.add(sol  , "uh" , 0, 1);
    // writer.add(solex, "uex", 0, 1);
  }

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {

    Rn u1(u0);
    Rn u2(Wh.NbDoF(), 0.);

    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, u1, tid);
    // uh += dt * u1;

    // THIRD ORDER RK
    // =================================================
    solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    u1 += dt * uh;
    u2 += 3./4*u0 + 1./4*u1;
    solve_problem(Wh, interface, u1, uh, Ah, Mh, tid+dt);
    u2 += 1./4 * dt * uh;
    solve_problem(Wh, interface, u2, uh, Ah, Mh, tid+0.5*dt);
    uh *= 2./3 * dt;
    uh += 1./3*u0 + 2./3*u2;

    u0 = uh;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // interpolate(Wh, u1, fun_solution, tid);
    // u1 -= uh;
    // Fun2_h femSolh(Wh, u1);
    // Expression2 femSol(femSolh, 0, op_id);

    R errU = 0;//sqrt(integral(femSol*femSol, Wh));
    errSum += errU;

    Fun2_h femSolh(Wh, uh);
    Expression2 femSol(femSolh, 0, op_id);
    R qu = integral(femSol, Wh);



    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%10 == 0) {
      // Fun2_h solex(Wh, fun_solution, tid);
      for(int i=0;i<uh.size();++i) {
        if(fabs(uh(i)) < 1e-16 ) uh(i) = 0.;
      }
      Fun2_h sol(Wh, uh);
      Paraview2 writer(Wh, levelSet, "peiFu_ex2_"+to_string(ifig++)+".vtk");
      writer.add(sol, "uh", 0, 1);
      // writer.add(solex, "uex", 0, 1);
    }

    std::cout << "Iteration " << i << " / " << niteration
              << " \t time = " << tid << " , || u-uex ||_2 = " << errU << std::endl;
    std::cout << setprecision(16) << "q(u) = " << qu
              << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;

  }

  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}




/*
int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef Interface2 Interface;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE MESH
  // =====================================================
  int nx = 160;
  int ny = 160;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = 0.2 * meshSize / 3 ;
  double tend = 1;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  FESpace Lh(Th, DataFE<Mesh>::P1);
  Fun2_h levelSet(Lh, fun_levelSet);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  Interface2 interface(Th, levelSet.v);
  FESpace2 Vh(Th, DataFE<Mesh>::P1dc);
  // FESpace2 Wh(Th, interface, DataFE<Mesh>::P1dc);   // try stab on artificial problem
  CutFESpace2 Wh(Vh, interface, {1,-1});


  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_solution, 0.);
  Rn uh(u0);

  // if(MPIcf::IamMaster()) {
  //   Fun2_h solex(Wh, fun_solution, tid);
  //   Fun2_h sol(Wh, uh);
  //   Paraview2 writer(Wh, levelSet, "peiFu_0.vtk");
  //   writer.add(sol  , "uh" , 0, 1);
  //   writer.add(solex, "uex", 0, 1);
  // }


  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {

    Rn u1(u0);
    Rn u2(Wh.NbDoF(), 0.);

    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, u1, tid);
    // uh += dt * u1;

    // THIRD ORDER RK
    // =================================================
    solve_problem(Wh, interface, u0, uh, tid);
    u1 += dt * uh;
    u2 += 3./4*u0 + 1./4*u1;
    solve_problem(Wh, interface, u1, uh, tid+dt);
    u2 += 1./4 * dt * uh;
    solve_problem(Wh, interface, u2, uh, tid+0.5*dt);
    uh *= 2./3 * dt;
    uh += 1./3*u0 + 2./3*u2;

    u0 = uh;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    interpolate(Wh, u1, fun_solution, tid);
    u1 -= uh;
    Fun2_h femSolh(Wh, u1);
    Expression2 femSol(femSolh, 0, op_id);

    R errU = sqrt(integral(femSol*femSol, Wh));
    errSum += errU;


    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%5 == 0) {
    //   Fun2_h solex(Wh, fun_solution, tid);
    //   Fun2_h sol(Wh, uh);
    //   Paraview2 writer(Wh, levelSet, "peiFu_"+to_string(ifig++)+".vtk");
    //   writer.add(sol, "uh", 0, 1);
    //   writer.add(solex, "uex", 0, 1);
    // }

    std::cout << "Iteration " << i << " / " << niteration
              << " \t time = " << tid << " , || u ||_2 = " << errU << std::endl;
  }

  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}
*/

/*
void solve_problem(const FESpace2& Wh, const Interface2& interface, const Rn& u0, Rn& uh, double tn) {

  double t0 = MPIcf::Wtime();

  CutFEM<Mesh2> problem(Wh);
  // FEM<Mesh2> problem(Wh);
  CutFEM_R2 beta("beta", R2(3,0), R2(2,0));
  CutFEM_Parameter lambda("lambda",0., 1.);

  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter lambdaE("lambdaE",3, 2);

  double lambdaB = 3;
  double Cstab = 1e-2;

  Fun2_h gh(Wh, fun_solution, tn);
  Expression2 gx(gh, 0, op_id);

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);

  // BUILDING A
  // =====================================================
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
  );

  problem.addBilinear(
    // - innerProduct(average(beta*u*n), jump(v))
    // - innerProduct(0.5*lambdaE*jump(u)   , jump(v))
    - jump(innerProduct(beta*u,v*n))
    - innerProduct(jump(beta*u*n), jump(lambda*v))
    , interface
  );
  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
  );

  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addEdgeIntegral(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)   , jump(v))
  );

  problem.addBilinearFormBorder(
     -innerProduct(beta*u*n, v)
    , {2}          // label other boundary
  );
  problem.addBilinearFormBorder(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , {4}          // label left boundary
  );
  problem.addLinearFormBorder(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
    , {4}          // label left boundary
  );

  // matlab::Export(problem.mat, "matEdge.dat");

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  double tt0 = CPUtime();
  problem.addMatMul(u0);
  // std::cout << " Time  A*u0 \t\t" << CPUtime() - tt0 << std::endl;
  problem.cleanMatrix();



  // BUILDING (M + S)
  // =====================================================
  problem.addBilinear(
    innerProduct(u,v)
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
  );

  std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;
  // SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve();

  uh = problem.rhs;

}
*/
