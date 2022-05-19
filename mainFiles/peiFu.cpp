#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <omp.h>
#include "../util/cputime.h"
#include "cfmpi.hpp"


#include "baseProblem.hpp"
#include "paraview.hpp"

typedef std::map<std::pair<int,int>,R> MatMap;
typedef Mesh2 Mesh;
typedef FESpace2 Space;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;

#define EXAMPLE3

// Example to check convergence
#ifdef EXAMPLE1
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

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(3,1), R2(2,1)});
  CutFEMParameter lambda(0., 1.);

  const MeshParameter& h(Parameter::h);
  CutFEMParameter lambdaE(3, 2);  // max B = sqrt(10)/sqrt(5)

  double lambdaB = 3;
  double Cstab = 5e-1;
  double Cstabt = 5e-1;

  Normal n;
  TestFunction2 u(Wh,1), v(Wh,1);


  TestFunction2 Hessu = grad(grad(u)), Hessv = grad(grad(v));
  // BUILDING A
  // =====================================================
  problem.set_map(Ah);
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
    , Khi
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
      , Khi
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , innerFacet
  );

  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , boundary
    , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , boundary
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.set_map(Mh);
  problem.addBilinear(
    innerProduct(u,v)
    , Khi
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
      , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const Space& Wh, const Interface<Mesh>& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(3,1), R2(2,1)});
  double lambdaB = 3;

  Fun_h gh(Wh, fun_solution, tn);
  Expression2 gx(gh, 0, op_id);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // problem.mat = Ah;
  // problem.addMatMul(u0);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0, problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , boundary
    , {1,4}          // label left boundary
  );
  // std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve(Mh, problem.rhs_);

  uh = problem.rhs_;

}



int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // OUTPUT FILE
  // =====================================================
  std::ofstream outputData("output_test.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 50;
  int ny = 50;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.5;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;
  double qu0 = 0;
  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  Space Lh(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  InterfaceLevelSet<Mesh> interface(Th, levelSet);
  ActiveMesh<Mesh> Khi(Th, interface);
  CutSpace Wh(Khi, Vh);

  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_solution, 0.);
  Fun_h Un(Wh, u0);
  Rn uh(u0);

    // Fun2_h femSolh(Wh, uh);


  // if(MPIcf::IamMaster()) {
  //   Fun2_h solex(Wh, fun_solution, tid);
  //   Fun2_h sol(Wh, uh);
  //   Paraview2 writer(Wh, levelSet, "testDetector.vtk");
  //   writer.add(sol  , "uh" , 0, 1);
  //   // writer.add(levelSet, "levelSet", 0, 1);
  // }

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {


    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    // std::cout << uh << std::endl;
    //
    // uh *= dt;
    // uh += u0;

    Rn u1(u0);
    Rn u2(Wh.NbDoF(), 0.);
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
    interpolate(Wh, u1, fun_solution, tid);
    u1 -= uh;
    Fun_h femErrh(Wh, u1);
    Expression2 femErr(femErrh, 0, op_id);

    R errU = sqrt(integral(Khi, femErr*femErr));
    errSum += errU;

    // Expression2 femSol(femSolh, 0, op_id);
    // R qu = integral(Khi, femSol);



    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%10 == 0 || i+1 == niteration) {
      Fun2_h femSolh(Wh, uh);
      // Fun2_h solex(Wh, fun_solution, tid);
      for(int j=0;j<uh.size();++j) {
        if(fabs(uh(j)) < 1e-16 ) uh(j) = 0.;
      }
      Fun2_h sol(Wh, uh);
      Paraview<Mesh> writer(Khi, "hyperbolicLawExample1_"+to_string(ifig++)+".vtk");
      writer.add(sol, "uh", 0, 1);
      // writer.add(solex, "uex", 0, 1);
    }
    std::cout << "Iteration " << i << " / " << niteration
              << " \t time = " << tid << " , || u-uex ||_2 = " << errU << std::endl;
    // outputData << i << "\t"
    //            << tid << "\t"
    //            << setprecision(16) << qu << "\t"
    //            << setprecision(16) << fabs(qu-qu0) << "\t"
    //            << setprecision(5) << std::endl;

    // return 0;
  }


  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

// Example to check conservation
#ifdef EXAMPLE2
double c0 = 0.25;
R fun_levelSet(const R2 P, const int i) {
  return -P.x - P.y + c0;
}
R fun_initial(const R2 P, int elementComp, int domain) {

  if(domain == 1) return 0;
  double xs = -0.3, ys = -0.3;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0;
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(5));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);
  // double h = 2./40;

  double lambdaB = sqrt(10);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;

  Normal n;
  FunTest u(Wh,1), v(Wh,1);

  // BUILDING A
  // =====================================================
  problem.set_map(Ah);
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
    ,Khi
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
      ,Khi
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , innerFacet
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , boundary
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , boundary
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.set_map(Mh);
  problem.addBilinear(
    innerProduct(u,v)
    , Khi
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
      , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space& Wh, const Interface<Mesh>& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  const ActiveMesh<Mesh> Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  double lambdaB = sqrt(10);


  Normal n;
  FunTest u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,uh);

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve(Mh, uh);//problem.rhs_);

  // uh = problem.rhs_;
}


int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // OUTPUT FILE
  // =====================================================
  std::ofstream outputData("output_test.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 200;
  int ny = 200;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.4;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;
  double qu0 = 0;
  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  Space Lh(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  InterfaceLevelSet<Mesh> interface(Th, levelSet);
  ActiveMesh<Mesh> Khi(Th, interface);
  CutSpace Wh(Khi, Vh);

  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_initial);
  Fun_h Un(Wh, u0);
  Rn uh(u0);


  Fun_h fun_uh(Wh, uh);
  // Expression2 femSol(femSolh, 0, op_id);
  qu0 = integral(Khi, fun_uh, 0);


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
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    R qu = integral(Khi, fun_uh, 0);


    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%5 == 0 || i+1 == niteration) {
    //   // Fun2_h solex(Wh, fun_solution, tid);
    //   for(int j=0;j<uh.size();++j) {
    //     if(fabs(uh(j)) < 1e-16 ) uh(j) = 0.;
    //   }
    //   // Fun2_h sol(Wh, uh);
    //   Paraview<Mesh> writer(Khi, "conservation_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uh", 0, 1);
    //   // writer.add(solex, "uex", 0, 1);
    // }
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(5) << std::endl;
  }



  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif


// Example to check macro extension P0
#ifdef EXAMPLE3
double c0 = 0.25;
R fun_levelSet(const R2 P, const int i) {
  return -P.x - P.y + c0;
}
R fun_initial(const R2 P, int elementComp, int domain) {

  if(domain == 1) return 0;
  double xs = -0.15, ys = -0.15;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0;
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  // CutFEM_R2 beta({R2(1,2), R2(1,2)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(5));  // max B = sqrt(10)/sqrt(5)
  // CutFEMParameter lambdaE(sqrt(5), sqrt(5));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);
  // double h = 2./40;

  double lambdaB = sqrt(10);
  // double lambdaB = sqrt(5);

  double Cstab = 1e-2;
  double Cstabt = 1e-2;

  Normal n;
  FunTest u(Wh,1), v(Wh,1);

  // BUILDING A
  // =====================================================
  problem.set_map(Ah);
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
    ,Khi
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
      ,Khi
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , innerFacet
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , boundary
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , boundary
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.set_map(Mh);
  problem.addBilinear(
    innerProduct(u,v)
    , Khi
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
      , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space& Wh, const Interface<Mesh>& interface,const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  const ActiveMesh<Mesh> Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  double lambdaB = sqrt(10);
  // CutFEM_R2 beta({R2(1,2), R2(1,2)});
  // double lambdaB = sqrt(5);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);


  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,uh);

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.solve(Mh, uh);//problem.rhs_);

  // uh = problem.rhs_;
}


int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // OUTPUT FILE
  // =====================================================
  std::ofstream outputData("output_test.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 200;
  int ny = 200;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.3;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  Space Lh(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  InterfaceLevelSet<Mesh> interface(Th, levelSet);
  ActiveMesh<Mesh> Khi(Th, interface);
  CutSpace Wh(Khi, Vh);

  MacroElement<Mesh> macro(Khi, 0.5);

  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);
  interpolate(Wh, u0, fun_initial);
  Fun_h Un(Wh, u0);
  Rn uh(u0);
  Rn u1(u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_u1(Wh, u1);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "peiII.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
  //   writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge2.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge2.vtk");
  // }


  double qu0 = integral(Khi, fun_uh, 0);
  double min_u0 = uh.min();
  double max_u0 = uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {

    // Rn u1(Wh.NbDoF(), 0.);

    // EULER METHOD
    // =================================================
    solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    uh *= dt;
    uh += u0;


    u1 = uh;
    // limiter::extendToMacro(fun_uh, u1, macro);

    u0 = u1;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = integral(Khi, fun_u1, 0);
    // double min_uh = uh.min();
    // double max_uh = uh.max();
    // double min_u1 = u1.min();
    // double max_u1 = u1.max();
    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    //   std::cout << " Maximum principle satified! " <<  std::endl;
    //   std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
    //   std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
    //   std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    // else{
    //   std::cout << " Maximum principle not satified! " << std::endl;
    //   std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
    //   std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
    //   std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;
    //
    //   // if(MPIcf::IamMaster()) {
    //   //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
    //   //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   //   writer.add(fun_u1, "uhLimiter", 0, 1);
    //   // }
    //   return 0;
    // }
    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%1000 == 0 || i+1 == niteration) {

      Paraview<Mesh> writer(Khi, "conservationP0_"+to_string(ifig++)+".vtk");
      writer.add(fun_uh, "uhNoLimiter", 0, 1);
      writer.add(fun_u1, "uhLimiter", 0, 1);
    }

    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(5) << std::endl;

    // if(i==10) return 0;
  }



  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif





/*
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
  TestFunction2 Hessu = grad(grad(u)), Hessv = grad(grad(v));


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
  CutFEM_R2 beta("beta2", R2(3,0), R2(2,0));
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
  MatriceMap<double> mAh(problem.nDoF, problem.nDoF, Ah);
  mAh.addMatMul(u0, problem.rhs);
  // problem.mat = Ah;
  // problem.addMatMul(u0);
  problem.addLinearFormBorder(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
    , {4}          // label left boundary
  );


  // problem.cleanMatrix();
  // problem.mat = Mh;
  std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;
  // getchar();
// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  problem.clearMatrix = false;
  problem.solve(Mh, problem.rhs);
  std::cout << Mh.size() << "\t" << problem.mat.size() << std::endl;
  uh = problem.rhs;

}

#endif
#ifdef EXAMPLE4
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
void assembly(const CutSpace& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());
  CutFEM<Mesh> problem(Wh);
  // std::vector beta_v = {R2(3,1), R2(1,2)};
  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(5));

  const MeshParameter& h(Parameter::h);
  double lambdaB = sqrt(10);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;

  Normal n;
  FunTest u(Wh,1), v(Wh,1);


  // FunTest Hessu = grad(grad(u)), Hessv = grad(grad(v));
  problem.set_map(Ah);
  // problem.pmat = &Ah;
  problem.addBilinear(
    innerProduct(beta * u, grad(v))
    , Khi
  );

  problem.addBilinear(
    - jump(innerProduct(beta*u,v*n))
    - innerProduct(jump(beta*u*n), jump(lambda*v))
    , interface
  );
  // //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
      , Khi
  );

  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , innerFacet
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , boundary
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , boundary
    , {1, 4}          // label left boundary
  );

  // BUILDING (M + S)
  // =====================================================
  problem.set_map(Mh);
  // problem.pmat = &Mh;
  problem.addBilinear(
    innerProduct(u,v)
    , Khi
  );
  problem.addFaceStabilization(
        innerProduct(h*jump(u), Cstab*jump(v))
      + innerProduct((h^3)*jump(grad(u)), Cstab*jump(grad(v)))
      // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
      , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const CutSpace& Wh, const Interface<Mesh>& interface, Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {
  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh> Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;

  CutFEM<Mesh2> problem(Wh, optionProblem);
  CutFEM_R2 beta({R2(3,1), R2(1,2)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(5));  // max B = sqrt(10)/sqrt(5)
  double lambdaB = sqrt(10);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  Fun_h Un(Wh, u0);
  // Expression2 fun_Un(Un, 0, op_id);
  // FunTest f_Un(Wh, fun_Un);
  // FunTest flux_Un = beta*f_Un;

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // problem.addLinear(
  //   innerProduct(flux_Un, grad(v))
  // );
  //
  // problem.addLinear(
  //   - jump(innerProduct(flux_Un,v*n))
  //   - innerProduct(jump(flux_Un*n), jump(lambda*v))
  //   , interface
  // );
  //
  // problem.addLinear(
  //   - innerProduct(average(flux_Un*n), jump(v))
  //   - innerProduct(0.5*lambdaE*jump(f_Un)  , jump(v))
  //   , innerEdge
  // );
  // problem.addLinear(
  //    -innerProduct(flux_Un*n, v)
  //    , boundary
  //    , {2, 3}          // label other boundary
  // );
  //
  // problem.addLinear(
  //   - innerProduct(f_Un, beta*(0.5*v)*n)
  //   - innerProduct(f_Un, lambdaB*0.5*v)
  //   , boundary
  //   , {1, 4}          // label left boundary
  // );

  // matlab::Export(problem.rhs, "linearRHS.dat");


  // MatriceMap<double> mAh(problem.nDoF, problem.nDoF, Ah);
  // mAh.addMatMul(u0, problem.rhs);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,uh);//problem.rhs_);
  // matlab::Export(problem.rhs, "referenceRHS.dat");

  // problem.addLinear(
  //    - innerProduct(gx, beta*  (0.5*v)*n)
  //    + innerProduct(gx, lambdaB*0.5*v)
  //    , Khi
  //    , boundary
  //   , {1,4}          // label left boundary
  // );
//   // problem.cleanMatrix();
//
//   // problem.mat = Mh;
  // std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

// SOLVING  (M+S)E'(t) = rhs
  // =====================================================

  problem.solve(Mh, uh);//problem.rhs_);
  // uh = problem.rhs_;
}


#endif
*/
//
// int main(int argc, char** argv ) {
//
//   MPIcf cfMPI(argc,argv);
//   const double cpubegin = MPIcf::Wtime();
//
//   // OUTPUT FILE
//   // =====================================================
//   std::ofstream outputData("output_test.txt");
//
//   // DEFINITION OF THE MESH and SPACE
//   // =====================================================
//   int nx = 100;
//   int ny = 100;
//   Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
//   Space Vh(Th, DataFE<Mesh2>::P1dc);
//
//   // DEFINITION OF SPACE AND TIME PARAMETERS
//   // =====================================================
//   double tid = 0;
//   double meshSize = 2./nx;
//   double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
//   double tend = 0.5;
//   int niteration = tend / dt;
//   dt = tend / niteration;
//   double errSum = 0;
//   double qu0 = 0;
//   std::cout << "Mesh size h = \t" << meshSize << std::endl;
//   std::cout << "Time step dt = \t" << dt << std::endl;
//
//   // DEFINITION OF THE LEVELSET
//   // =====================================================
//   Space Lh(Th, DataFE<Mesh2>::P1);
//   Fun_h levelSet(Lh, fun_levelSet);
//
//   // CONSTRUCTION INTERFACE AND CUTSPACE
//   // =====================================================
//   InterfaceLevelSet<Mesh> interface(Th, levelSet);
//   ActiveMesh<Mesh> Khi(Th, interface);
//   CutSpace Wh(Khi, Vh);
//
//   // DECLARATION OF THE VECTOR CONTAINING THE solution
//   // =====================================================
//   Rn u0(Wh.NbDoF(), 0.);
//   interpolate(Wh, u0, fun_solution, 0.);
//   Fun_h Un(Wh, u0);
//   Rn uh(u0);
//
//   {
//     Fun2_h femSolh(Wh, uh);
//     Expression2 femSol(femSolh, 0, op_id);
//     qu0 = integral(Khi, femSol);
//   }
//
//   // if(MPIcf::IamMaster()) {
//   //   Fun2_h solex(Wh, fun_solution, tid);
//   //   Fun2_h sol(Wh, uh);
//   //   Paraview2 writer(Wh, levelSet, "testDetector.vtk");
//   //   writer.add(sol  , "uh" , 0, 1);
//   //   // writer.add(levelSet, "levelSet", 0, 1);
//   // }
//
//   // ASSEMBLY THE CONSTANT PART
//   // ==================================================
//   MatMap Ah, Mh;
//   assembly(Wh, interface, Ah, Mh);
//
//   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
//   // ==================================================
//   int ifig = 1;
//   for(int i=0;i<niteration;++i) {
//
//
//     // EULER METHOD
//     // =================================================
//     // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
//     // std::cout << uh << std::endl;
//     //
//     // uh *= dt;
//     // uh += u0;
//
//     Rn u1(u0);
//     Rn u2(Wh.NbDoF(), 0.);
//     // THIRD ORDER RK
//     // =================================================
//     solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
//     u1 += dt * uh;
//     u2 += 3./4*u0 + 1./4*u1;
//     solve_problem(Wh, interface, u1, uh, Ah, Mh, tid+dt);
//     u2 += 1./4 * dt * uh;
//     solve_problem(Wh, interface, u2, uh, Ah, Mh, tid+0.5*dt);
//     uh *= 2./3 * dt;
//     uh += 1./3*u0 + 2./3*u2;
//
//
//     u0 = uh;
//     tid += dt;
//
//     // COMPUTATION OF THE L2 ERROR
//     // =================================================
//     // Rn u1(u0);
//     // Rn u1(Wh.NbDoF(), 0.);
//     interpolate(Wh, u1, fun_solution, tid);
//     u1 -= uh;
//     Fun2_h femErrh(Wh, u1);
//     Expression2 femErr(femErrh, 0, op_id);
//
//     R errU = sqrt(integral(Khi, femErr*femErr));
//     errSum += errU;
//
//     Fun2_h femSolh(Wh, uh);
//     Expression2 femSol(femSolh, 0, op_id);
//     R qu = integral(Khi, femSol);
//
//
//
//     // PLOT THE SOLUTION
//     // ==================================================
//     // if(MPIcf::IamMaster() && i%5 == 0 || i+1 == niteration) {
//     //   // Fun2_h solex(Wh, fun_solution, tid);
//     //   for(int j=0;j<uh.size();++j) {
//     //     if(fabs(uh(j)) < 1e-16 ) uh(j) = 0.;
//     //   }
//     //   Fun2_h sol(Wh, uh);
//     //   Paraview<Mesh> writer(Khi, "peiFu_"+to_string(ifig++)+".vtk");
//     //   writer.add(sol, "uh", 0, 1);
//     //   // writer.add(solex, "uex", 0, 1);
//     // }
//     std::cout << "Iteration " << i << " / " << niteration
//               << " \t time = " << tid << " , || u-uex ||_2 = " << errU << std::endl;
//     std::cout << setprecision(16) << "q(u) = " << qu
//               << std::endl;
//     std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
//               << setprecision(6) << std::endl;
//
//     outputData << i << "\t"
//                << tid << "\t"
//                << setprecision(16) << qu << "\t"
//                << setprecision(16) << fabs(qu-qu0) << "\t"
//                << setprecision(5) << std::endl;
//
//     // return 0;
//   }
//
//
//
//   // KN<const GTypeOfFE<Mesh2>* > arrayFE_Flux(2);
//   // arrayFE_Flux(0) = &DataFE<Mesh2>::P1dc;
//   // arrayFE_Flux(1) = &DataFE<Mesh2>::P1dc;
//   // GTypeOfFESum<Mesh2> peiFE_Flux(arrayFE_Flux);
//   // Space Fh(Th, peiFE_Flux);
//   // CutFESpace2 WFh(Fh, interface, {1, -1});
//
//
//   // CutFEM_Parameter betaX("betaX", 3,2);
//   // CutFEM_Parameter betaY("betaY", 0,0);
//   // Expression2   expr_Un(Un, 0, op_id);
//   // Fun2_h flux_Un(WFh, betaX*expr_Un, betaY*expr_Un);
//   //
//   // Limiter limiter;
//   // limiter.KXRCF_indicator(Un, flux_Un);
//   // std::cout << limiter.indicator << std::endl;
//
//
//
//   // if(MPIcf::IamMaster()) {
//   //   Paraview2 writer(Wh, levelSet, "testFlux.vtk");
//   //   writer.add(flux_Un, "flux", 0, 2);
//   //   // writer.add(fun_Ex, "u_ex", 0, 1);
//   // }
//
//
//   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
//   std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
//   return 0;
// }
