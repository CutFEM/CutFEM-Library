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
#include "../num/matlab.hpp"


typedef std::map<std::pair<int,int>,R> MatMap;
typedef Mesh2 Mesh;
typedef FESpace2 Space;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;
typedef ExpressionFunFEM<Mesh2> Expression;

// #define FEM_SMOOTH_SOLUTION
// #define CUTFEM_SMOOTH_SOLUTION
// #define FEM_DISCONTINUOUS_SOLUTION
// #define CUTFEM_DISCONTINUOUS_SOLUTION
// #define P1dc_FEM_LIMITER
// #define P1dc_CutFEM_LIMITER
// #define P0_CutFEM_LIMITER

// Test for the paper
// #define CUTFEM_LIMITER_ACCURACY_TEST
#define CUTFEM_LIMITER_BURGER_TEST


#ifdef FEM_SMOOTH_SOLUTION

R fun_initial(const R2 P, int elementComp) {
  return 0.5+0.5*sin(pi*(P.x+P.y));
}
R fun_boundary(const R2 P, int elementComp, double t) {
    return 0.5+0.5*sin(pi*(P.x + P.y- 4*t));
}
void assembly(const Space& Wh, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const Mesh& Khi(Wh.Th);
  FEM<Mesh2> problem(Wh);
  R2 beta(3,1);
  double lambdaE = sqrt(10);
  const MeshParameter& h(Parameter::h);
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


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , INTEGRAL_INNER_EDGE_2D
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , INTEGRAL_BOUNDARY
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , INTEGRAL_BOUNDARY
    , {1, 4}          // label left boundary
  );


  // BUILDING (M + S)
  // =====================================================
  problem.set_map(Mh);
  problem.addBilinear(
    innerProduct(u,v)
    , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space& Wh, const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  const Mesh& Khi(Wh.Th);

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  FEM<Mesh2> problem(Wh, optionProblem);
  R2 beta(3,1);
  double lambdaB = sqrt(10);
  Normal n;
  FunTest u(Wh,1), v(Wh,1);

  Fun_h gh(Wh, fun_boundary, tn);
  Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , INTEGRAL_BOUNDARY
    , {1,4}          // label left boundary
  );

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
  std::ofstream outputData("output_testFEM_P0.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 40;
  int ny = 40;
  Mesh Kh(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Kh, DataFE<Mesh2>::P0);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;


  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Vh.NbDoF(), 0.);
  interpolate(Vh, u0, fun_initial);

  Rn uh(u0);
  Rn uh_tild(u0);
  Fun_h fun_u0(Vh, u0);
  Fun_h fun_uh(Vh, uh);
  Fun_h fun_uh_tild(Vh, uh_tild);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "peiII.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
  //   writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  // }

  double qu0 = 0.;//integral(Khi, fun_uh, 0);
  double min_u0 = 0.;//uh.min();
  double max_u0 = 1.;//uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Vh,Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    // solve_problem(Vh, u0, uh, Ah, Mh, tid);
    // uh *= dt;
    // uh += u0;
    // limiter::FEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0);


    // THIRD ORDER RK
    // =================================================
    Rn u1(u0);
    Rn u2(Vh.NbDoF(), 0.);
    // Rn u2_tild(Vh.NbDoF(), 0.);
    // Rn u1_tild(Vh.NbDoF(), 0.);
    Fun_h fun_u1(Vh, u1);
    Fun_h fun_u2(Vh, u2);


    // if()
    // solve_problem(Vh, u0, uh, Ah, Mh, tid);
    // u1 += dt * uh;
    // limiter::FEM::limiter_Pei_P1(fun_u1, u1_tild, min_u0, max_u0);
    // u2 += 3./4*u0 + 1./4*u1_tild;
    // solve_problem(Vh, u1_tild, uh, Ah, Mh, tid+dt);
    // u2 += 1./4 * dt * uh;
    // limiter::FEM::limiter_Pei_P1(fun_u2, u2_tild, min_u0, max_u0);
    // solve_problem(Vh, u2_tild, uh, Ah, Mh, tid+0.5*dt);
    // uh *= 2./3 * dt;
    // uh += 1./3*u0 + 2./3*u2_tild;
    // limiter::FEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0);


    solve_problem(Vh, u0, uh, Ah, Mh, tid);
    u1 += dt * uh;
    // limiter::CutFEM::extendToMacroP0(fun_u1, u1_tild, macro);
    // u1_tild = u1;
    u2 += 3./4*u0 + 1./4*u1;
    solve_problem(Vh, u1, uh, Ah, Mh, tid+dt);
    u2 += 1./4 * dt * uh;
    // limiter::CutFEM::extendToMacroP0(fun_u2, u2_tild, macro);
    // u2_tild = u2;
    solve_problem(Vh, u2, uh, Ah, Mh, tid+0.5*dt);
    uh *= 2./3 * dt;
    uh += 1./3*u0 + 2./3*u2;
    // limiter::CutFEM::extendToMacroP0(fun_uh, uh_tild, macro);
    uh_tild = uh;


    // solve_problem(Vh, u0, uh, Ah, Mh, tid);
    // u1 += dt * uh;
    // u2 += 3./4*u0 + 1./4*u1;
    // solve_problem(Vh, u1, uh, Ah, Mh, tid+dt);
    // u2 += 1./4 * dt * uh;
    // solve_problem(Vh, u2, uh, Ah, Mh, tid+0.5*dt);
    // uh *= 2./3 * dt;
    // uh += 1./3*u0 + 2./3*u2;
    // uh_tild = uh;

    // initialize solution
    // =================================================
    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = 0;//integral(Khi, fun_u1, 0);
    double min_uh, max_uh;
    limiter::FEM::minmaxP1(fun_uh, min_uh, max_uh);
    double min_uh_tild, max_uh_tild;
    limiter::FEM::minmaxP1(fun_uh_tild, min_uh_tild, max_uh_tild);

    limiter::FEM::check_mean_value(fun_uh_tild, 0., 1.);

    // if((i==24) && MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    //   writer.add(fun_uM, "macroExtend", 0, 1);
    // }

    // if(MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Kh, "maxPrincipleFEM_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    // }
    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    if( min_u0 <= min_uh_tild+Epsilon && max_uh_tild+Epsilon <= max_u0) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h}  < " << max_uh <<  std::endl;
      std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h}  < " << max_uh <<  std::endl;
      std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild <<  std::endl;

      // std::cout << min_u0 - min_uh_tild << std::endl;
      // std::cout << max_u0 - max_uh_tild << std::endl;


      // return 0;
    }

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    Rn usol(Vh.NbDoF(), 0.);
    interpolate(Vh, usol, fun_boundary, tid);
    Rn uerr(usol);
    uerr -= uh_tild;
    Fun_h femErrh(Vh, uerr);
    Fun_h fun_ex (Vh, usol);

    Expression2 femErr(femErrh, 0, op_id);
    R errU = sqrt(integral(Kh, femErr*femErr));
    errSum += errU;

    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%10 == 0 || i+1 == niteration) {

      Paraview<Mesh> writer(Kh, "smooth_FEM_P0_"+to_string(ifig++)+".vtk");
      writer.add(fun_uh, "uhNoLimiter", 0, 1);
      writer.add(fun_u1, "uhLimiter", 0, 1);
      writer.add(fun_ex , "u_exact", 0, 1);

    }

    std::cout << " || u-uex ||_2 = " << errU << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;


    if(i == 0) {
      outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU " << std::endl;

    }
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << errU << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(16) << min_uh_tild << "\t"
               << setprecision(16) << max_uh_tild << "\t"
               << setprecision(5) << std::endl;

    // return 0;
    // if(i==10) return 0;
  }


  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

#ifdef CUTFEM_SMOOTH_SOLUTION
// double c0 = 0.77;//0.25;
// R fun_levelSet(const R2 P, const int i) {
//   return -0.5*P.x - P.y + 0.67*c0;
// }
// R fun_initial(const R2 P, int elementComp, int domain) {
//   return 0.5+0.5*sin(pi*(P.x+P.y));
// }
// R fun_boundary(const R2 P, int elementComp, int domain, double t) {
//     return 0.5+0.5*sin(pi*(P.x + P.y- 4*t));
// }
// R fun_solution(const R2 P, int elementComp, int domain, double t) {
//     return 0.5+0.5*sin(pi*(P.x + P.y- 4*t));
// }
R fun_levelSet(const R2 P, const int i) {
  return P.y - 0.5*P.x - 1./4;
}

R fun_initial (const R2 P, int elementComp, int domain) {
  return sin(pi*(P.x+P.y));
}
R fun_solution (const R2 P, int elementComp, int domain, double t) {
  return sin(pi*(P.x+P.y - 4*t));
}
R fun_boundary (const R2 P, int elementComp, int domain, double t) {
  return sin(pi*(P.x+P.y - 4*t));
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(3,1), R2(3,1)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);

  double lambdaB = sqrt(10);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;
  // Fun_h beta(Wh, fun_velocity);
  // Expression2 beta(fun_beta, 0, op_id);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);



  // BUILDING A
  // =====================================================
  problem.set_map(Ah);
  problem.addBilinear(
    innerProduct(beta*u, grad(v))
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
    , INTEGRAL_INNER_EDGE_2D
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , INTEGRAL_BOUNDARY
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , INTEGRAL_BOUNDARY
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

  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(3,1), R2(3,1)});
  double lambdaB = sqrt(10);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  Fun_h gh(Wh, fun_boundary, tn);
  Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , INTEGRAL_BOUNDARY
    , {1,4}          // label left boundary
  );

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
  std::ofstream outputData("output_testP1.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 80;
  int ny = 80;
  // Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Mesh2 Th(nx, ny, 0., 0., 2., 2.);

  Space Vh(Th, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./(nx+1);
  double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  Space Lh(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  Lagrange2 FE_beta(1);
  Space Uh(Th, FE_beta);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  InterfaceLevelSet<Mesh> interface(Th, levelSet);
  ActiveMesh<Mesh> Khi(Th, interface);
  // Khi.truncate(interface, 1);
  CutSpace Wh(Khi, Vh);
  CutSpace Qh(Khi, Uh);

  MacroElement<Mesh> macro(Khi, 0.5);


  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.); interpolate(Wh, u0, fun_initial);
  Rn uh(u0);
  Rn uh_tild(u0);

  Fun_h fun_u0(Wh, u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_uh_tild(Wh, uh_tild);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "test_accuracy_mesh.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge1.vtk");
  //   // writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
  //   // writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  // }
  // return 0;

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
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    // uh *= dt;
    // uh += u0;
    // // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0, macro);
    // limiter::CutFEM::extendToMacroP0(fun_uh, uh_tild, macro);


    // std::map<int, double> u_mean;
    // uM = u0;
    // uh = u0;
    // limiter::extendToMacroP1(fun_uh, uM, u_mean, macro);
    // limiter::check_maximum_principle<Mesh>(u_mean, 0, 1);
    // The U_mean should satisfy the maximum principle
    // u0 = uh_tild;

    // THIRD ORDER RK
    // =================================================
    {
      Rn u1(u0);
      Rn u2(Wh.NbDoF(), 0.);
      Rn u2_tild(Wh.NbDoF(), 0.);
      Rn u1_tild(Wh.NbDoF(), 0.);
      Fun_h fun_u1(Wh, u1);
      Fun_h fun_u2(Wh, u2);

      if(Wh.basisFctType == BasisFctType::P1dc) {
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild, min_u0, max_u0, macro);
        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild, min_u0, max_u0, macro);
        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild, min_u0, max_u0, macro);
      }
      else if (Wh.basisFctType == BasisFctType::P0){

        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
      }
      else {
        assert(0);
      }
    }

    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = integral(Khi, fun_uh_tild, 0);
    double min_uh, max_uh;
    limiter::CutFEM::findMinAndMaxValue(fun_uh, min_uh, max_uh);
    double min_u1, max_u1;
    limiter::CutFEM::findMinAndMaxValue(fun_uh_tild, min_u1, max_u1);

    // if(MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "smoothSol_CutFEM_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    //   // writer.add(fun_uM, "macroExtend", 0, 1);
    // }

    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    if( min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0+Epsilon) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;

      std::cout << min_u0 - min_u1 << std::endl;
      std::cout << max_u0 - max_u1 << std::endl;


      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }
      // return 0;
    }


    // COMPUTATION OF THE L2 ERROR
    // =================================================
    Rn usol(Wh.NbDoF(), 0.);
    interpolate(Wh, usol, fun_solution, tid);
    Rn uerr(usol);
    uerr -= uh_tild;
    Fun_h femErrh(Wh, uerr);
    Fun_h fun_ex (Wh, usol);

    Expression2 femErr(femErrh, 0, op_id);
    R errU = sqrt(integral(Khi, femErr*femErr));
    errSum += errU;

    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
    //
    //   Paraview<Mesh> writer(Khi, "smoothSolCutFEM_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    //   writer.add(femErrh, "error"  , 0, 1);
    //   writer.add(fun_ex , "u_exact", 0, 1);
    //
    // }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << " || u-uex ||_2 = " << errU << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << errU << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(5) << std::endl;

  }


  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif


#ifdef P1dc_FEM_LIMITER

R fun_initial(const R2 P, int elementComp) {
  double xs = -0.1, ys = -0.1;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0.;
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return 0.;
}
void assembly(const Space& Wh, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const Mesh& Khi(Wh.Th);
  FEM<Mesh2> problem(Wh);
  R2 beta(1,1);
  double lambdaE = sqrt(10);
  const MeshParameter& h(Parameter::h);
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

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space& Wh, const Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  const Mesh& Khi(Wh.Th);

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  FEM<Mesh2> problem(Wh, optionProblem);
  R2 beta(1,1);
  double lambdaB = sqrt(10);
  Normal n;
  FunTest u(Wh,1), v(Wh,1);

  Fun_h gh(Wh, fun_boundary, tn);
  Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , boundary
    , {1,4}          // label left boundary
  );

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
  std::ofstream outputData("output_testP1.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 100;
  int ny = 100;
  Mesh Kh(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Kh, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(2) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;


  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Vh.NbDoF(), 0.);
  interpolate(Vh, u0, fun_initial);

  Rn uh(u0);
  Rn uh_tild(u0);
  Fun_h fun_u0(Vh, u0);
  Fun_h fun_uh(Vh, uh);
  Fun_h fun_uh_tild(Vh, uh_tild);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "peiII.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
  //   writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  // }

  double qu0 = 0.;//integral(Khi, fun_uh, 0);
  double min_u0 = uh.min();
  double max_u0 = uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Vh,Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    solve_problem(Vh, u0, uh, Ah, Mh, tid);
    uh *= dt;
    uh += u0;
    limiter::FEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0);


    // // THIRD ORDER RK
    // // =================================================
    // Rn u1(u0);
    // Rn u2(Vh.NbDoF(), 0.);
    // Rn u2_tild(Vh.NbDoF(), 0.);
    // Rn u1_tild(Vh.NbDoF(), 0.);
    // Fun_h fun_u1(Vh, u1);
    // Fun_h fun_u2(Vh, u2);
    //
    // solve_problem(Vh, u0, uh, Ah, Mh, tid);
    // u1 += dt * uh;
    // limiter::FEM::limiter_Pei_P1(fun_u1, u1_tild, min_u0, max_u0);
    //
    // u2 += 3./4*u0 + 1./4*u1_tild;
    // solve_problem(Vh, u1_tild, uh, Ah, Mh, tid+dt);
    // u2 += 1./4 * dt * uh;
    // limiter::FEM::limiter_Pei_P1(fun_u2, u2_tild, min_u0, max_u0);
    //
    // solve_problem(Vh, u2_tild, uh, Ah, Mh, tid+0.5*dt);
    // uh *= 2./3 * dt;
    // uh += 1./3*u0 + 2./3*u2_tild;
    // limiter::FEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0);
    //


    // initialize solution
    // =================================================
    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = 0;//integral(Khi, fun_u1, 0);
    double min_uh, max_uh;
    limiter::FEM::minmaxP1(fun_uh, min_uh, max_uh);
    double min_uh_tild, max_uh_tild;
    limiter::FEM::minmaxP1(fun_uh_tild, min_uh_tild, max_uh_tild);

    // if((i==24) && MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    //   writer.add(fun_uM, "macroExtend", 0, 1);
    // }

    if(MPIcf::IamMaster()) {
      Paraview<Mesh> writer(Kh, "maxPrincipleFEM_"+to_string(ifig++)+".vtk");
      writer.add(fun_uh, "uhNoLimiter", 0, 1);
      writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    }
    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    if( min_u0 <= min_uh_tild && max_uh_tild <= max_u0) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h}  < " << max_uh <<  std::endl;
      std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h}  < " << max_uh <<  std::endl;
      std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild <<  std::endl;

      std::cout << min_u0 - min_uh_tild << std::endl;
      std::cout << max_u0 - max_uh_tild << std::endl;


      // return 0;
    }


    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
    //
    //   Paraview<Mesh> writer(Khi, "conservationP1_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    // }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    // std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              // << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(5) << std::endl;

    // return 0;
    // if(i==10) return 0;
  }



  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

#ifdef P0_CutFEM_LIMITER
double c0 = 0.77;//0.25;
R fun_levelSet(const R2 P, const int i) {
  return -0.5*P.x - P.y + 0.67*c0;
}
R fun_initial(const R2 P, int elementComp, int domain) {
  // return 1.;//
  if(domain == 1) return 0;
  // double xs = -0.3, ys = -0.3;
  double xs = 0.1, ys = 0.1;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0.;
  // return 1+0.5*sin(pi*(P.x+P.y));
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return 1+0.5*sin(pi*(P.x+P.y));
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(1,1), R2(1,1)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);

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

  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(1,1), R2(1,1)});
  double lambdaB = sqrt(10);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  Fun_h gh(Wh, fun_boundary, tn);
  Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , boundary
    , {1,4}          // label left boundary
  );

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
  std::ofstream outputData("output_CUTFEM_nx100_P1_delta1_fixed.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 100;
  int ny = 100;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P0);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(2) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;
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
// getchar();
  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.); interpolate(Wh, u0, fun_initial);
  Rn uh(u0);
  Rn uh_tild(u0);

  Fun_h fun_u0(Wh, u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_uh_tild(Wh, uh_tild);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "peiII.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
  //   writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  // }


  double qu0 = integral(Khi, fun_uh, 0);
  double min_u0 = 0.;//uh.min();
  double max_u0 = 1.;//uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    uh *= dt;
    uh += u0;
    limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0, macro);



    // std::map<int, double> u_mean;
    // uM = u0;
    // uh = u0;
    // limiter::extendToMacroP1(fun_uh, uM, u_mean, macro);
    // limiter::check_maximum_principle<Mesh>(u_mean, 0, 1);
    // The U_mean should satisfy the maximum principle
    // u0 = uh_tild;

    // THIRD ORDER RK
    // =================================================
    // {
    //   Rn u1(u0);
    //   Rn u2(Wh.NbDoF(), 0.);
    //   Rn u2_tild(Wh.NbDoF(), 0.);
    //   Rn u1_tild(Wh.NbDoF(), 0.);
    //   Fun_h fun_u1(Wh, u1);
    //   Fun_h fun_u2(Wh, u2);
    //
    //
    //   solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    //   u1 += dt * uh;
    //   limiter::CutFEM::limiter_Pei_P0(fun_u1, u1_tild, min_u0, max_u0, macro);
    //
    //   u2 += 3./4*u0 + 1./4*u1_tild;
    //   solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
    //   u2 += 1./4 * dt * uh;
    //   limiter::CutFEM::limiter_Pei_P0(fun_u2, u2_tild, min_u0, max_u0, macro);
    //
    //   solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
    //   uh *= 2./3 * dt;
    //   uh += 1./3*u0 + 2./3*u2_tild;
    //   limiter::CutFEM::limiter_Pei_P0(fun_uh, uh_tild, min_u0, max_u0, macro);
    //
    // }

    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = integral(Khi, fun_uh_tild, 0);
    double min_uh = uh.min();
    double max_uh = uh.max();
    // limiter::CutFEM::minmaxP0(fun_uh, min_uh, max_uh);
    double min_u1 = uh_tild.min();
    double max_u1 = uh_tild.max();
    // limiter::CutFEM::minmaxP0(fun_uh_tild, min_u1, max_u1);



    // if((i==24) && MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    //   writer.add(fun_uM, "macroExtend", 0, 1);
    // }



    if(MPIcf::IamMaster()) {
      Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
      writer.add(fun_uh, "uhNoLimiter", 0, 1);
      writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      // writer.add(fun_uM, "macroExtend", 0, 1);
    }
    getchar();

    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    if( min_u0 <= min_u1 && max_u1 <= max_u0) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;


      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }
      // return 0;
    }
    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
    //
    //   Paraview<Mesh> writer(Khi, "conservationP1_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    // }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(16) << min_u1 << "\t"
               << setprecision(16) << max_u1 << "\t"
               << setprecision(5) << std::endl;
    // getchar();
    // return 0;
    // if(i==10) return 0;
  }



  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

#ifdef CUTFEM_DISCONTINUOUS_SOLUTION
double c0 = 0.77;//0.25;
R fun_levelSet(const R2 P, const int i) {
  return -0.5*P.x - P.y + 0.67*c0;
}
R fun_initial(const R2 P, int elementComp, int domain) {
  // return 1.;//
  if(domain == 1) return 0;
  // double xs = -0.3, ys = -0.3;
  double xs = 0.1, ys = 0.1;
  double r = 0.3;
  double v = (P.x-xs)*(P.x-xs) + (P.y-ys)*(P.y-ys);
  if(v < r*r) return 1;//exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
  else return 0.;
  // return 1+0.5*sin(pi*(P.x+P.y));
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return 0.;//1+0.5*sin(pi*(P.x+P.y));
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  CutFEM_R2 beta({R2(1,1), R2(1,1)});
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);

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
    , INTEGRAL_INNER_EDGE_2D
  );


  problem.addBilinear(
     -innerProduct(beta*u*n, v)
     , Khi
     , INTEGRAL_BOUNDARY
     , {2, 3}          // label other boundary
  );
  problem.addBilinear(
    - innerProduct(u, beta*(0.5*v)*n)
    - innerProduct(u, lambdaB*0.5*v)
    , Khi
    , INTEGRAL_BOUNDARY
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

  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  CutFEM_R2 beta({R2(1,1), R2(1,1)});
  double lambdaB = sqrt(10);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  Fun_h gh(Wh, fun_boundary, tn);
  Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  problem.addLinear(
     - innerProduct(gx, beta*  (0.5*v)*n)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , INTEGRAL_BOUNDARY
    , {1,4}          // label left boundary
  );

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
  std::ofstream outputData("output_CUTFEM_nx100_P1_delta1_fixed.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 20;
  int ny = 20;
  Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P1dc);
  Space Eh(Th, DataFE<Mesh2>::P0);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 5 / sqrt(2) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;
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
// getchar();
  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.); interpolate(Wh, u0, fun_initial);
  Rn uh(u0);
  Rn uh_tild(u0);

  Fun_h fun_u0(Wh, u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_uh_tild(Wh, uh_tild);

  // Plot the macro elements
  {
    Paraview<Mesh> writer(Khi, "peiII.vtk");
    writer.add(fun_uh, "uh", 0, 1);
    writer.add(levelSet, "levelSet", 0, 1);
    writer.writeMacroInnerEdge (macro, 0, "pei_macro_inner_edge1.vtk");
    writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
    writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
    writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  }
  // return 0;

  double qu0 = integral(Khi, fun_uh, 0);
  double min_u0 = 0.;//uh.min();
  double max_u0 = 1.;//uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    // uh *= dt;
    // uh += u0;
    // // std::cout << uh << std::endl;
    // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0, macro);



    // std::map<int, double> u_mean;
    // uM = u0;
    // uh = u0;
    // limiter::extendToMacroP1(fun_uh, uM, u_mean, macro);
    // limiter::check_maximum_principle<Mesh>(u_mean, 0, 1);
    // The U_mean should satisfy the maximum principle
    // u0 = uh_tild;

    // THIRD ORDER RK
    // =================================================
    {
      Rn u1(u0);
      Rn u2(Wh.NbDoF(), 0.);
      Rn u2_tild(Wh.NbDoF(), 0.);
      Rn u1_tild(Wh.NbDoF(), 0.);
      Fun_h fun_u1(Wh, u1);
      Fun_h fun_u2(Wh, u2);

      Rn u_mean(Wh.get_nb_element());
      std::map<int, double>  map_mean_value;
      Rn uM(u0);
      Fun_h fun_uM(Wh, uM);
      double min_u0=0., max_u0=1.;

      // std::cout << u0 << std::endl;
      // limiter::CutFEM::check_mean_value(fun_u0, macro, 0., 1., u_mean);
      // limiter::CutFEM::limiter_Pei_P1(fun_u0, uM, min_u0, max_u0, macro);
      // limiter::CutFEM::extendToMacroP1(fun_u0, uM, map_mean_value, macro);
      limiter::CutFEM::minmaxP1(fun_uM, min_u0, max_u0);
      std::cout << min_u0 << "\t" << max_u0 << std::endl;


      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
      //   writer.add(fun_u0, "u0", 0, 1);
      //   writer.add(fun_uM, "macroExtend", 0, 1);
      // }
      // limiter::CutFEM::limiter_Pei_P1(fun_u0, uM, min_u0, max_u0, macro);
      // limiter::CutFEM::check_mean_value(fun_u0, macro, 0., 1., u_mean);

      // getchar();
      u0 = uM;

      std::cout << " compute u1 " << std::endl;
      solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
      u1 += dt * uh;
      uM = u1;


      limiter::CutFEM::check_mean_value(fun_u1, macro, 0.,1., u_mean);
      Fun_h fun_u_mean(Eh, u_mean);

      limiter::CutFEM::extendToMacroP1(fun_u1, uM, map_mean_value, macro);
      std::cout << "hey " << std::endl;
      Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
      writer.add(fun_u1, "u0", 0, 1);
      writer.add(fun_uM, "macroextended", 0, 1);
      writer.add(fun_u_mean, "uMean", 0, 1);

      uM = u1;
      limiter::CutFEM::limiter_Pei_P1(fun_u1, uM, min_u0, max_u0, macro);
      limiter::CutFEM::minmaxP1(fun_uM, min_u0, max_u0);

      std::cout << min_u0 << "\t" << max_u0 << std::endl;
      // if(MPIcf::IamMaster()) {
        // Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
        writer.add(fun_uM, "macroLimited", 0, 1);
      // }

      getchar();



      std::cout << " compute u2 " << std::endl;
      u2 += 3./4*u0 + 1./4*u1_tild;
      solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
      u2 += 1./4 * dt * uh;
      limiter::CutFEM::limiter_Pei_P1(fun_u2, u2_tild, min_u0, max_u0, macro);
      std::cout << " compute u3 " << std::endl;

      solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
      uh *= 2./3 * dt;
      uh += 1./3*u0 + 2./3*u2_tild;
      limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0, macro);


    }

    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = integral(Khi, fun_uh_tild, 0);
    double min_uh, max_uh;
    limiter::CutFEM::minmaxP1(fun_uh, min_uh, max_uh);
    double min_u1, max_u1;
    limiter::CutFEM::minmaxP1(fun_uh_tild, min_u1, max_u1);

    // if((i==24) && MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    //   writer.add(fun_uM, "macroExtend", 0, 1);
    // }



    // if(MPIcf::IamMaster()) {
    //   Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    //   // writer.add(fun_uM, "macroExtend", 0, 1);
    // }

    // if( (min_u0-Epsilon) < min_u1 && max_u1< (max_u0+Epsilon)) {
    if( min_u0 <= min_u1 && max_u1 <= max_u0) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;


      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }
      // return 0;
    }
    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
    //
    //   Paraview<Mesh> writer(Khi, "conservationP1_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_u1, "uhLimiter", 0, 1);
    // }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(16) << min_u1 << "\t"
               << setprecision(16) << max_u1 << "\t"
               << setprecision(5) << std::endl;
    getchar();
    // return 0;
    // if(i==10) return 0;
  }



  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

#ifdef CUTFEM_LIMITER_ACCURACY_TEST
double c0 = 0.77;//0.25;
R fun_levelSet(const R2 P, const int i) {
  return P.x*P.x + P.y*P.y - 1;
}
R fun_initial1 (const R2 P, int elementComp, int domain) {
  double r0 = 0.4;
  double a  = 0.3;
  return 0.5*(1 - tanh( ((P.x-r0)*(P.x-r0)+P.y*P.y) / (a*a) -1));
}
R fun_theta(const R2 P) {
  double r = P.norm();
  int s_y = util::fsign(P.y);
  if(P.x < -Epsilon) return (s_y*pi + atan(P.y/P.x));
  else if (P.x > Epsilon) return atan(P.y/(P.x));
  else return s_y*pi/2;
  // return atan(P.y/(P.x+Epsilon));
}

R fun_solution(const R2 P, int elementComp, int domain, double t) {
  double r = P.norm();
  double theta = fun_theta(P);
  R2 Q(r*cos(theta-2*pi*t), r*sin(theta-2*pi*t));
  return fun_initial1(Q,0,0);
}

R fun_initial (const R2 P, int elementComp, int domain) {
  // double r0 = 0.4;
  // double a  = 0.2;
  // return 0.5*(1 - tanh( ((P.x-r0)*(P.x-r0)+P.y*P.y) / (a*a) -1));
  return fun_solution(P,elementComp,domain,0);
}
R fun_boundary(const R2 P, int elementComp, double t) {
  return 0.;
}
R fun_velocity(const R2 P, int elementComp, int domain) {
  return (elementComp==0)? -2*pi*P.y : 2*pi*P.x;
}

void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh, const Fun_h& beta){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());
  CutFEM<Mesh2> problem(Wh);
  // CutFEM_R2 beta({R2(1,1), R2(1,1)});
  CutFEMParameter lambda(0, 1.);
  // CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
  double lambdaE = 1;
  const MeshParameter& h(Parameter::h);

  double lambdaB = sqrt(10);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;
  // Fun_h beta(Wh, fun_velocity);
  // Expression2 beta(fun_beta, 0, op_id);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);



  // BUILDING A
  // =====================================================
  problem.set_map(Ah);
  double tt = MPIcf::Wtime();

  problem.addBilinear(
    innerProduct(beta.expression(2) * u, grad(v))
    ,Khi
  );
  // matlab::Export(Ah, "mat1.dat");
  // std::cout << " Time assembly bilinear \t" << MPIcf::Wtime() - tt << std::endl;
  // getchar();

  // problem.addBilinear(
  //   - jump(innerProduct(beta*u,v*n))
  //   - innerProduct(jump(beta*u*n), jump(lambda*v))
  //   , interface
  // );
  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
      ,Khi
  );


  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(average(beta.expression(2)*u*n), jump(v))
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    , Khi
    , INTEGRAL_INNER_EDGE_2D
  );


  // problem.addBilinear(
  //    -innerProduct(beta.expression(2)*u*n, v)
  //    , interface
  //    // , Khi
  //    // , boundary
  //    // , {2, 3}          // label other boundary
  // );
  // problem.addBilinear(
  //   // - innerProduct(u, beta.expression(2)*(0.5*v)*n)
  //   - innerProduct(u, lambdaB*(0.5*v))
  //   , interface
  //   // , Khi
  //   // , INTEGRAL_BOUNDARY
  //   // , {1, 4}          // label left boundary
  // );
  // problem.addBilinear(
  //   - innerProduct(beta.expression(2)*u,v*n)
  //   - innerProduct(beta.expression(2)*u*n, lambdaB*v)
  //   , interface
  // );


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

  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  // CutFEM_R2 beta({R2(1,1), R2(1,1)});
  // double lambdaB = sqrt(10);

  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  // Fun_h gh(Wh, fun_boundary, tn);
  // Expression2 gx(gh, 0, op_id);

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  // mAh.addMatMul(u0, problem.rhs_);
  int N = problem.get_nb_dof();
  multiply(N,N,Ah,u0,problem.rhs_);

  // problem.addLinear(
  //    - innerProduct(gx, beta*  (0.5*v)*n)
  //    + innerProduct(gx, lambdaB*0.5*v)
  //    , Khi
  //    , boundary
  //   , {1,4}          // label left boundary
  // );

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
  std::ofstream outputData("output_smooth2_nx80_P0ooo.txt");

  // DEFINITION OF THE MESH and SPACE
  // ====================================================
  int nx = 40;//160;
  int ny = 40;//160;
  Mesh Th(nx, ny, -1.0075, -1.0075, 2.015, 2.015);   // [-1,1]*[-1,1]
  Space Vh(Th, DataFE<Mesh2>::P1dc);

  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./nx;
  double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 1;
  int niteration = tend / dt;
  dt = tend / niteration;
  double errSum = 0;

  std::cout << "Mesh size h = \t" << meshSize << std::endl;
  std::cout << "Time step dt = \t" << dt << std::endl;

  // DEFINITION OF THE LEVELSET
  // =====================================================
  Space Lh(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  Lagrange2 FE_beta(1);
  Space Uh(Th, FE_beta);

  // CONSTRUCTION INTERFACE AND CUTSPACE
  // =====================================================
  InterfaceLevelSet<Mesh> interface(Th, levelSet);
  ActiveMesh<Mesh> Khi(Th);
  Khi.truncate(interface, 1);
  CutSpace Wh(Khi, Vh);
  CutSpace Qh(Khi, Uh);
  Fun_h beta(Qh, fun_velocity);

  MacroElement<Mesh> macro(Khi, 0.5);

  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.); interpolate(Wh, u0, fun_initial);
  Rn uh(u0);
  Rn uh_tild(u0);

  Fun_h fun_u0(Wh, u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_uh_tild(Wh, uh_tild);

  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "test_accuracy_mesh.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge1.vtk");
  //   // writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
  //   // writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
  // }
  // return 0;

  double qu0 = integral(Khi, fun_uh, 0);
  double min_u0 = uh.min();
  double max_u0 = uh.max();

  // ASSEMBLY THE CONSTANT PART
  // ==================================================
  MatMap Ah, Mh;
  assembly(Wh, interface, Ah, Mh, beta);

  // RESOLUTION OF THE PROBLEM_MIXED_DARCY
  // ==================================================
  int ifig = 1;
  for(int i=0;i<niteration;++i) {
    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;

    // EULER METHOD
    // =================================================
    // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
    // uh *= dt;
    // uh += u0;
    // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0, macro);
    //


    // std::map<int, double> u_mean;
    // uM = u0;
    // uh = u0;
    // limiter::extendToMacroP1(fun_uh, uM, u_mean, macro);
    // limiter::check_maximum_principle<Mesh>(u_mean, 0, 1);
    // The U_mean should satisfy the maximum principle
    // u0 = uh_tild;

    // THIRD ORDER RK
    // =================================================
    {
      Rn u1(u0);
      Rn u2(Wh.NbDoF(), 0.);
      Rn u2_tild(Wh.NbDoF(), 0.);
      Rn u1_tild(Wh.NbDoF(), 0.);
      Fun_h fun_u1(Wh, u1);
      Fun_h fun_u2(Wh, u2);

      if(Wh.basisFctType == BasisFctType::P1dc) {
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild, min_u0, max_u0, macro);

        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild, min_u0, max_u0, macro);

        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild, min_u0, max_u0, macro);
      }
      else if (Wh.basisFctType == BasisFctType::P0){

        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::extendToMacroP0(fun_u1, u1_tild, macro);
        // u1_tild = u1;
        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::extendToMacroP0(fun_u2, u2_tild, macro);
        // u2_tild = u2;
        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::extendToMacroP0(fun_uh, uh_tild, macro);
        // uh_tild = uh;
      }
      else {
        assert(0);
      }
    }

    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double qu = integral(Khi, fun_uh_tild, 0);
    double min_uh, max_uh;
    limiter::CutFEM::findMinAndMaxValue(fun_uh, min_uh, max_uh);
    double min_u1, max_u1;
    limiter::CutFEM::findMinAndMaxValue(fun_uh_tild, min_u1, max_u1);

    if( min_u0 <= min_u1+Epsilon && max_u1 <= max_u0+Epsilon) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;

    }


    // COMPUTATION OF THE L2 ERROR
    // =================================================
    Rn usol(Wh.NbDoF(), 0.);
    interpolate(Wh, usol, fun_solution, tid);
    Rn uerr(usol);
    uerr -= uh_tild;
    Fun_h femErrh(Wh, uerr);
    Fun_h fun_ex (Wh, usol);

    Expression2 femErr(femErrh, 0, op_id);
    R errU = sqrt(integral(Khi, femErr*femErr));
    errSum += errU;

    // PLOT THE SOLUTION
    // ==================================================
    if(MPIcf::IamMaster() && i%100000 == 0 || i+1 == niteration) {
      // Fun_h fun_thet (Wh, fun_theta);
      Paraview<Mesh> writer(Khi, "test_accuracyP0_"+to_string(ifig++)+".vtk");
      writer.add(fun_uh, "uhNoLimiter", 0, 1);
      writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      writer.add(femErrh, "error"  , 0, 1);
      writer.add(fun_ex , "u_exact", 0, 1);
      // writer.add(fun_thet , "yx", 0, 1);

    }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << " || u-uex ||_2 = " << errU << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;


    if(i == 0) {
      outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU " << std::endl;

    }
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << errU << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(16) << min_u1 << "\t"
               << setprecision(16) << max_u1 << "\t"
               << setprecision(5) << std::endl;

  }


  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif

#ifdef CUTFEM_LIMITER_BURGER_TEST

// #define fluxFctX(I) 0.5*I*I
// #define fluxFctY(I) 0.5*I*I
// #define fluxFctN(I) (0.5*I*I*N.x + 0.5*I*I*N.y)

#define fluxFctX(I) burgerFlux(I)
#define fluxFctY(I) burgerFlux(I)
#define fluxFctN(I) burgerFlux(I,N)

// #define fluxFctX(I) 3*I
// #define fluxFctY(I) I
// #define fluxFctN(I) (3*I*N.x + I*N.y)

R fun_levelSet(const R2 P, const int i) {
  return P.y - 0.5*P.x - 1./4;
}
R fun_initial (const R2 P, int elementComp, int domain) {
  return sin(pi*(P.x+P.y));
}
// R fun_solution (const R2 P, int elementComp, int domain, double t) {
//   return sin(pi*(P.x+P.y - 4*t));
// }
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

    if(fabs(fm) < 1e-12) break;
    x0 = 0.5*(x1+x2);

    if(i == 50) {
      std::cout << " f(x0) = " << fm << std::endl;
      // assert(0);
      std::cout << " Burger Newton not converged" << std::endl;
      getchar();
    }
  }
  return sin(pi*x0);
}


void assembly(const Space& Wh, const Interface<Mesh>& interface, MatMap& Ah, MatMap& Mh){

  double t0 = MPIcf::Wtime();
  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());

  CutFEM<Mesh2> problem(Wh);
  CutFEMParameter lambda(0, 1.);
  CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
  const MeshParameter& h(Parameter::h);
  double Cstab = 1e-2;
  double Cstabt = 1e-2;
  Normal n;
  FunTest u(Wh,1), v(Wh,1);
  // CutFEM_R2 beta({R2(3,1), R2(3,1)});
  double lambdaB = sqrt(10);


  // BUILDING A
  // =====================================================
  problem.set_map(Ah);



  // F(u)_e = {B.u} - lambda_e/2 [u]
  problem.addBilinear(
    - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
    // - innerProduct(average(beta*u*n), jump(v))
    , Khi
    , INTEGRAL_INNER_EDGE_2D
  );

  // problem.addBilinear(
  //   - jump(innerProduct(beta*u,v*n))
  //   - innerProduct(jump(beta*u*n), jump(lambda*v))
  //   , interface
  // );

  //[ u \beta \cdot n]=3u_1-2u_2=0
  problem.addFaceStabilization(
      - innerProduct(jump(u), Cstab*jump(v))
      - innerProduct((h^2)*jump(grad(u)), Cstab*jump(grad(v)))
      // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
      ,Khi
  );

  // problem.addBilinear(
  //   - innerProduct(u, beta*(0.5*v)*n)
  //   - innerProduct(u, lambdaB*0.5*v)
  //   , Khi
  //   , INTEGRAL_BOUNDARY
  //   , {1, 4}          // label left boundary
  // );



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
      , Khi
  );

  std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space& Wh, const Interface<Mesh>& interface,Rn& u0,  Rn& uh, MatMap& Ah, MatMap& Mh, double tn) {

  double t0 = MPIcf::Wtime();

  const ActiveMesh<Mesh>& Khi(Wh.get_mesh());
  ProblemOption optionProblem;
  optionProblem.solver_name_ = "umfpack";
  optionProblem.clear_matrix_ = false;
  CutFEM<Mesh2> problem(Wh, optionProblem);

  // CutFEM_R2 beta({R2(1,1), R2(1,1)});
  double lambdaB = sqrt(10);
  CutFEMParameter lambda(0, 1.);

  Normal N;
  FunTest v(Wh,1);
  Fun_h fun_Un(Wh, u0);
  Expression   Un(fun_Un, 0, op_id);

  Fun_h gh(Wh, fun_solution, tn);
  Expression gx(gh, 0, op_id);

  CutFEM_R2 beta({R2(3,1), R2(3,1)});

  // MULTIPLYING A * u_0  => rhs
  // =====================================================
  MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
  mAh.addMatMul(u0, problem.rhs_);

  // matlab::Export(problem.rhs_, "rhs0.dat");
  // problem.rhs_ = 0.;
  // CONSTRUCT THE RHS
  double tt = MPIcf::Wtime();

  problem.addLinear(
        innerProduct(fluxFctX(Un)  , dx(v))
      + innerProduct(fluxFctY(Un)  , dy(v))
      , Khi
  );

  problem.addLinear(
    -innerProduct(average(fluxFctN(Un)), jump(v))
    , Khi
    , INTEGRAL_INNER_EDGE_2D
  );


  problem.addLinear(
    // - innerProduct(average(0.5*Un*Un*N.x), jump(v))
    // - innerProduct(average(0.5*Un*Un*N.y), jump(v))
    // - jump(innerProduct(fluxFctN(Un),v))
      // - jump(innerProduct(Un,v))
      // - innerProduct(average(fluxFctN(Un)), jump(v))
      // - innerProduct(0.5*jump(Un)  , jump(v))
      - jump(innerProduct(fluxFctN(Un),v))
      - innerProduct(jump(fluxFctN(Un)), jump(lambda*v))
    , interface
  );
//   matlab::Export(problem.rhs_, "rhs1.dat");
// std::cout << "hey " << std::endl;
//   getchar();
  // BOUNDARY CONDITION IN RHS
  problem.addLinear(
    - innerProduct(fluxFctN(Un), 0.5*v)
    - innerProduct(Un    , lambdaB*0.5*v)
    , Khi
    , INTEGRAL_BOUNDARY
    // , {1,4}          // boundary in
  );

  // problem.addLinear(
  //    - innerProduct(fluxFctN(Un), v)
  //   , Khi
  //   , INTEGRAL_BOUNDARY
  //   , {2,3}        // boundary out
  // );

  // matlab::Export(problem.rhs_, "rhs1.dat");
  // std::cout << "hey " << std::endl;
  // getchar();

  problem.addLinear(
    - innerProduct(fluxFctN(gx), 0.5*v)
     + innerProduct(gx, lambdaB*0.5*v)
     , Khi
     , INTEGRAL_BOUNDARY
    // , {1,4}          // label left boundary
  );
  // std::cout << "Time build rhs \t" << MPIcf::Wtime()-tt << std::endl;
  // getchar();
// SOLVING  (M+S)E'(t) = rhs
  // =====================================================
  tt = MPIcf::Wtime();
  problem.solve(Mh, problem.rhs_);
  // std::cout << "Time solver \t" << MPIcf::Wtime()-tt << std::endl;

  uh = problem.rhs_;
  problem.rhs_ = 0.;
}

int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // OUTPUT FILE
  // =====================================================
  std::ofstream outputData("output_burgerCutFEM_nx160_P0.txt");

  // DEFINITION OF THE MESH and SPACE
  // =====================================================
  int nx = 160;
  int ny = 160;
  Mesh2 Th(nx, ny, 0., 0., 2., 2.);
  Space Vh(Th, DataFE<Mesh2>::P0);

  // LagrangeDC2 FE_Flux(1);
  // Space Fh_0(Th, FE_Flux);


  // DEFINITION OF SPACE AND TIME PARAMETERS
  // =====================================================
  double tid = 0;
  double meshSize = 2./(nx+1);
  double dt = meshSize / 3 / sqrt(10) * 0.5;//0.3 * meshSize / 3 ;  // h /10
  double tend = 0.1;//1/pi;
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
  // CutSpace Fh(Khi, Fh_0);
  std::cout << "Number of dof = \t" << Wh.NbDoF() << std::endl;

  MacroElement<Mesh> macro(Khi, 0.5);

  // DECLARATION OF THE VECTOR CONTAINING THE solution
  // =====================================================
  Rn u0(Wh.NbDoF(), 0.);interpolate(Wh, u0, fun_initial);
  Rn uh(u0);
  Rn uh_tild(u0);

  Fun_h fun_u0(Wh, u0);
  Fun_h fun_uh(Wh, uh);
  Fun_h fun_uh_tild(Wh, uh_tild);


  // Plot the macro elements
  // {
  //   Paraview<Mesh> writer(Khi, "test_burger_mesh.vtk");
  //   writer.add(fun_uh, "uh", 0, 1);
  //   writer.add(levelSet, "levelSet", 0, 1);
  //   writer.writeMacroInnerEdge (macro, 0, "test_burger_macro_inner_edge1.vtk");
  //   writer.writeMacroOutterEdge(macro, 0, "test_burger_macro_outter_edge1.vtk");
  //   writer.writeMacroInnerEdge (macro, 1, "test_burger_macro_inner_edge2.vtk");
  //   writer.writeMacroOutterEdge(macro, 1, "test_burger_macro_outter_edge2.vtk");
  // }

  double qu0 = 0.;//integral(Khi, fun_uh, 0);
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
    double tt_iter = MPIcf::Wtime();

    std::cout << " ------------------------------ " << std::endl;
    std::cout << "Iteration " << i+1 << " / " << niteration
              << " \t time = " << tid << std::endl;


    // THIRD ORDER RK
    // =================================================
    {
      Rn u1(u0);
      Rn u2(Wh.NbDoF(), 0.);
      Rn u2_tild(Wh.NbDoF(), 0.);
      Rn u1_tild(Wh.NbDoF(), 0.);
      Fun_h fun_u1(Wh, u1);
      Fun_h fun_u2(Wh, u2);

      if(Wh.basisFctType == BasisFctType::P1dc) {
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild, min_u0, max_u0, macro);
        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild, min_u0, max_u0, macro);
        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild, min_u0, max_u0, macro);
      }
      else if (Wh.basisFctType == BasisFctType::P0){
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
        u2 += 3./4*u0 + 1./4*u1_tild;
        solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid+dt);
        u2 += 1./4 * dt * uh;
        limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
        solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid+0.5*dt);
        uh *= 2./3 * dt;
        uh += 1./3*u0 + 2./3*u2_tild;
        limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
      }

    }

    u0 = uh_tild;
    tid += dt;

    // COMPUTATION OF THE L2 ERROR
    // =================================================
    // Fun2_h femSolh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    double tt = MPIcf::Wtime();
    double qu = integral(Khi, fun_uh_tild, 0);
    double min_uh, max_uh;
    limiter::CutFEM::findMinAndMaxValue(fun_uh, min_uh, max_uh);
    double min_u1, max_u1;
    limiter::CutFEM::findMinAndMaxValue(fun_uh_tild, min_u1, max_u1);
    std::cout << " time min max \t" << MPIcf::Wtime() - tt << std::endl;


    if( min_u0 <= min_u1+Epsilon && max_u1 <= max_u0+Epsilon) {

      std::cout << " Maximum principle satified! " <<  std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;    }
    else{
      std::cout << " Maximum principle not satified! " << std::endl;
      std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
      std::cout << min_uh << " < u_{h,M}  < " << max_uh <<  std::endl;
      std::cout << min_u1 << " < u1_{h,M} < " << max_u1 <<  std::endl;

      // std::cout << min_u0 - min_u1 << std::endl;
      // std::cout << max_u0 - max_u1 << std::endl;


      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrincipleCutFEM_"+to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }
      // return 0;
    }


    // COMPUTATION OF THE L2 ERROR
    // =================================================
    Rn usol(Wh.NbDoF(), 0.);
    interpolate(Wh, usol, fun_solution, tid);
    Rn uerr(usol);
    uerr -= uh_tild;
    Fun_h femErrh(Wh, uerr);
    Fun_h fun_ex (Wh, usol);

    Expression2 femErr(femErrh, 0, op_id);
    R errU = sqrt(integral(Khi, femErr*femErr));
    errSum += errU;
    //
    // PLOT THE SOLUTION
    // ==================================================
    // if(MPIcf::IamMaster() && i%5 == 0 || i+1 == niteration) {
    //   // Fun_h fun_thet (Wh, fun_theta);
    //   Paraview<Mesh> writer(Khi, "test_burger_P1_"+to_string(ifig++)+".vtk");
    //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
    //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
    //   writer.add(femErrh, "error"  , 0, 1);
    //   writer.add(fun_ex , "u_exact", 0, 1);
    //   // writer.add(fun_thet , "yx", 0, 1);
    //
    // }

    // std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
    std::cout << " || u-uex ||_2 = " << errU << std::endl;
    std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu-qu0)
              << setprecision(6) << std::endl;


    if(i == 0) {
      outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU " << std::endl;

    }
    outputData << i << "\t"
               << tid << "\t"
               << setprecision(16) << errU << "\t"
               << setprecision(16) << qu << "\t"
               << setprecision(16) << fabs(qu-qu0) << "\t"
               << setprecision(16) << min_u1 << "\t"
               << setprecision(16) << max_u1 << "\t"
               << setprecision(5) << std::endl;

std::cout << " Iteration computation time \t" << MPIcf::Wtime() - tt_iter << std::endl;
  }

  std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
  std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
  return 0;
}


#endif
