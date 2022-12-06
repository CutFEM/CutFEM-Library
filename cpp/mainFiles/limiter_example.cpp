/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../tool.hpp"

using namespace globalVariable;

typedef std::map<std::pair<int, int>, R> MatMap;
typedef Mesh2 Mesh;
typedef FESpace2 Space;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;
typedef ExpressionFunFEM<Mesh2> Expression;

// #define FEM_SMOOTH_SOLUTION
// #define CUTFEM_SMOOTH_SOLUTION

#define CUTFEM_LIMITER_ACCURACY_TEST

// #define FEM_LIMITER_BURGER_TEST
// #define CUTFEM_LIMITER_BURGER_TEST

// #define FEM_DISCONTINUOUS_SOLUTION
// #define CUTFEM_DISCONTINUOUS_SOLUTION

#ifdef FEM_SMOOTH_SOLUTION

R fun_initial(double *P, int elementComp) {
   return 0.5 + 0.5 * sin(pi * (P[0] + P[1]));
}
R fun_boundary(double *P, int elementComp, double t) {
   return 0.5 + 0.5 * sin(pi * (P[0] + P[1] - 4 * t));
}
void assembly(const Space &Wh, MatMap &Ah, MatMap &Mh) {

   double t0 = getTime();
   const Mesh &Khi(Wh.Th);
   FEM<Mesh2> problem(Wh);
   R2 beta(3, 1);
   double lambdaE = sqrt(10);
   const MeshParameter &h(Parameter::h);
   double lambdaB = sqrt(10);

   double Cstab  = 1e-2;
   double Cstabt = 1e-2;

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);
   problem.addBilinear(innerProduct(beta * u, grad(v)), Khi);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(average(beta * u * n), jump(v)) -
                           innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                       Khi, INTEGRAL_INNER_EDGE_2D);

   problem.addBilinear(-innerProduct(beta * u * n, v), Khi, INTEGRAL_BOUNDARY,
                       {2, 3} // label other boundary
   );
   problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) -
                           innerProduct(u, lambdaB * 0.5 * v),
                       Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // BUILDING (M + S)
   // =====================================================
   problem.set_map(Mh);
   problem.addBilinear(innerProduct(u, v), Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Rn &u0, Rn &uh, MatMap &Ah,
                   MatMap &Mh, double tn) {

   double t0 = getTime();

   const Mesh &Khi(Wh.Th);

   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   FEM<Mesh2> problem(Wh, optionProblem);
   R2 beta(3, 1);
   double lambdaB = sqrt(10);
   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   Fun_h gh(Wh, fun_boundary, tn);
   Expression2 gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
   // mAh.addMatMul(u0, problem.rhs_);
   int N = problem.get_nb_dof();
   multiply(N, N, Ah, u0, problem.rhs_);

   problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   problem.solve(Mh, problem.rhs_);

   uh = problem.rhs_;
}
int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_testFEM_P0.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 40;
   int ny = 40;
   Mesh Kh(nx, ny, -1., -1., 2., 2.); // [-1,1]*[-1,1]
   Space Vh(Kh, DataFE<Mesh2>::P0);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / nx;
   double dt = meshSize / 5 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend    = 0.1;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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

   double qu0    = 0.; // integral(Khi, fun_uh, 0);
   double min_u0 = 0.; // uh.min();
   double max_u0 = 1.; // uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Vh, Ah, Mh);

   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
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
      u2 += 3. / 4 * u0 + 1. / 4 * u1;
      solve_problem(Vh, u1, uh, Ah, Mh, tid + dt);
      u2 += 1. / 4 * dt * uh;
      // limiter::CutFEM::extendToMacroP0(fun_u2, u2_tild, macro);
      // u2_tild = u2;
      solve_problem(Vh, u2, uh, Ah, Mh, tid + 0.5 * dt);
      uh *= 2. / 3 * dt;
      uh += 1. / 3 * u0 + 2. / 3 * u2;
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
      double qu = 0; // integral(Khi, fun_u1, 0);
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
      //   Paraview<Mesh> writer(Kh,
      //   "maxPrincipleFEM_"+std::to_string(ifig++)+".vtk"); writer.add(fun_uh,
      //   "uhNoLimiter", 0, 1); writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      // }
      // if( (min_u0-Epsilon) < min_u1 && max_u1<
      // (max_u0+Epsilon)) {
      if (min_u0 <= min_uh_tild + Epsilon && max_uh_tild + Epsilon <= max_u0) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h}  < " << max_uh << std::endl;
         std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild
                   << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h}  < " << max_uh << std::endl;
         std::cout << min_uh_tild << " < tild{u_{h}} < " << max_uh_tild
                   << std::endl;

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
      Fun_h fun_ex(Vh, usol);

      Expression2 femErr(femErrh, 0, op_id);
      R errU = sqrt(integral(Kh, femErr * femErr));
      errSum += errU;

      // PLOT THE SOLUTION
      // ==================================================
      if (MPIcf::IamMaster() && i % 10 == 0 || i + 1 == niteration) {

         Paraview<Mesh> writer(Kh, "smooth_FEM_P0_" + std::to_string(ifig++) +
                                       ".vtk");
         writer.add(fun_uh, "uhNoLimiter", 0, 1);
         writer.add(fun_u1, "uhLimiter", 0, 1);
         writer.add(fun_ex, "u_exact", 0, 1);
      }

      std::cout << " || u-uex ||_2 = " << errU << std::endl;
      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU
                 << "\t" << std::setprecision(16) << qu << "\t"
                 << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_uh_tild << "\t"
                 << std::setprecision(16) << max_uh_tild << "\t"
                 << std::setprecision(5) << std::endl;

      // return 0;
      // if(i==10) return 0;
   }

   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration
             << std::endl;
   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef CUTFEM_SMOOTH_SOLUTION
// double c0 = 0.77;//0.25;
// R fun_levelSet(double* P, const int i) {
//   return -0.5*P[0] - P[1] + 0.67*c0;
// }
// R fun_initial(double* P, int elementComp, int domain) {
//   return 0.5+0.5*sin(pi*(P[0]+P[1]));
// }
// R fun_boundary(double* P, int elementComp, int domain, double t) {
//     return 0.5+0.5*sin(pi*(P[0] + P[1]- 4*t));
// }
// R fun_solution(double* P, int elementComp, int domain, double t) {
//     return 0.5+0.5*sin(pi*(P[0] + P[1]- 4*t));
// }
R fun_levelSet(double *P, const int i) {
   // return P[1] - 0.5*P[0] - 1./4;
   return -(P[1] + 0.5 * P[0] - 0.73);
}

R fun_initial(double *P, int elementComp, int domain) {
   return sin(pi * (P[0] + P[1]));
}
R fun_solution(double *P, int elementComp, int domain, double t) {
   return sin(pi * (P[0] + P[1] - 4 * t));
}
R fun_boundary(double *P, int elementComp, int domain, double t) {
   return sin(pi * (P[0] + P[1] - 4 * t));
}

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah,
              MatMap &Mh) {

   double t0 = getTime();
   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
   CutFEM<Mesh2> problem(Wh);
   CutFEM_R2 beta({R2(3, 1), R2(3, 1)});
   CutFEMParameter lambda(0, 1.);
   CutFEMParameter lambdaE(sqrt(10), sqrt(10)); // max B = sqrt(10)/sqrt(5)
   const MeshParameter &h(Parameter::h);

   double lambdaB = sqrt(10);
   double Cstab   = 1e-2;
   double Cstabt  = 1e-2;
   // Fun_h beta(Wh, fun_velocity);
   // Expression2 beta(fun_beta, 0, op_id);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);
   problem.addBilinear(innerProduct(beta * u, grad(v)), Khi);

   problem.addBilinear(-jump(innerProduct(beta * u, v * n)) -
                           innerProduct(jump(beta * u * n), jump(lambda * v)),
                       interface);
   //[ u \beta \cdot n]=3u_1-2u_2=0
   problem.addFaceStabilization(
       -innerProduct(jump(u), Cstab * jump(v)) -
           innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
       // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(average(beta * u * n), jump(v)) -
                           innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                       Khi, INTEGRAL_INNER_EDGE_2D);

   problem.addBilinear(-innerProduct(beta * u * n, v), Khi, INTEGRAL_BOUNDARY,
                       {2, 3} // label other boundary
   );
   problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) -
                           innerProduct(u, lambdaB * 0.5 * v),
                       Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // BUILDING (M + S)
   // =====================================================
   problem.set_map(Mh);
   problem.addBilinear(innerProduct(u, v), Khi);
   problem.addFaceStabilization(
       innerProduct(h * jump(u), Cstab * jump(v)) +
           innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v)))
       // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Interface<Mesh> &interface,
                   const Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh, double tn) {

   double t0 = getTime();

   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   CutFEM<Mesh2> problem(Wh, optionProblem);

   CutFEM_R2 beta({R2(3, 1), R2(3, 1)});
   double lambdaB = sqrt(10);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   Fun_h gh(Wh, fun_boundary, tn);
   Expression2 gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
   // mAh.addMatMul(u0, problem.rhs_);
   int N = problem.get_nb_dof();
   multiply(N, N, Ah, u0, problem.rhs_);

   problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   problem.solve(Mh, problem.rhs_);

   uh = problem.rhs_;
}

int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_testP1.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 20;
   int ny = 20;
   // Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
   Mesh2 Th(nx, ny, 0., 0., 2., 2.);

   Space Vh(Th, DataFE<Mesh2>::P1dc);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / (nx + 1);
   double dt = meshSize / 3 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend    = 0.1;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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
   Rn u0(Wh.NbDoF(), 0.);
   interpolate(Wh, u0, fun_initial);
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

   double qu0    = integral(Khi, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Wh, interface, Ah, Mh);
   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
                << " \t time = " << tid << std::endl;

      // EULER METHOD
      // =================================================
      // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
      // uh *= dt;
      // uh += u0;
      // // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0,
      // macro); limiter::CutFEM::extendToMacroP0(fun_uh, uh_tild, macro);

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

         if (Wh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild,
                                                         min_u0, max_u0, macro);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild,
                                                         min_u0, max_u0, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild,
                                                         min_u0, max_u0, macro);
         } else if (Wh.basisFctType == BasisFctType::P0) {

            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
         } else {
            assert(0);
         }
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      // Fun2_h femSolh(Wh, uh);
      // Expression2 femSol(femSolh, 0, op_id);
      double qu             = integral(Khi, fun_uh_tild, 0);
      // double min_uh, max_uh;
      // limiter::CutFEM::findMinAndMaxValue(fun_uh, min_uh, max_uh);
      // double min_u1, max_u1;
      // limiter::CutFEM::findMinAndMaxValue(fun_uh_tild, min_u1, max_u1);
      auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
      auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);
      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi,
      //   "smoothSol_CutFEM_"+std::to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1); writer.add(fun_uh_tild,
      //   "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }

      // if( (min_u0-Epsilon) < min_u1 && max_u1<
      // (max_u0+Epsilon)) {
      if (min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0 + Epsilon) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

         std::cout << min_u0 - min_u1 << std::endl;
         std::cout << max_u0 - max_u1 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi,
         //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
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
      Fun_h fun_ex(Wh, usol);

      Expression2 femErr(femErrh, 0, op_id);
      R errU = sqrt(integral(Khi, femErr * femErr));
      errSum += errU;

      // PLOT THE SOLUTION
      // ==================================================
      if (MPIcf::IamMaster() && i % 5 == 0 || i + 1 == niteration) {

         Paraview<Mesh> writer(Khi, "smoothSolCutFEM_" +
                                        std::to_string(ifig++) + ".vtk");
         writer.add(fun_uh, "uhNoLimiter", 0, 1);
         writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         writer.add(femErrh, "error", 0, 1);
         writer.add(fun_ex, "u_exact", 0, 1);
      }

      // std::cout << std::setprecision(16) << "q(u) = " << qu << std::endl;
      std::cout << " || u-uex ||_2 = " << errU << std::endl;
      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU
                 << "\t" << std::setprecision(16) << qu << "\t"
                 << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(5) << std::endl;
   }

   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration
             << std::endl;
   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef CUTFEM_LIMITER_ACCURACY_TEST
double c0 = 0.77; // 0.25;
R fun_levelSet(double *P, const int i) { return P[0] * P[0] + P[1] * P[1] - 1; }
R fun_initial1(double *P, int elementComp, int domain) {
   double r0 = 0.4;
   double a  = 0.3;
   return 0.5 *
          (1 - tanh(((P[0] - r0) * (P[0] - r0) + P[1] * P[1]) / (a * a) - 1));
}
R fun_theta(double *P) {
   double r = Norme2((R2)(P)); // P.norm();
   int s_y  = util::fsign(P[1]);
   if (P[0] < -Epsilon)
      return (s_y * pi + atan(P[1] / P[0]));
   else if (P[0] > Epsilon)
      return atan(P[1] / (P[0]));
   else
      return s_y * pi / 2;
   // return atan(P[1]/(P[0]+Epsilon));
}

R fun_solution(double *P, int elementComp, int domain, double t) {
   // double r = P.norm();
   double r     = Norme2((R2)(P));
   double theta = fun_theta(P);
   R2 Q(r * cos(theta - 2 * pi * t), r * sin(theta - 2 * pi * t));
   return fun_initial1(Q, 0, 0);
}

R fun_initial(double *P, int elementComp, int domain) {
   // double r0 = 0.4;
   // double a  = 0.2;
   // return 0.5*(1 - tanh( ((P[0]-r0)*(P[0]-r0)+P[1]*P[1]) / (a*a) -1));
   return fun_solution(P, elementComp, domain, 0);
}
R fun_boundary(double *P, int elementComp, double t) { return 0.; }
R fun_velocity(double *P, int elementComp, int domain) {
   return (elementComp == 0) ? -2 * pi * P[1] : 2 * pi * P[0];
}

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah,
              MatMap &Mh, const Fun_h &beta) {

   double t0 = getTime();
   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
   CutFEM<Mesh2> problem(Wh);
   CutFEMParameter lambda(0, 1.);
   double lambdaE = 1;
   const MeshParameter &h(Parameter::h);

   double lambdaB = sqrt(10);
   double Cstab   = 1e-2;
   double Cstabt  = 1e-2;
   // Fun_h beta(Wh, fun_velocity);
   // Expression2 beta(fun_beta, 0, op_id);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);
   double tt = getTime();

   problem.addBilinear(innerProduct(beta.exprList(2) * u, grad(v)), Khi);

   // problem.addBilinear(
   //   - jump(innerProduct(beta*u,v*n))
   //   - innerProduct(jump(beta*u*n), jump(lambda*v))
   //   , interface
   // );
   //[ u \beta \cdot n]=3u_1-2u_2=0
   problem.addFaceStabilization(
       -innerProduct(jump(u), Cstab * jump(v)) -
           innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
       // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(
       -innerProduct(average(beta.exprList(2) * u * n), jump(v)) -
           innerProduct(0.5 * lambdaE * jump(u), jump(v)),
       Khi, INTEGRAL_INNER_EDGE_2D);

   // problem.addBilinear(
   //    -innerProduct(beta.expression(2)*u*n, v)
   //    , interface
   //    // , Khi
   //    // , boundary
   //    // , {2, 3}          // label other boundary
   // );
   // problem.addBilinear(
   //   - innerProduct(u, beta.expression(2)*(0.5*v)*n)
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
   problem.addBilinear(innerProduct(u, v), Khi);
   problem.addFaceStabilization(
       innerProduct(h * jump(u), Cstab * jump(v)) +
           innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v)))
       // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Interface<Mesh> &interface,
                   const Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh, double tn) {

   double t0 = getTime();

   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   CutFEM<Mesh2> problem(Wh, optionProblem);

   // CutFEM_R2 beta({R2(1,1), R2(1,1)});
   // double lambdaB = sqrt(10);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   // Fun_h gh(Wh, fun_boundary, tn);
   // Expression2 gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
   // mAh.addMatMul(u0, problem.rhs_);
   int N = problem.get_nb_dof();
   multiply(N, N, Ah, u0, problem.rhs_);

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

int main(int argc, char **argv) {

#ifdef USE_MPI
   MPIcf cfMPI(argc, argv);
#endif
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_smooth2_nx80_P0ooo.txt");

   // DEFINITION OF THE MESH and SPACE
   // ====================================================
   int nx = 20;                                     // 160;
   int ny = 20;                                     // 160;
   Mesh Th(nx, ny, -1.0075, -1.0075, 2.015, 2.015); // [-1,1]*[-1,1]
   Space Vh(Th, DataFE<Mesh2>::P1dc);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / nx;
   double dt = meshSize / 5 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend    = 0.1;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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

   MacroElement<Mesh> macro(Khi, 1.);

   // DECLARATION OF THE VECTOR CONTAINING THE solution
   // =====================================================
   Rn u0(Wh.NbDoF(), 0.);
   interpolate(Wh, u0, fun_initial);
   Rn uh(u0);
   Rn uh_tild(u0);

   Fun_h fun_u0(Wh, u0);
   Fun_h fun_uh(Wh, uh);
   Fun_h fun_uh_tild(Wh, uh_tild);

   // Plot the macro elements
   {
      Paraview<Mesh> writer(Khi, "test_accuracy_mesh.vtk");
      writer.add(fun_uh, "uh", 0, 1);
      // writer.add(levelSet, "levelSet", 0, 1);
      writer.writeMacroInnerEdge(macro, 0, "pei_macro_inner_edge1.vtk");
      writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
      // writer.writeMacroInnerEdge (macro, 1, "pei_macro_inner_edge2.vtk");
      // writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
   }
   // return 0;

   double qu0    = integral(Khi, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Wh, interface, Ah, Mh, beta);

   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
                << " \t time = " << tid << std::endl;

      // EULER METHOD
      // =================================================
      // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
      // uh *= dt;
      // uh += u0;
      // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0,
      // macro);
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

         if (Wh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild,
                                                         min_u0, max_u0, macro);

            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild,
                                                         min_u0, max_u0, macro);

            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild,
                                                         min_u0, max_u0, macro);
         } else if (Wh.basisFctType == BasisFctType::P0) {

            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
            // u1_tild = u1;
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
            // u2_tild = u2;
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
            // uh_tild = uh;
         } else {
            assert(0);
         }
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      // Fun2_h femSolh(Wh, uh);
      // Expression2 femSol(femSolh, 0, op_id);
      double qu             = integral(Khi, fun_uh_tild, 0);
      auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
      auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);

      if (min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0 + Epsilon) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      }

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      Rn usol(Wh.NbDoF(), 0.);
      interpolate(Wh, usol, fun_solution, tid);
      Rn uerr(usol);
      uerr -= uh_tild;
      Fun_h femErrh(Wh, uerr);
      Fun_h fun_ex(Wh, usol);

      // Expression2 femErr(femErrh, 0, op_id);
      R errU = sqrt(integral(Khi, femErrh.expr() * femErrh.expr()));
      errSum += errU;

      // PLOT THE SOLUTION
      // ==================================================
#ifdef USE_MPI
      if (MPIcf::IamMaster() && i % 1 == 0 || i + 1 == niteration) {
#else
      if (i % 1 == 0 || i + 1 == niteration) {
#endif
         // Fun_h fun_thet (Wh, fun_theta);
         Paraview<Mesh> writer(Khi, "test_accuracyP0_" +
                                        std::to_string(ifig++) + ".vtk");
         writer.add(fun_uh, "uhNoLimiter", 0, 1);
         writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         writer.add(femErrh, "error", 0, 1);
         writer.add(fun_ex, "u_exact", 0, 1);
         // writer.add(fun_thet , "yx", 0, 1);
      }
      if (i == 10)
         return 0;
      // std::cout << std::setprecision(16) << "q(u) = " << qu << std::endl;
      std::cout << " || u-uex ||_2 = " << errU << std::endl;
      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU
                 << "\t" << std::setprecision(16) << qu << "\t"
                 << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_u1 << "\t"
                 << std::setprecision(16) << max_u1 << "\t"
                 << std::setprecision(5) << std::endl;
   }

   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration
             << std::endl;
   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef FEM_LIMITER_BURGER_TEST
#define fluxFctX(I) burgerFlux(I)
#define fluxFctY(I) burgerFlux(I)
#define fluxFctN(I) burgerFlux(I, N)

R fun_initial(double *P, int elementComp) {
   return 1 + 0.5 * sin(pi * (P[0] + P[1]));
}
// R fun_solution(double* P, int elementComp, int domain, double t) {
//   double C = P[0]+P[1];
//   double x1 = 0;
//   double x2 = 4;
//   double x0 = 0.5*(x1+x2);
//   double f1 = x1 + 2.*sin(pi*x1)*t - C;
//   if(fabs(f1) < 1e-10) {
//     return sin(pi*x1);
//   }
//   double f2 = x2 + 2.*sin(pi*x2)*t - C;
//   if(fabs(f2) < 1e-10) {
//     return sin(pi*x2);
//   }
//
//   for(int i=0;i<=50;++i) {
//     double fa = x1 + 2.*sin(pi*x1)*t - C;
//     double fm = x0 + 2.*sin(pi*x0)*t - C;
//
//     if(fa*fm < 0) x2 = x0;
//     else          x1 = x0;
//
//     if(fabs(fm) < 1e-12) break;
//     x0 = 0.5*(x1+x2);
//
//     if(i == 50) {
//       std::cout << " f(x0) = " << fm << std::endl;
//       // assert(0);
//       std::cout << " Burger Newton not converged" << std::endl;
//       getchar();
//     }
//   }
//   return sin(pi*x0);
// }
R fun_scalar(double x) { return 1 + 0.5 * sin(pi * x); }
R fun_solution(double *P, int elementComp, int domain, double t) {
   double x0 = P[0] + P[1];
   double xy = P[0] + P[1];
   for (int i = 0; i <= 50; ++i) {

      double x1 = xy - 2 * t * fun_scalar(x0);
      if (fabs(x1 - x0) < 1e-12)
         return fun_scalar(x1);
      x0 = x1;
      if (i == 50) {
         std::cout << " Burger Newton not converged" << std::endl;
         getchar();
      }
   }
}

void assembly(const Space &Wh, MatMap &Ah, MatMap &Mh) {

   double t0 = getTime();
   const Mesh &Khi(Wh.Th);

   FEM<Mesh2> problem(Wh);
   // CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
   const MeshParameter &h(Parameter::h);
   double Cstab  = 1e-2;
   double Cstabt = 1e-2;
   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   // CutFEM_R2 beta({R2(3,1), R2(3,1)});
   double lambdaB = sqrt(10);
   double lambdaE = sqrt(10);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(0.5 * lambdaE * jump(u), jump(v))
                       // - innerProduct(average(beta*u*n), jump(v))
                       ,
                       Khi, INTEGRAL_INNER_EDGE_2D);

   // BUILDING (M + S)
   // =====================================================
   problem.set_map(Mh);
   problem.addBilinear(innerProduct(u, v), Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh,
                   double tn) {

   double t0 = getTime();

   const Mesh &Khi(Wh.Th);
   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   FEM<Mesh2> problem(Wh, optionProblem);

   // CutFEM_R2 beta({R2(1,1), R2(1,1)});
   double lambdaB = sqrt(10);
   Normal N;
   FunTest v(Wh, 1);
   Fun_h fun_Un(Wh, u0);
   Expression Un(fun_Un, 0, op_id);

   Fun_h gh(Wh, fun_solution, tn);
   Expression gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
   mAh.addMatMul(u0, problem.rhs_);

   // CONSTRUCT THE RHS
   double tt = getTime();

   problem.addLinear(innerProduct(fluxFctX(Un), dx(v)) +
                         innerProduct(fluxFctY(Un), dy(v)),
                     Khi);

   problem.addLinear(-innerProduct(average(fluxFctN(Un)), jump(v)), Khi,
                     INTEGRAL_INNER_EDGE_2D);

   // BOUNDARY CONDITION IN RHS
   problem.addLinear(-innerProduct(fluxFctN(Un), 0.5 * v) -
                         innerProduct(Un, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY
                     // , {1,4}          // boundary in
   );

   // problem.addLinear(
   //    - innerProduct(fluxFctN(Un), v)
   //   , Khi
   //   , INTEGRAL_BOUNDARY
   //   , {2,3}        // boundary out
   // );

   problem.addLinear(-innerProduct(fluxFctN(gx), 0.5 * v) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY
                     // , {1,4}          // label left boundary
   );

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   tt = getTime();
   problem.solve(Mh, problem.rhs_);
   // std::cout << "Time solver \t" << getTime()-tt << std::endl;

   uh           = problem.rhs_;
   problem.rhs_ = 0.;
}

int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_burgerFEM_nx40_P1.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 40;
   int ny = 40;
   Mesh2 Th(nx, ny, 0., 0., 2., 2.);
   Space Vh(Th, DataFE<Mesh2>::P1dc);

   // LagrangeDC2 FE_Flux(1);
   // Space Fh_0(Th, FE_Flux);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / (nx + 1);
   double dt = meshSize / 3 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend    = 0.1;                      // 1/pi;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

   std::cout << "Mesh size h = \t" << meshSize << std::endl;
   std::cout << "Time step dt = \t" << dt << std::endl;

   // DEFINITION OF THE LEVELSET
   // =====================================================
   Space Lh(Th, DataFE<Mesh2>::P1);

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
   //   Paraview<Mesh> writer(Khi, "test_burger_mesh.vtk");
   //   writer.add(fun_uh, "uh", 0, 1);
   //   writer.add(levelSet, "levelSet", 0, 1);
   //   writer.writeMacroInnerEdge (macro, 0,
   //   "test_burger_macro_inner_edge1.vtk"); writer.writeMacroOutterEdge(macro,
   //   0, "test_burger_macro_outter_edge1.vtk"); writer.writeMacroInnerEdge
   //   (macro, 1, "test_burger_macro_inner_edge2.vtk");
   //   writer.writeMacroOutterEdge(macro, 1,
   //   "test_burger_macro_outter_edge2.vtk");
   // }

   double qu0    = 0.; // integral(Khi, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Vh, Ah, Mh);
   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      double tt_iter = getTime();

      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
                << " \t time = " << tid << std::endl;

      // THIRD ORDER RK
      // =================================================
      {
         Rn u1(u0);
         Rn u2(Vh.NbDoF(), 0.);
         Fun_h fun_u1(Vh, u1);
         Fun_h fun_u2(Vh, u2);

         if (Vh.basisFctType == BasisFctType::P1dc) {
            Rn u2_tild(Vh.NbDoF(), 0.);
            Rn u1_tild(Vh.NbDoF(), 0.);
            solve_problem(Vh, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::FEM::applyBoundPreservingLimiter(fun_u1, u1_tild, min_u0,
                                                      max_u0);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Vh, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::FEM::applyBoundPreservingLimiter(fun_u2, u2_tild, min_u0,
                                                      max_u0);
            solve_problem(Vh, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::FEM::applyBoundPreservingLimiter(fun_uh, uh_tild, min_u0,
                                                      max_u0);

         } else if (Vh.basisFctType == BasisFctType::P0) {
            solve_problem(Vh, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            u2 += 3. / 4 * u0 + 1. / 4 * u1;
            solve_problem(Vh, u1, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            solve_problem(Vh, u2, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2;
            uh_tild = uh;
         }
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      // Fun2_h femSolh(Wh, uh);
      // Expression2 femSol(femSolh, 0, op_id);
      double tt = getTime();
      double qu = integral(Th, fun_uh_tild, 0);
      double min_uh, max_uh;
      limiter::FEM::findMinAndMaxValue(fun_uh, min_uh, max_uh);
      double min_u1, max_u1;
      limiter::FEM::findMinAndMaxValue(fun_uh_tild, min_u1, max_u1);
      std::cout << " time min max \t" << getTime() - tt << std::endl;

      if (min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0 + Epsilon) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

         // std::cout << min_u0 - min_u1 << std::endl;
         // std::cout << max_u0 - max_u1 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi,
         //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
         //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
         //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         //   // writer.add(fun_uM, "macroExtend", 0, 1);
         // }
         // return 0;
      }

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      Rn usol(Vh.NbDoF(), 0.);
      interpolate(Vh, usol, fun_solution, tid);
      Rn uerr(usol);
      uerr -= uh_tild;
      Fun_h femErrh(Vh, uerr);
      Fun_h fun_ex(Vh, usol);

      Expression2 femErr(femErrh, 0, op_id);
      R errU = sqrt(integral(Th, femErr * femErr));
      errSum += errU;
      //
      // PLOT THE SOLUTION
      // ==================================================
      // if(MPIcf::IamMaster() && i%5 == 0 || i+1 == niteration) {
      //   // Fun_h fun_thet (Wh, fun_theta);
      //   Paraview<Mesh> writer(Th,
      //   "test_burgerFEM_P01_"+std::to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   // writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      //   writer.add(femErrh, "error"  , 0, 1);
      //   writer.add(fun_ex , "u_exact", 0, 1);
      //   // writer.add(fun_thet , "yx", 0, 1);
      //
      // }

      // std::cout << std::setprecision(16) << "q(u) = " << qu << std::endl;
      std::cout << " || u-uex ||_2 = " << errU << std::endl;
      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU
                 << "\t" << std::setprecision(16) << qu << "\t"
                 << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_u1 << "\t"
                 << std::setprecision(16) << max_u1 << "\t"
                 << std::setprecision(5) << std::endl;

      std::cout << " Iteration computation time \t" << getTime() - tt_iter
                << std::endl;
   }

   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration
             << std::endl;
   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef CUTFEM_LIMITER_BURGER_TEST
#define fluxFctX(I) burgerFlux(I)
#define fluxFctY(I) burgerFlux(I)
#define fluxFctN(I) burgerFlux(I, N)

// #define fluxFctX(I) (3*I)
// #define fluxFctY(I) (I)
// #define fluxFctN(I) (3.*I*N.x + I*N.y)

R fun_levelSet(double *P, const int i) { return -(P[1] - 0.5 * P[0] - 1. / 4); }
R fun_initial(double *P, int elementComp, int domain) {
   return 1 + 0.5 * sin(pi * (P[0] + P[1]));
}
R fun_scalar(double x) { return 1 + 0.5 * sin(pi * x); }
R fun_solution(double *P, int elementComp, int domain, double t) {
   double x0 = P[0] + P[1];
   double xy = x0;
   for (int i = 0; i <= 50; ++i) {

      double x1 = xy - 2 * t * fun_scalar(x0);
      if (fabs(x1 - x0) < 1e-10)
         return fun_scalar(x1);
      x0 = x1;
      if (i == 100) {
         std::cout << " Burger Newton not converged \t"
                   << fabs(xy - 2 * t * fun_scalar(x0) - x0) << std::endl;

         getchar();
      }
   }
}

// R fun_initial (double* P, int elementComp, int domain) {
//   return sin(pi*(P[0]+P[1]));
// }
// R fun_solution (double* P, int elementComp, int domain, double t) {
//   return sin(pi*(P[0]+P[1] - 4*t));
// }

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah,
              MatMap &Mh) {

   double t0 = getTime();
   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

   CutFEM<Mesh2> problem(Wh);
   CutFEMParameter lambda(0, 1.);
   CutFEMParameter lambdaE(sqrt(10), sqrt(10)); // max B = sqrt(10)/sqrt(5)
   // CutFEMParameter lambdaE(3., 3.);  // max B = sqrt(10)/sqrt(5)

   const MeshParameter &h(Parameter::h);
   double Cstab  = 1e-2;
   double Cstabt = 1e-2;
   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   // double lambdaB = 2.;//sqrt(10);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(0.5 * lambdaE * jump(u), jump(v))
                       // - innerProduct(average(beta*u*n), jump(v))
                       ,
                       Khi, INTEGRAL_INNER_EDGE_2D);

   // problem.addBilinear(
   //   - jump(innerProduct(beta*u,v*n))
   //   - innerProduct(jump(beta*u*n), jump(lambda*v))
   //   , interface
   // );

   //[ u \beta \cdot n]=3u_1-2u_2=0
   problem.addFaceStabilization(
       -innerProduct(jump(u), Cstab * jump(v)) -
           innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
       // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

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
   problem.addBilinear(innerProduct(u, v), Khi);
   problem.addFaceStabilization(
       innerProduct(h * jump(u), Cstab * jump(v)) +
           innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v))),
       Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Interface<Mesh> &interface, Rn &u0,
                   Rn &uh, MatMap &Ah, MatMap &Mh, double tn) {

   double t0 = getTime();

   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   CutFEM<Mesh2> problem(Wh, optionProblem);

   // CutFEM_R2 beta({R2(1,1), R2(1,1)});
   double lambdaB = sqrt(10);
   // CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
   // double lambdaB = 3.;//sqrt(10);
   CutFEMParameter lambda(0, 1.);

   Normal N;
   FunTest v(Wh, 1);
   Fun_h fun_Un(Wh, u0);
   Expression Un(fun_Un, 0, op_id);

   Fun_h gh(Wh, fun_solution, tn);
   Expression gx(gh, 0, op_id);

   CutFEM_R2 beta({R2(3, 1), R2(3, 1)});

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
   mAh.addMatMul(u0, problem.rhs_);

   // matlab::Export(problem.rhs_, "rhs0.dat");
   // problem.rhs_ = 0.;
   // CONSTRUCT THE RHS
   double tt = getTime();

   problem.addLinear(innerProduct(fluxFctX(Un), dx(v)) +
                         innerProduct(fluxFctY(Un), dy(v)),
                     Khi);

   problem.addLinear(
       // - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
       -innerProduct(average(fluxFctN(Un)), jump(v)), Khi,
       INTEGRAL_INNER_EDGE_2D);

   problem.addLinear(-jump(innerProduct(fluxFctN(Un), v)) -
                         innerProduct(jump(fluxFctN(Un)), jump(lambda * v)),
                     interface);
   // matlab::Export(problem.rhs_, "rhs1.dat");

   // BOUNDARY CONDITION IN RHS
   problem.addLinear(-innerProduct(fluxFctN(Un), 0.5 * v) -
                         innerProduct(Un, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY
                     // , {1,4}          // boundary in
   );
   problem.addLinear(-innerProduct(fluxFctN(gx), 0.5 * v) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY
                     // , {1,4}          // label left boundary
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

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   tt = getTime();
   problem.solve(Mh, problem.rhs_);
   // std::cout << "Time solver \t" << getTime()-tt << std::endl;

   uh           = problem.rhs_;
   problem.rhs_ = 0.;
}

int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_burgerCutFEM_nx40_P1.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 40;
   int ny = 40;
   Mesh2 Th(nx, ny, 0., 0., 2., 2.);
   Space Vh(Th, DataFE<Mesh2>::P1dc);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / (nx + 1);
   double dt = meshSize / 3 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend    = 0.1;                      // 1/pi;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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
   std::cout << "Number of dof = \t" << Wh.NbDoF() << std::endl;

   MacroElement<Mesh> macro(Khi, 0.5);

   // DECLARATION OF THE VECTOR CONTAINING THE solution
   // =====================================================
   Rn u0(Wh.NbDoF(), 0.);
   interpolate(Wh, u0, fun_initial);
   Rn uh(u0);
   Rn uh_tild(u0);

   Fun_h fun_u0(Wh, u0);
   Fun_h fun_uh(Wh, uh);
   Fun_h fun_uh_tild(Wh, uh_tild);

   // Plot the macro elements
   //   {
   //     Paraview<Mesh> writer(Khi, "test_burger_mesh.vtk");
   //     writer.add(fun_uh, "uh", 0, 1);
   //     writer.add(levelSet, "levelSet", 0, 1);
   //     writer.writeMacroInnerEdge (macro, 0,
   //     "test_burger_macro_inner_edge1.vtk");
   //     writer.writeMacroOutterEdge(macro, 0,
   //     "test_burger_macro_outter_edge1.vtk"); writer.writeMacroInnerEdge
   //     (macro, 1, "test_burger_macro_inner_edge2.vtk");
   //     writer.writeMacroOutterEdge(macro, 1,
   //     "test_burger_macro_outter_edge2.vtk");
   //   }
   // return 0;
   double qu0    = 0.; // integral(Khi, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Wh, interface, Ah, Mh);

   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      double tt_iter = getTime();

      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
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

         if (Wh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild,
                                                         min_u0, max_u0, macro);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild,
                                                         min_u0, max_u0, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild,
                                                         min_u0, max_u0, macro);
         } else if (Wh.basisFctType == BasisFctType::P0) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
         }
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      double tt             = getTime();
      double qu             = integral(Khi, fun_uh_tild, 0);
      auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
      auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);
      std::cout << " time min max \t" << getTime() - tt << std::endl;

      if (min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0 + Epsilon) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

         // std::cout << min_u0 - min_u1 << std::endl;
         // std::cout << max_u0 - max_u1 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi,
         //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
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
      Fun_h fun_ex(Wh, usol);

      Expression2 femErr(femErrh, 0, op_id);
      R errU = sqrt(integral(Khi, femErr * femErr));
      errSum += errU;

      // PLOT THE SOLUTION
      // ==================================================
      if (MPIcf::IamMaster() && i % 5 == 0 || i + 1 == niteration) {
         // Fun_h fun_thet (Wh, fun_theta);
         Paraview<Mesh> writer(Khi, "test_burger_P1_" + std::to_string(ifig++) +
                                        ".vtk");
         writer.add(fun_uh, "uhNoLimiter", 0, 1);
         writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         writer.add(femErrh, "error", 0, 1);
         writer.add(fun_ex, "u_exact", 0, 1);
         // writer.add(fun_thet , "yx", 0, 1);
      }

      std::cout << " || u-uex ||_2 = " << errU << std::endl;
      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU
                 << "\t" << std::setprecision(16) << qu << "\t"
                 << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_u1 << "\t"
                 << std::setprecision(16) << max_u1 << "\t"
                 << std::setprecision(5) << std::endl;

      std::cout << " Iteration computation time \t" << getTime() - tt_iter
                << std::endl;
   }

   std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration
             << std::endl;
   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef FEM_DISCONTINUOUS_SOLUTION
R fun_initial(double *P, int elementComp, int domain) {
   double xs = 0.9, ys = 0.9;
   double r = 0.3;
   double v = (P[0] - xs) * (P[0] - xs) + (P[1] - ys) * (P[1] - ys);
   if (v < r * r)
      return 1;
   else
      return 0.;
}
R fun_boundary(double *P, int elementComp, int domain, double t) { return 0.; }

void assembly(const Space &Wh, MatMap &Ah, MatMap &Mh) {

   double t0 = getTime();
   const Mesh &Khi(Wh.Th);
   FEM<Mesh> problem(Wh);
   R2 beta(1, 1);
   double lambdaE = sqrt(10);
   const MeshParameter &h(Parameter::h);

   double lambdaB = sqrt(10);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);

   problem.addBilinear(innerProduct(beta * u, grad(v)), Khi);

   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(average(beta * u * n), jump(v)) -
                           innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                       Khi, INTEGRAL_INNER_EDGE_2D);

   problem.addBilinear(-innerProduct(beta * u * n, v), Khi, INTEGRAL_BOUNDARY,
                       {2, 3} // label other boundary
   );
   problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) -
                           innerProduct(u, lambdaB * 0.5 * v),
                       Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // BUILDING (M + S)
   // =====================================================
   problem.set_map(Mh);
   problem.addBilinear(innerProduct(u, v), Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}

void solve_problem(const Space &Wh, const Rn &u0, Rn &uh, MatMap &Ah,
                   MatMap &Mh, double tn) {

   double t0 = getTime();

   const Mesh &Khi(Wh.Th);

   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   FEM<Mesh> problem(Wh, optionProblem);

   R2 beta(1, 1);
   double lambdaB = sqrt(10);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   Fun_h gh(Wh, fun_boundary, tn);
   Expression2 gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   int N = problem.get_nb_dof();
   multiply(N, N, Ah, u0, problem.rhs_);

   problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   problem.solve(Mh, problem.rhs_);

   uh = problem.rhs_;
}

int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_FEM_nx20_DiscontinuousSol_P1dc.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 20;
   int ny = 20;
   // Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
   Mesh2 Th(nx, ny, 0., 0., 2., 2.);

   Space Vh(Th, DataFE<Mesh2>::P1dc);
   Space Eh(Th, DataFE<Mesh2>::P0);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / nx;
   double dt   = meshSize / 5 / sqrt(2) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend = 0.1;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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

   double qu0    = integral(Th, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Vh, Ah, Mh);

   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
                << " \t time = " << tid << std::endl;

      // THIRD ORDER RK
      // =================================================
      {
         Rn u1(u0);
         Rn u2(Vh.NbDoF(), 0.);
         Rn u2_tild(Vh.NbDoF(), 0.);
         Rn u1_tild(Vh.NbDoF(), 0.);
         Fun_h fun_u1(Vh, u1);
         Fun_h fun_u2(Vh, u2);

         if (Vh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Vh, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::FEM::applyBoundPreservingLimiter(fun_u1, u1_tild, min_u0,
                                                      max_u0);
            std::vector<double> indicator =
                limiter::FEM::compute_trouble_cell_indicator(fun_u1);

            uh                    = u1;
            uh_tild               = u1_tild;
            auto [min_uh, max_uh] = limiter::FEM::findMinAndMaxValue(fun_uh);
            auto [min_u1, max_u1] =
                limiter::FEM::findMinAndMaxValue(fun_uh_tild);
            std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                      << std::endl;
            std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
            std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

            Fun_h fun_indicator(Eh, indicator);

            if (MPIcf::IamMaster()) {
               Paraview<Mesh> writer(Th, "checkErrorFEMDiscontinuous.vtk");
               writer.add(fun_uh, "uhNoLimiter", 0, 1);
               writer.add(fun_uh_tild, "uhLimiter", 0, 1);
               writer.add(fun_indicator, "indicator", 0, 1);
            }
            return 0;

            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Vh, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::FEM::applyBoundPreservingLimiter(fun_u2, u2_tild, min_u0,
                                                      max_u0);
            solve_problem(Vh, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::FEM::applyBoundPreservingLimiter(fun_uh, uh_tild, min_u0,
                                                      max_u0);
         } else if (Vh.basisFctType == BasisFctType::P0) {
            solve_problem(Vh, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            u2 += 3. / 4 * u0 + 1. / 4 * u1;
            solve_problem(Vh, u1, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            solve_problem(Vh, u2, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2;
            uh_tild = uh;
         }

         // limiter::CutFEM::check_mean_value(fun_u1, macro, 0.,1., u_mean);
         // Fun_h fun_u_mean(Eh, u_mean);

         // limiter::CutFEM::extendToMacroP1(fun_u1, uM, map_mean_value, macro);
         // std::cout << "hey " << std::endl;
         // Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
         // writer.add(fun_u1, "u0", 0, 1);
         // writer.add(fun_uM, "macroextended", 0, 1);
         // writer.add(fun_u_mean, "uMean", 0, 1);

         // uM = u1;
         // limiter::CutFEM::limiter_Pei_P1(fun_u1, uM, min_u0, max_u0, macro);
         // limiter::CutFEM::minmaxP1(fun_uM, min_u0, max_u0);

         // std::cout << min_u0 << "\t" << max_u0 << std::endl;
         // // if(MPIcf::IamMaster()) {
         //   // Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
         //   writer.add(fun_uM, "macroLimited", 0, 1);
         // // }

         // getchar();
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      double tt             = getTime();
      double qu             = integral(Th, fun_uh_tild, 0);
      auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
      auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);
      std::cout << " time min max \t" << getTime() - tt << std::endl;

      // if((i==24) && MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_u1, "uhLimiter", 0, 1);
      //   writer.add(fun_uM, "macroExtend", 0, 1);
      // }

      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi,
      //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1); writer.add(fun_uh_tild,
      //   "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }

      // if( (min_u0-Epsilon) < min_u1 && max_u1<
      // (max_u0+Epsilon)) {
      if (min_u0 <= min_u1 && max_u1 <= max_u0) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi,
         //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
         //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
         //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         //   // writer.add(fun_uM, "macroExtend", 0, 1);
         // }
         // return 0;
      }
      // PLOT THE SOLUTION
      // ==================================================
      // if(MPIcf::IamMaster() && i%1 == 0 || i+1 == niteration) {
      //   Paraview<Mesh> writer(Khi,
      //   "conservationP1_"+std::to_string(ifig++)+".vtk"); writer.add(fun_uh ,
      //   "uhNoLimiter", 0, 1); writer.add(fun_uh_tild, "uhLimiter"  , 0, 1);
      // }

      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << qu
                 << "\t" << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_u1 << "\t"
                 << std::setprecision(16) << max_u1 << "\t"
                 << std::setprecision(5) << std::endl;
      // getchar();
      // return 0;
      // if(i==10) return 0;
   }

   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif

#ifdef CUTFEM_DISCONTINUOUS_SOLUTION
double c0 = 0.77; // 0.25;
R fun_levelSet(double *P, const int i) { return -(P[1] + 0.5 * P[0] - 1.73); }
R fun_initial(double *P, int elementComp, int domain) {
   double xs = 0.9, ys = 0.9;
   double r = 0.3;
   double v = (P[0] - xs) * (P[0] - xs) + (P[1] - ys) * (P[1] - ys);
   if (v < r * r)
      return 1; // exp(-pow(P[0] - xs,2)/r - pow(P[1] - ys,2)/r);
   else
      return 0.;
}
R fun_boundary(double *P, int elementComp, int domain, double t) {
   return 0.;
   // return sin(pi*(P[0]+P[1] - 4*t));
}

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah,
              MatMap &Mh) {

   double t0 = getTime();
   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
   CutFEM<Mesh2> problem(Wh);
   CutFEM_R2 beta({R2(1, 1), R2(1, 1)});
   CutFEMParameter lambda(0, 1.);
   CutFEMParameter lambdaE(sqrt(10), sqrt(10)); // max B = sqrt(10)/sqrt(5)
   const MeshParameter &h(Parameter::h);

   double lambdaB = sqrt(10);
   double Cstab   = 1e-2;
   double Cstabt  = 1e-2;

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);

   // BUILDING A
   // =====================================================
   problem.set_map(Ah);
   problem.addBilinear(innerProduct(beta * u, grad(v)), Khi);

   problem.addBilinear(-jump(innerProduct(beta * u, v * n)) -
                           innerProduct(jump(beta * u * n), jump(lambda * v)),
                       interface);
   // F(u)_e = {B.u} - lambda_e/2 [u]
   problem.addBilinear(-innerProduct(average(beta * u * n), jump(v)) -
                           innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                       Khi, INTEGRAL_INNER_EDGE_2D);

   //[ u \beta \cdot n]=3u_1-2u_2=0
   problem.addFaceStabilization(
       -innerProduct(jump(u), Cstab * jump(v)) -
           innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
       // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   problem.addBilinear(-innerProduct(beta * u * n, v), Khi, INTEGRAL_BOUNDARY,
                       {2, 3} // label other boundary
   );
   problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) -
                           innerProduct(u, lambdaB * 0.5 * v),
                       Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // BUILDING (M + S)
   // =====================================================
   problem.set_map(Mh);
   problem.addBilinear(innerProduct(u, v), Khi);
   problem.addFaceStabilization(
       innerProduct(h * jump(u), Cstab * jump(v)) +
           innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v)))
       // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
       ,
       Khi);

   std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Interface<Mesh> &interface,
                   const Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh, double tn) {

   double t0 = getTime();

   const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

   ProblemOption optionProblem;
   optionProblem.solver_name_  = "umfpack";
   optionProblem.clear_matrix_ = false;
   CutFEM<Mesh2> problem(Wh, optionProblem);

   CutFEM_R2 beta({R2(1, 1), R2(1, 1)});
   double lambdaB = sqrt(10);

   Normal n;
   FunTest u(Wh, 1), v(Wh, 1);
   Fun_h gh(Wh, fun_boundary, tn);
   Expression2 gx(gh, 0, op_id);

   // MULTIPLYING A * u_0  => rhs
   // =====================================================
   int N = problem.get_nb_dof();
   multiply(N, N, Ah, u0, problem.rhs_);

   problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) +
                         innerProduct(gx, lambdaB * 0.5 * v),
                     Khi, INTEGRAL_BOUNDARY, {1, 4} // label left boundary
   );

   // SOLVING  (M+S)E'(t) = rhs
   // =====================================================
   problem.solve(Mh, problem.rhs_);

   uh = problem.rhs_;
}

int main(int argc, char **argv) {

   MPIcf cfMPI(argc, argv);
   const double cpubegin = getTime();

   // OUTPUT FILE
   // =====================================================
   std::ofstream outputData("output_nx20_DiscontinuousSol_P1dc.txt");

   // DEFINITION OF THE MESH and SPACE
   // =====================================================
   int nx = 20;
   int ny = 20;
   // Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
   Mesh2 Th(nx, ny, 0., 0., 2., 2.);

   Space Vh(Th, DataFE<Mesh2>::P1dc);
   Space Eh(Th, DataFE<Mesh2>::P0);

   // DEFINITION OF SPACE AND TIME PARAMETERS
   // =====================================================
   double tid      = 0;
   double meshSize = 2. / nx;
   double dt   = meshSize / 5 / sqrt(2) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
   double tend = 0.1;
   int niteration = tend / dt;
   dt             = tend / niteration;
   double errSum  = 0;

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

   MacroElement<Mesh> macro(Khi, 0.25);
   // getchar();
   // DECLARATION OF THE VECTOR CONTAINING THE solution
   // =====================================================
   Rn u0(Wh.NbDoF(), 0.);
   interpolate(Wh, u0, fun_initial);
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
      writer.writeMacroInnerEdge(macro, 0, "pei_macro_inner_edge1.vtk");
      writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
      writer.writeMacroInnerEdge(macro, 1, "pei_macro_inner_edge2.vtk");
      writer.writeMacroOutterEdge(macro, 1, "pei_macro_outter_edge2.vtk");
   }
   // return 0;

   double qu0    = integral(Khi, fun_uh, 0);
   double min_u0 = uh.min();
   double max_u0 = uh.max();

   // ASSEMBLY THE CONSTANT PART
   // ==================================================
   MatMap Ah, Mh;
   assembly(Wh, interface, Ah, Mh);

   // RESOLUTION OF THE PROBLEM_MIXED_DARCY
   // ==================================================
   int ifig = 1;
   for (int i = 0; i < niteration; ++i) {
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "Iteration " << i + 1 << " / " << niteration
                << " \t time = " << tid << std::endl;

      // EULER METHOD
      // =================================================
      // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
      // uh *= dt;
      // uh += u0;
      // // std::cout << uh << std::endl;
      // limiter::CutFEM::limiter_Pei_P1(fun_uh, uh_tild, min_u0, max_u0,
      // macro);

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

         // Rn u_mean(Wh.get_nb_element());
         // std::map<int, double>  map_mean_value;
         // Rn uM(u0);
         // Fun_h fun_uM(Wh, uM);
         // double min_u0=0., max_u0=1.;

         // std::cout << u0 << std::endl;
         // limiter::CutFEM::check_mean_value(fun_u0, macro, 0., 1., u_mean);
         // limiter::CutFEM::limiter_Pei_P1(fun_u0, uM, min_u0, max_u0, macro);
         // limiter::CutFEM::extendToMacroP1(fun_u0, uM, map_mean_value, macro);
         // limiter::CutFEM::minmaxP1(fun_uM, min_u0, max_u0);
         // std::cout << min_u0 << "\t" << max_u0 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
         //   writer.add(fun_u0, "u0", 0, 1);
         //   writer.add(fun_uM, "macroExtend", 0, 1);
         // }
         // limiter::CutFEM::limiter_Pei_P1(fun_u0, uM, min_u0, max_u0, macro);
         // limiter::CutFEM::check_mean_value(fun_u0, macro, 0., 1., u_mean);

         // getchar();
         // u0 = uM;

         // std::cout << " compute u1 " << std::endl;
         // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
         // u1 += dt * uh;
         // uM = u1;

         if (Wh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, u1_tild,
                                                         min_u0, max_u0, macro);

            std::vector<double> indicator =
                limiter::CutFEM::compute_trouble_cell_indicator(fun_u1);

            uh                    = u1;
            uh_tild               = u1_tild;
            auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
            auto [min_u1, max_u1] =
                limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);
            std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                      << std::endl;
            std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
            std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

            Fun_h fun_indicator(Eh, u2);

            if (MPIcf::IamMaster()) {
               Paraview<Mesh> writer(Khi, "checkErrorDiscontinuous.vtk");
               writer.add(fun_uh, "uhNoLimiter", 0, 1);
               writer.add(fun_uh_tild, "uhLimiter", 0, 1);
            }
            return 0;

            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_u2, u2_tild,
                                                         min_u0, max_u0, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, uh_tild,
                                                         min_u0, max_u0, macro);
         } else if (Wh.basisFctType == BasisFctType::P0) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            u1 += dt * uh;
            limiter::CutFEM::extendToMacro(fun_u1, u1_tild, macro);
            u2 += 3. / 4 * u0 + 1. / 4 * u1_tild;
            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            u2 += 1. / 4 * dt * uh;
            limiter::CutFEM::extendToMacro(fun_u2, u2_tild, macro);
            solve_problem(Wh, interface, u2_tild, uh, Ah, Mh, tid + 0.5 * dt);
            uh *= 2. / 3 * dt;
            uh += 1. / 3 * u0 + 2. / 3 * u2_tild;
            limiter::CutFEM::extendToMacro(fun_uh, uh_tild, macro);
         }

         // limiter::CutFEM::check_mean_value(fun_u1, macro, 0.,1., u_mean);
         // Fun_h fun_u_mean(Eh, u_mean);

         // limiter::CutFEM::extendToMacroP1(fun_u1, uM, map_mean_value, macro);
         // std::cout << "hey " << std::endl;
         // Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
         // writer.add(fun_u1, "u0", 0, 1);
         // writer.add(fun_uM, "macroextended", 0, 1);
         // writer.add(fun_u_mean, "uMean", 0, 1);

         // uM = u1;
         // limiter::CutFEM::limiter_Pei_P1(fun_u1, uM, min_u0, max_u0, macro);
         // limiter::CutFEM::minmaxP1(fun_uM, min_u0, max_u0);

         // std::cout << min_u0 << "\t" << max_u0 << std::endl;
         // // if(MPIcf::IamMaster()) {
         //   // Paraview<Mesh> writer(Khi, "maxPrincipleU0.vtk");
         //   writer.add(fun_uM, "macroLimited", 0, 1);
         // // }

         // getchar();
      }

      u0 = uh_tild;
      tid += dt;

      // COMPUTATION OF THE L2 ERROR
      // =================================================
      double tt             = getTime();
      double qu             = integral(Khi, fun_uh_tild, 0);
      auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
      auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);
      std::cout << " time min max \t" << getTime() - tt << std::endl;

      // if((i==24) && MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi, "maxPrinciple.vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
      //   writer.add(fun_u1, "uhLimiter", 0, 1);
      //   writer.add(fun_uM, "macroExtend", 0, 1);
      // }

      // if(MPIcf::IamMaster()) {
      //   Paraview<Mesh> writer(Khi,
      //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
      //   writer.add(fun_uh, "uhNoLimiter", 0, 1); writer.add(fun_uh_tild,
      //   "uhLimiter", 0, 1);
      //   // writer.add(fun_uM, "macroExtend", 0, 1);
      // }

      // if( (min_u0-Epsilon) < min_u1 && max_u1<
      // (max_u0+Epsilon)) {
      if (min_u0 <= min_u1 && max_u1 <= max_u0) {

         std::cout << " Maximum principle satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
      } else {
         std::cout << " Maximum principle not satified! " << std::endl;
         std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]"
                   << std::endl;
         std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
         std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;

         // if(MPIcf::IamMaster()) {
         //   Paraview<Mesh> writer(Khi,
         //   "maxPrincipleCutFEM_"+std::to_string(ifig++)+".vtk");
         //   writer.add(fun_uh, "uhNoLimiter", 0, 1);
         //   writer.add(fun_uh_tild, "uhLimiter", 0, 1);
         //   // writer.add(fun_uM, "macroExtend", 0, 1);
         // }
         // return 0;
      }
      // PLOT THE SOLUTION
      // ==================================================
      if (MPIcf::IamMaster() && i % 1 == 0 || i + 1 == niteration) {
         Paraview<Mesh> writer(Khi, "conservationP1_" + std::to_string(ifig++) +
                                        ".vtk");
         writer.add(fun_uh, "uhNoLimiter", 0, 1);
         writer.add(fun_uh_tild, "uhLimiter", 0, 1);
      }

      std::cout << std::setprecision(16)
                << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                << std::endl;

      if (i == 0) {
         outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU "
                    << std::endl;
      }
      outputData << i << "\t" << tid << "\t" << std::setprecision(16) << qu
                 << "\t" << std::setprecision(16) << fabs(qu - qu0) << "\t"
                 << std::setprecision(16) << min_u1 << "\t"
                 << std::setprecision(16) << max_u1 << "\t"
                 << std::setprecision(5) << std::endl;
      // getchar();
      // return 0;
      // if(i==10) return 0;
   }

   std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
   return 0;
}

#endif
