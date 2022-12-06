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
#include "../util/cputime.h"
#ifdef USE_MPI
#include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "../num/matlab.hpp"

// #include "../num/gnuplot.hpp"

// #define PROBLEM_STOKES_DIV
#define POISSEUILLE_EXAMPLE
// #define STATIC_DROP_EXAMPLE
// #define FICTITIOUS_STOKES
// #define DYNAMIC_DROP_EXAMPLE

#ifdef PROBLEM_STOKES_DIV

/* REMAINING ISSUES & CURRENT KNOWLEDGE:
= Normal must point outward from domain
= Normal must point outward from + element in jump term
= Adding terms on the border should not be done when sol is zero there
= The BC is ignored when sigma >> penparam
= Big sigma seems to steer the velocity SE to incorrect sols - eg sigma*100
gives vel/100
== Get correct sol (u,p) for sigma=1, pen>1e3/meshsize, with BC normal pen & BC
tangent natural
== Same sol if include all natural BC (ie avg and jumps)
== Same sol if remove BC tangent natural, but then the pressure gets the NE SW
corner issue " For example, in the H(div) method, the normal component of the
velocity is set as an essential boundary condition and is strongly enforced, but
the tangential component of the velocity is treated as a natural boundary
condition and is weakly imposed. "


- Velocity is incorrect (weird symmetry along y=-x+c)

+ Pressure gets corrected by adding interior edge terms also on the boundary
*/
namespace Erik_Data_StokesDiv {

// R fun_rhs(const R2 P, int i) {
//   return 0;
// }
// R fun_exact(const R2 P, int i) {
//   if(i==0) return 0;
//   else if(i==1) return P.x;
//   else return 0;
// }

// R fun_rhs(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//    // if(i==0) return -2*y;
//    // else return 2*x;
//   if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) -
//   12*(x-1)*(x-1)*x*x*(2*y - 1); else if(i==1) return 4*x*(2*x*x - 3 *x +
//   1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y; else return 0;
// }
// R fun_exact(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return 2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   else return 0;
// }
// R fun_exact_u(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   // if(i==0)      return  x*x*y;
//   // else if(i==1) return -x*y*y;
// }
// R fun_exact_p(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }
// R fun_rhs(const R2 P, int i) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1)
//   return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
//   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }

R fun_div(const R2 P, int i) {
   R x = P.x;
   R y = P.y;
   return 0;
}

R fun_rhs(const R2 P, int i) {
   // R mu=1;
   R x = P.x;
   R y = P.y;
   if (i == 0)
      return 40 * x * (x * x - y * y) - (16 * y - 8);
   else if (i == 1)
      return (16 * x - 8) - 40 * y * (x * x - y * y);
   else
      return 0;
}
R fun_exact_u(const R2 P, int i) {
   R x = P.x;
   R y = P.y;
   if (i == 0)
      return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
   else
      return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}
R fun_exact_p(const R2 P, int i) {
   R x = P.x;
   R y = P.y;
   return (10 * (x * x - y * y) *
           (x * x - y * y)); // - 670.171/(pi*interfaceRad*interfaceRad);
}
} // namespace Erik_Data_StokesDiv
using namespace Erik_Data_StokesDiv;

int main(int argc, char **argv) {
   typedef TestFunction<2> FunTest;
   typedef FunFEM<Mesh2> Fun_h;
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;

   const double cpubegin = CPUtime();
   MPIcf cfMPI(argc, argv);

   int nx = 10;
   int ny = 10;
   // int d = 2;
   ProblemOption optionStokes;
   optionStokes.order_space_element_quadrature_ = 9;
   vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

   for (int i = 0; i < 3; ++i) { // i<3

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Th(nx, ny, 0., 0., 1., 1.);

      Lagrange2 FEvelocity(4);
      Space Uh(Th, FEvelocity);
      // Space Qh(Th, DataFE<Mesh>::P1);

      Space Vh(Th, FEvelocity);
      // Space Vh(Th, DataFE<Mesh2>::RT1); // REMOVE THESE TWO LATER
      Space Qh(Th, DataFE<Mesh2>::P3dc); // FOR MIXEDSPACE

      const R hi           = Th[0].lenEdge(0);
      const R penaltyParam = 1e3 / pow(hi, 1);
      const R sigma        = 0.5; // HIGHER 1e2,1e3,etc WORSENS THE SOL [??]

      Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h gh(Uh, fun_exact_u);
      // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

      FEM<Mesh2> stokes(Vh, optionStokes);
      stokes.add(Qh);

      Normal n;
      Tangent t;
      double kappa1 = 0.5, kappa2 = 0.5;
      /* Syntax:
      FunTest (fem space, #components, place in space)
      */
      FunTest u(Vh, 2, 0), p(Qh, 1, 0), v(Vh, 2, 0), q(Qh, 1, 0);

      R mu = 1;
      // const MeshParameter& h(Parameter::h);
      // const MeshParameter& invlEdge(Parameter::invmeas);

      // a_h(u,v)
      stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) -
                             innerProduct(p, div(v)) + innerProduct(div(u), q),
                         Th);

      stokes.addBilinear(
          -innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2),
                        jump(v * t)) +
              innerProduct(jump(u * t),
                           average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
              innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
          Th, INTEGRAL_INNER_EDGE_2D);

      stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                             + innerProduct(u, 2 * mu * Eps(v) * n) +
                             innerProduct(p, v * n) // natural
                             // - innerProduct(u*n, q)  // essential
                             + innerProduct(penaltyParam * u, v),
                         Th, INTEGRAL_BOUNDARY);

      // l(v)_Omega
      stokes.addLinear(+innerProduct(gh.expression(2), penaltyParam * v) +
                           innerProduct(gh.expression(2), 2 * mu * Eps(v) * n)
                       // - innerProduct(gh.expression(2), q*n)
                       ,
                       Th, INTEGRAL_BOUNDARY);

      stokes.addLinear(innerProduct(fh.expression(2), v), Th);

      stokes.addLagrangeMultiplier(innerProduct(1., p), 0., Th);

      stokes.solve();

      // EXTRACT SOLUTION
      int idx0_s  = Vh.get_nb_dof();
      Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
      Rn_ data_ph = stokes.rhs_(SubArray(Qh.get_nb_dof(), idx0_s));
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Qh, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      {
         Fun_h soluErr(Vh, fun_exact_u);
         Fun_h soluh(Vh, fun_exact_u);
         soluErr.v -= uh.v;
         soluErr.v.map(fabs);
         // Fun_h divSolh(Wh, fun_div);
         // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

         Paraview<Mesh> writer(Th, "stokes_" + to_string(i) + ".vtk");
         writer.add(uh, "velocity", 0, 2);
         writer.add(ph, "pressure", 0, 1);
         writer.add(femSol_0dx + femSol_1dy, "divergence");
         writer.add(soluh, "velocityExact", 0, 2);
         writer.add(soluErr, "velocityError", 0, 2);
         // writer.add(solh, "velocityError" , 0, 2);

         // writer.add(fabs(femDiv, "divergenceError");
      }

      // Rn solExVec(mixedSpace.NbDoF());
      // interpolate(mixedSpace, solExVec, fun_exact);
      // R errU = L2norm(sol,fun_exact,0,2);
      // R errP = L2norm(sol,fun_exact,2,1);

      R errU      = L2norm(uh, fun_exact_u, 0, 2);
      R errP      = L2norm(ph, fun_exact_p, 0, 1);
      R errDiv    = L2norm(femSol_0dx + femSol_1dy, Th);
      R maxErrDiv = maxNorm(femSol_0dx + femSol_1dy, Th);

      // // solExVec -= stokesDiv.rhs;
      // for(int i=0;i<solExVec.size();++i){
      //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
      // }
      //
      // Fun_h solEx(mixedSpace, solExVec);
      // writerS.add(solEx,"uh_err", 0, 2);

      // writerS.add(solEx,"uh_ex", 0, 2);
      // writerS.add(solEx,"ph_ex", 2, 1);

      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2.push_back(errDiv);
      divmax.push_back(maxErrDiv);
      h.push_back(1. / nx);
      if (i == 0) {
         convu.push_back(0);
         convp.push_back(0);
      } else {
         convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
         convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
      }

      nx *= 2;
      ny *= 2;
   }
   std::cout << "\n"
             << std::left << std::setw(10) << std::setfill(' ') << "h"
             << std::setw(15) << std::setfill(' ') << "err_p u" << std::setw(15)
             << std::setfill(' ') << "conv p" << std::setw(15)
             << std::setfill(' ') << "err u" << std::setw(15)
             << std::setfill(' ') << "conv u" << std::setw(15)
             << std::setfill(' ')
             << "err divu"
             // << std::setw(15) << std::setfill(' ') << "conv divu"
             // << std::setw(15) << std::setfill(' ') << "err_new divu"
             // << std::setw(15) << std::setfill(' ') << "convLoc divu"
             << std::setw(15) << std::setfill(' ')
             << "err maxdivu"
             // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
             << "\n"
             << std::endl;
   for (int i = 0; i < h.size(); ++i) {
      std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i]
                << std::setw(15) << std::setfill(' ') << pl2[i] << std::setw(15)
                << std::setfill(' ') << convp[i] << std::setw(15)
                << std::setfill(' ') << ul2[i] << std::setw(15)
                << std::setfill(' ') << convu[i] << std::setw(15)
                << std::setfill(' ')
                << divl2[i]
                // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                << std::setw(15) << std::setfill(' ')
                << divmax[i]
                // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                << std::endl;
   }
}
#endif

#ifdef DYNAMIC_DROP_EXAMPLE
namespace DynamicDropData {
R2 shift(0, 0);
double sigma = 1. / 3;
double mu1   = 10;
double mu2   = 1;

R fun_levelSet(const R2 P, int i) {
   return sqrt((P.x - shift.x) * (P.x - shift.x) +
               (P.y - shift.y) * (P.y - shift.y)) -
          2. / 3;
}

static R falpha1(const R r) {
   // const R MU1 = 1e-1;
   // const R MU2 = 1e-3;
   const R MU1 = mu1;
   const R MU2 = mu2;
   const R r0  = 2. / 3;
   return 1. / MU1 + (1. / MU2 - 1. / MU1) * exp(r * r - r0 * r0);
}
static R falpha2(const R r) {
   R MU2 = mu2; // 1e-3;
   return 1. / MU2;
}

static R falpha(const R r) {
   const R r0 = 2. / 3;
   return (r < r0) ? falpha2(r) : falpha1(r);
}

static R fun_rhs(const R2 P, int i) {
   const R r2 = Norme2_2(P - shift);
   const R s  = exp(-r2);
   R2 R(4 * s * (r2 - 2) * P.y + 3 * P.x * P.x, -4 * s * (r2 - 2) * P.x);
   return (i < 2) ? R[i] : 0;
}
static R fun_boundary(const R2 P, const int i) {
   const R r = Norme2(P - shift);
   R2 R(-P.y, P.x);
   R = falpha(r) * exp(-r * r) * R;
   return (i < 2) ? R[i] : 0;
}

static R2 fun_velocity1(const R2 P) {
   R r = Norme2(P - shift);
   R2 R(-P.y, P.x);
   R = falpha1(r) * exp(-r * r) * R;
   return R;
}
static R2 fun_velocity2(const R2 P) {
   R r = Norme2(P - shift);
   R2 R(-P.y, P.x);
   R = falpha2(r) * exp(-r * r) * R;
   return R;
}
static R fun_pressure1(const R2 P) { return pow(P.x, 3); }
static R fun_pressure2(const R2 P) {
   // R sigma = 1;//24.5;//700;
   return pow(P.x, 3) + sigma * 3. / 2.;
}

// static R fun_solution(const R2 P, int ci, int domain) {
//   if(domain == 0) return (ci < 2)? fun_velocity1(P)[ci] : fun_pressure1(P);
//   else return (ci < 2)? fun_velocity2(P)[ci] : fun_pressure2(P);
// }
static R fun_exact_u(const R2 P, int ci, int domain) {
   if (domain == 0)
      return fun_velocity1(P)[ci];
   else
      return fun_velocity2(P)[ci];
}
static R fun_exact_p(const R2 P, int ci, int domain) {
   if (domain == 0)
      return fun_pressure1(P);
   else
      return fun_pressure2(P);
}

R2 fparam(double t) {
   return R2(2. / 3 * cos(t + 1. / 3), 2. / 3 * sin(t + 1. / 3));
}
R fun_1(const R2 P, const int i) { return 1; }
} // namespace DynamicDropData
using namespace DynamicDropData;
namespace Erik_Data_StokesDiv {
double mu1 = 1;
double mu2 = 1;
R2 shift(0., 0.);
double sigma = 1;
// R interfaceRad = 0.250001;//2./3; // not exactly 1/4 to avoid interface
// cutting exaclty a vertex R fun_levelSet(const R2 P, const int i) {
//   return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) -
//   interfaceRad;
// }
R fun_levelSet(const R2 P, int i) {
   return sqrt((P.x - shift.x) * (P.x - shift.x) +
               (P.y - shift.y) * (P.y - shift.y)) -
          2. / 3;
}

// [Ex: divu = 0, inhomog. BC]
R fun_rhs(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   if (i == 0)
      return -2 * y;
   else
      return 2 * x;
   // if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) -
   // 12*(x-1)*(x-1)*x*x*(2*y - 1); else if(i==1) return 4*x*(2*x*x - 3 *x +
   // 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y; else return 0;
}
R fun_exact_u(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
   // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
   if (i == 0)
      return x * x * y;
   else if (i == 1)
      return -x * y * y;
}
R fun_exact_p(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   return 0.;
}
R fun_div(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   return 0;
}

// R fun_rhs(const R2 P, int i) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1)
//   return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
//   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
} // namespace Erik_Data_StokesDiv
// using namespace Erik_Data_StokesDiv;

int main(int argc, char **argv) {
   typedef TestFunction<2> FunTest;
   typedef FunFEM<Mesh2> Fun_h;
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;

   const double cpubegin = CPUtime();
   MPIcf cfMPI(argc, argv);

   int nx = 10;
   int ny = 10;
   // int d = 2;

   ProblemOption optionStokes;
   optionStokes.order_space_element_quadrature_ = 5;

   vector<double> ul2, pl2, divmax0, divmax1, divl2_0, divl2_1, h, convu, convp;
   int niter = 1;
   for (int i = 0; i < niter; ++i) { // i<3
      double t_iter_begin = MPIcf::Wtime();
      std::cout << "\n ------------------------------------- " << std::endl;
      std::cout << " Iteration \t" << i + 1 << "  / " << niter << std::endl;
      Mesh Kh(nx, ny, -1., -1., 2., 2.);
      const R hi           = 2. / (nx - 1);
      const R penaltyParam = 1e1 / hi;
      // const R sigma = 1;
      double uPenParam     = 1e-2;
      double pPenParam     = 1e0;
      double jumpParam     = 1e0;

      // CutFEMParameter mu(1e-1,1e-3);
      // CutFEMParameter invmu(1e1,1e3);
      CutFEMParameter mu(mu1, mu2);
      // double sigma = 0;//24.5;//700;
      double kappa1 = 0.5; // mu(1)/(mu(0)+mu(1));
      double kappa2 = 0.5; // mu(0)/(mu(0)+mu(1));

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);

      Lagrange2 FEvelocity(2);
      Space UUh(Kh, FEvelocity);

      Space Wh(Kh, DataFE<Mesh2>::RT1);
      Space Qh(Kh, DataFE<Mesh2>::P1dc);

      ActiveMesh<Mesh> Khi(Kh, interface);

      // ActiveMesh<Mesh> Kh0(Kh);
      // Kh0.createSurfaceMesh(interface);
      // ActiveMesh<Mesh> Kh1(Kh);
      // Kh1.truncate(interface, -1);
      // ActiveMesh<Mesh> Kh2(Kh);
      // Kh2.truncate(interface, 1);

      // Khi.truncate(interface, 1);
      CutSpace Vh(Khi, Wh);
      CutSpace Ph(Khi, Qh);
      CutSpace Uh(Khi, UUh);

      std::cout << " nb dof u \t" << Vh.get_nb_dof() << std::endl;
      std::cout << " nb dof p \t" << Ph.get_nb_dof() << std::endl;

      Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h gh(Uh, fun_exact_u);
      Fun_h ghRT(Vh, fun_exact_u);
      // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

      CutFEM<Mesh2> stokes(Vh, optionStokes);
      stokes.add(Ph);

      Normal n;
      Tangent t;
      FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0),
          p1(Ph, 1, 0, 0);
      double t_assembly_begin = MPIcf::Wtime();

      stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) -
                             innerProduct(p, div(v)) + innerProduct(div(u), q),
                         Khi);

      // stokes.addBilinear(
      //   innerProduct(div(u),q)
      //   , Khi
      //   , INTEGRAL_EXTENSION
      //   , 1.
      // );
      // matlab::Export(stokes.mat_, "matA.dat");

      // stokes.addBilinear(
      //   innerProduct(div(u),q)
      //   , Khi
      //   // , INTEGRAL_EXTENSION
      //   // , 1.
      // );
      //     matlab::Export(stokes.mat_, "matB.dat");

      // l(v)_Omega
      stokes.addLinear(innerProduct(fh.expression(2), v), Khi
                       // , INTEGRAL_EXTENSION
                       // , 1.
      );

      stokes.addBilinear(
          -innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2),
                        jump(v * t)) +
              innerProduct(jump(u * t),
                           average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
              innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
          Khi, INTEGRAL_INNER_EDGE_2D);
      stokes.addBilinear(
          innerProduct(jump(u), -2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
              innerProduct(-2 * mu * average(Eps(u) * n, kappa1, kappa2),
                           jump(v)) +
              innerProduct(1. / hi / hi * jump(u), jump(v)) +
              innerProduct(average(p, kappa1, kappa2), jump(v * n))
          // - innerProduct(jump(u*n), average(q,0.5,0.5))
          ,
          interface);
      // stokes.addBilinear(
      //   innerProduct(jump(u), -2*mu*average(grad(v)*n, 0.5,0.5))
      //   + innerProduct(-2*mu*average(grad(u)*n,0.5,0.5), jump(v))
      //   + innerProduct(10*jump(u), jump(v))
      //   + innerProduct(average(p,0.5,0.5), jump(v*n))
      //   - innerProduct(jump(u*n), average(q,0.5,0.5))
      //   , interface
      // );
      stokes.addLinear(innerProduct(3. / 2, average(v * n, kappa2, kappa1)) *
                           sigma,
                       interface);

      stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                             + innerProduct(p, v * n)          // natural

                             + innerProduct(u, 2 * mu * Eps(v) * n) +
                             innerProduct(penaltyParam * u, v)

                         // Weakly imposing the tangential component
                         // - innerProduct(2*mu*Eps(u)*n*t, v*t)  //natural
                         // + innerProduct(p, v*n)  // natural
                         //
                         // + innerProduct(u*t, 2*mu*Eps(v)*t*n)
                         // + innerProduct(penaltyParam*u*t, v*t)

                         // - innerProduct(u*n, q)  // essential
                         ,
                         Khi, INTEGRAL_BOUNDARY);

      stokes.addLinear(+innerProduct(gh.expression(2), penaltyParam * v) +
                           innerProduct(gh.expression(2), 2 * mu * Eps(v) * n)

                       // impose tangential component weakly
                       // + innerProduct(gh*t, penaltyParam*v*t)
                       // + innerProduct(gh*t, 2*mu*Eps(v)*t*n)

                       // - innerProduct(gh.expression(2), q*n)
                       ,
                       Khi, INTEGRAL_BOUNDARY);

      FunTest grad2un = grad(grad(u) * n) * n;
      stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
          innerProduct(uPenParam * pow(hi, 1) * jump(u),
                       jump(v)) // [Method 1: Remove jump in vel]
              + innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n),
                             jump(grad(v) * n)) +
              innerProduct(uPenParam * pow(hi, 5) * jump(grad2un),
                           jump(grad2un))
              // //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
              // //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)),
              // jump(grad(q)))
              // //
              // //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), mu*jump(v))
              // // [Method 1: Remove jump in vel]
              // //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n),
              // mu*jump(grad(v)*n))
              // //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un),
              // jump(grad2un))
              // //
              - innerProduct(pPenParam * hi * jump(p), jump(div(v))) +
              innerProduct(pPenParam * hi * jump(div(u)), jump(q)) -
              innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)),
                           jump(grad(div(v)))) +
              innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))),
                           jump(grad(q)))
          // //
          // //   // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)),
          // mu*jump(div(v))) // [Method 1: Remove jump in vel]
          // //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))),
          // mu*jump(grad(div(v))))
          ,
          Khi
          // //   // , macro
      );
      std::cout << "Time assembly " << MPIcf::Wtime() - t_assembly_begin
                << std::endl;

      // matlab::Export(stokes.mat_, "mat1.dat");

      // // // Sets uniqueness of the pressure
      double tt0 = MPIcf::Wtime();
      // stokes.setDirichlet(ghRT, Khi);

      int N      = stokes.get_nb_dof();
      R2 P       = Ph[0].Pt(0);
      int dfp    = Ph[0].loc2glb(0) + Vh.get_nb_dof();
      double val = fun_exact_p(P, 0, 0);
      // // std::cout << Ph[0].PtHat(0) << std::endl;
      // // std::cout << P << "\t" << val << std::endl;
      // // df2rm.insert({dfu, val});
      // // stokes.mat_[std::make_pair(dfu,dfp)] += 1e15;
      // // stokes.mat_[std::make_pair(Vh.get_nb_dof(),Vh.get_nb_dof())] +=
      // 1e15;
      // // stokes.rhs_(Vh.get_nb_dof()) = -10;
      // // stokes.mat_[std::make_pair(dfp,dfp)] += 1e15;
      // // stokes.rhs_(dfu) = val;
      // eraseAndSetRow(N, stokes.mat_, stokes.rhs_, dfp, dfp, val);
      //   matlab::Export(stokes.mat_, "mat0.dat");
      // return 0;
      // // std::cout << " rm dof \t" << dfp << std::endl;
      // // std::cout << Vh.get_nb_dof()  << std::endl;
      std::cout << "Time setting boundary " << MPIcf::Wtime() - tt0
                << std::endl;
      // getchar();

      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,p1), 0., Khi
      //   // , INTEGRAL_BOUNDARY
      //   // , INTEGRAL_EXTENSION
      //   // , 1.
      // );

      // Fun_h fun1(Wh, fun_1);
      // double areaBubble   = integral(Khi, fun1, 0, 1, 0) ;
      // int nndf = stokes.rhs_.size();
      // stokes.mat_[std::make_pair(N,N)] = -4+areaBubble;
      // stokes.mat_[std::make_pair(N,N+1)] = 1;
      // stokes.rhs_.resize(nndf+1);
      // stokes.rhs_(nndf) = 1.;
      double t_solving = MPIcf::Wtime();
      stokes.solve();
      std::cout << "Time solver " << MPIcf::Wtime() - t_solving << std::endl;

      // std::cout << stokes.rhs_(dfp) << std::endl;
      // nx *= 2;
      // ny *= 2;
      // continue;
      // ------------------------------------------------
      // ------------------------------------------------
      // POST PROCESS SOLUTION
      // Lagrange multiplier
      int ndf    = Ph.get_nb_dof() + Vh.get_nb_dof();
      int idx0_s = Vh.get_nb_dof();
      // R data_lagrange = stokes.rhs_[ndf];
      // std::cout << " lambda = " << data_lagrange << std::endl;
      // std::cout << " gamma = " << stokes.rhs_[ndf+1] << std::endl;

      // Exrtact solution
      Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
      Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Ph, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      // std::cout << "mean P \t" << integral(Khi,ph, 0, 0, 0)/(4-areaBubble) <<
      // std::endl;

      {
         Fun_h soluErr(Vh, fun_exact_u);
         Fun_h soluh(Vh, fun_exact_u);
         Fun_h solph(Ph, fun_exact_p);

         soluErr.v -= uh.v;
         soluErr.v.map(fabs);
         // Fun_h divSolh(Wh, fun_div);
         // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

         Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
         writer.add(uh, "velocity", 0, 2);
         writer.add(ph, "pressure", 0, 1);
         writer.add(solph, "pressureExact", 0, 1);
         writer.add(fabs(femSol_0dx + femSol_1dy), "divergence");
         writer.add(soluh, "velocityExact", 0, 2);
         writer.add(soluErr, "velocityError", 0, 2);
      }

      R errU       = L2normCut(uh, fun_exact_u, 0, 2);
      R errP       = L2normCut(ph, fun_exact_p, 0, 1);
      R errDiv1    = L2normCut(femSol_0dx + femSol_1dy, Khi, 1);
      R errDiv0    = L2normCut(femSol_0dx + femSol_1dy, Khi, 0);
      R maxErrDiv0 = maxNormCut(femSol_0dx + femSol_1dy, Khi, 0);
      R maxErrDiv1 = maxNormCut(femSol_0dx + femSol_1dy, Khi, 1);

      // // solExVec -= stokesDiv.rhs;
      // for(int i=0;i<solExVec.size();++i){
      //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
      // }
      //
      // Fun_h solEx(mixedSpace, solExVec);
      // writerS.add(solEx,"uh_err", 0, 2);

      // writerS.add(solEx,"uh_ex", 0, 2);
      // writerS.add(solEx,"ph_ex", 2, 1);

      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2_0.push_back(errDiv0);
      divl2_1.push_back(errDiv1);
      divmax0.push_back(maxErrDiv0);
      divmax1.push_back(maxErrDiv1);

      // divmax.push_back(maxErrDiv);
      h.push_back(1. / nx);
      if (i == 0) {
         convu.push_back(0);
         convp.push_back(0);
      } else {
         convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
         convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
      }

      std::cout << "Time iteration " << MPIcf::Wtime() - t_iter_begin
                << std::endl;

      nx *= 2;
      ny *= 2;
   }
   std::cout << "\n"
             << std::left << std::setw(10) << std::setfill(' ') << "h"
             << std::setw(15) << std::setfill(' ') << "err_p u" << std::setw(15)
             << std::setfill(' ') << "conv p" << std::setw(15)
             << std::setfill(' ') << "err u" << std::setw(15)
             << std::setfill(' ') << "conv u" << std::setw(15)
             << std::setfill(' ') << "err divu0" << std::setw(15)
             << std::setfill(' ') << "err divu1" << std::setw(15)
             << std::setfill(' ') << "err max divu0" << std::setw(15)
             << std::setfill(' ')
             << "err max divu1"
             // << std::setw(15) << std::setfill(' ') << "conv divu"
             // << std::setw(15) << std::setfill(' ') << "err_new divu"
             // << std::setw(15) << std::setfill(' ') << "convLoc divu"
             // << std::setw(15) << std::setfill(' ') << "err maxdivu"
             // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
             << "\n"
             << std::endl;
   for (int i = 0; i < h.size(); ++i) {
      std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i]
                << std::setw(15) << std::setfill(' ') << pl2[i] << std::setw(15)
                << std::setfill(' ') << convp[i] << std::setw(15)
                << std::setfill(' ') << ul2[i] << std::setw(15)
                << std::setfill(' ') << convu[i] << std::setw(15)
                << std::setfill(' ') << divl2_0[i] << std::setw(15)
                << std::setfill(' ')
                << divl2_1[i]
                // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                << std::setw(15) << std::setfill(' ') << divmax0[i]
                << std::setw(15) << std::setfill(' ')
                << divmax1[i]
                // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                << std::endl;
   }
}
#endif

#ifdef STATIC_DROP_EXAMPLE

R2 shift(0, 0);
double sigma = 1;
double mu1   = 1;
double mu2   = 1;

R fun_levelSet(const R2 P, int i) {
   return sqrt((P.x - shift.x) * (P.x - shift.x) +
               (P.y - shift.y) * (P.y - shift.y)) -
          0.50001;
}

int main(int argc, char **argv) {
   typedef TestFunction<2> FunTest;
   typedef FunFEM<Mesh2> Fun_h;
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;

   const double cpubegin = CPUtime();
   MPIcf cfMPI(argc, argv);

   ProblemOption optionStokes;
   optionStokes.order_space_element_quadrature_ = 5;

   vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

   for (int i = 0; i < 3; ++i) { // i<3

      std::cout << "\n ------------------------------------- " << std::endl;
      //
      int nx = 10, ny = nx;
      Mesh Kh(nx, ny, -1., -1., 2., 2.);
      // Mesh Kh0("../mesh/FittedMesh.msh");

      const R hi           = 2. / (nx - 1);
      const R penaltyParam = 1e1 / hi;
      // const R sigma = 1;
      double uPenParam     = 1e0;
      double pPenParam     = 1e0;
      double jumpParam     = 1e0;

      // CutFEMParameter mu(1e-1,1e-3);
      // CutFEMParameter invmu(1e1,1e3);
      CutFEMParameter mu(mu1, mu2);
      // double sigma = 0;//24.5;//700;
      double kappa1 = 0.5; // mu(1)/(mu(0)+mu(1));
      double kappa2 = 0.5; // mu(0)/(mu(0)+mu(1));

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);

      // Space Lh0(Kh0, DataFE<Mesh2>::P1);
      // Fun_h levelSet0(Lh0, fun_levelSet);
      // Paraview<Mesh> writer2(Kh0, "fittedMesh.vtk");
      // writer2.add(levelSet0, "levelSet", 0, 1);
      // return 0;

      Lagrange2 FEvelocity(2);
      Space UUh(Kh, FEvelocity);

      // Space Qh(Kh, DataFE<Mesh>::P1);

      Space Wh(Kh, DataFE<Mesh2>::RT1); // REMOVE KhESE TWO LATER
      // Space Wh(Kh, FEvelocity);
      Space Qh(Kh, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE

      ActiveMesh<Mesh> Khi(Kh, interface);
      // Khi.truncate(interface, 1);
      CutSpace Vh(Khi, Wh);
      CutSpace Ph(Khi, Qh);

      CutFEM<Mesh2> stokes(Vh, optionStokes);
      stokes.add(Ph);

      Normal n;
      Tangent t;
      FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0),
          p1(Ph, 1, 0, 0);

      stokes.addBilinear(
          // contractProduct(mu*grad(u),grad(v))
          contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) +
              innerProduct(div(u), q),
          Khi);

      stokes.addBilinear(
          -innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2),
                        jump(v * t)) +
              innerProduct(jump(u * t),
                           average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
              innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
          Khi, innerFacet);
      stokes.addBilinear(
          innerProduct(jump(u), -2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
              innerProduct(-2 * mu * average(Eps(u) * n, kappa1, kappa2),
                           jump(v)) +
              innerProduct(1. / hi / hi * jump(u), jump(v)) +
              innerProduct(average(p, kappa1, kappa2), jump(v * n))
          // - innerProduct(jump(u*n), average(q,0.5,0.5))
          ,
          interface);
      // stokes.addBilinear(
      //   innerProduct(jump(u), -2*mu*average(grad(v)*n, 0.5,0.5))
      //   + innerProduct(-2*mu*average(grad(u)*n,0.5,0.5), jump(v))
      //   + innerProduct(10*jump(u), jump(v))
      //   + innerProduct(average(p,0.5,0.5), jump(v*n))
      //   - innerProduct(jump(u*n), average(q,0.5,0.5))
      //   , interface
      // );
      stokes.addLinear(innerProduct(2., average(v * n, kappa2, kappa1)) * sigma,
                       interface);

      stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                             + innerProduct(u, 2 * mu * Eps(v) * n) +
                             innerProduct(p, v * n) // natural
                             + innerProduct(penaltyParam * u, v)
                         // - innerProduct(u*n, q)  // essential
                         ,
                         Khi, boundary);

      FunTest grad2un = grad(grad(u) * n) * n;
      stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
          innerProduct(uPenParam * pow(hi, 1) * jump(u),
                       jump(v)) // [Method 1: Remove jump in vel]
              + innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n),
                             jump(grad(v) * n)) +
              innerProduct(uPenParam * pow(hi, 5) * jump(grad2un),
                           jump(grad2un))
              // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
              // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

              // +innerProduct(uPenParam*pow(hi,1)*jump(u), mu*jump(v)) //
              // [Method 1: Remove jump in vel]
              // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n),
              // mu*jump(grad(v)*n))
              // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))

              -
              innerProduct(pPenParam * hi * jump(p), (1. / mu) * jump(div(v))) +
              innerProduct(pPenParam * hi * jump(div(u)), (1. / mu) * jump(q)) -
              innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)),
                           (1. / mu) * jump(grad(div(v)))) +
              innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))),
                           (1. / mu) * jump(grad(q)))

          // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)), mu*jump(div(v))) //
          // [Method 1: Remove jump in vel]
          // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))),
          // mu*jump(grad(div(v))))
          ,
          Khi
          // , macro
      );

      // Sets uniqueness of the pressure
      // double tt0 = MPIcf::Wtime();
      // int N = stokes.get_nb_dof();
      // std::map<int, double> df2rm;
      // R2 P = Qh[0].Pt(0);
      // double val = fun_exact_p(P, 0, 0);
      // df2rm.insert({Vh.get_nb_dof(), val});
      // eraseRow(N, stokes.mat_, stokes.rhs_, df2rm);
      // std::cout << "Time setting p val " << MPIcf::Wtime() - tt0 <<
      // std::endl; stokes.addLagrangeMultiplier(
      //   innerProduct(1.,p1), 0., Khi
      // );

      stokes.solve();

      // EXTRACT SOLUTION
      int idx0_s  = Vh.get_nb_dof();
      Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
      Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Ph, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      {
         // Fun_h soluErr(Vh, fun_exact_u);
         // Fun_h soluh(Vh, fun_exact_u);
         // soluErr.v -= uh.v;
         // soluErr.v.map(fabs);
         // Fun_h divSolh(Wh, fun_div);
         // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

         Paraview<Mesh> writer(Khi, "staticDrop.vtk");
         writer.add(levelSet, "levelSet", 0, 1);
         writer.add(uh, "velocity", 0, 2);
         writer.add(ph, "pressure", 0, 1);
         writer.add(fabs(femSol_0dx + femSol_1dy), "divergence");
         // writer.add(soluh, "velocityExact" , 0, 2);
         // writer.add(soluErr, "velocityError" , 0, 2);
         // writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");

         // writer.add(solh, "velocityError" , 0, 2);

         // writer.add(fabs(femDiv, "divergenceError");
      }

      // R errU      = L2normCut(uh,fun_exact_u,0,2);
      // R errP      = L2normCut(ph,fun_exact_p,0,1);
      // R errDiv1    = L2normCut (femSol_0dx+femSol_1dy,Khi, 1);
      // R errDiv0    = L2normCut (femSol_0dx+femSol_1dy,Khi, 0);
      // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,Khi);

      // ul2.push_back(errU);
      // pl2.push_back(errP);
      // divl2.push_back(errDiv0);
      // divmax.push_back(errDiv1);
      //
      // // divmax.push_back(maxErrDiv);
      // h.push_back(1./nx);
      // if(i==0) {convu.push_back(0); convp.push_back(0);}
      // else {
      //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
      // }

      nx *= 2;
      ny *= 2;
   }
   // std::cout << "\n" << std::left
   // << std::setw(10) << std::setfill(' ') << "h"
   // << std::setw(15) << std::setfill(' ') << "err_p u"
   // << std::setw(15) << std::setfill(' ') << "conv p"
   // << std::setw(15) << std::setfill(' ') << "err u"
   // << std::setw(15) << std::setfill(' ') << "conv u"
   // << std::setw(15) << std::setfill(' ') << "err divu0"
   // << std::setw(15) << std::setfill(' ') << "err divu1"
   //
   // // << std::setw(15) << std::setfill(' ') << "conv divu"
   // // << std::setw(15) << std::setfill(' ') << "err_new divu"
   // // << std::setw(15) << std::setfill(' ') << "convLoc divu"
   // // << std::setw(15) << std::setfill(' ') << "err maxdivu"
   // // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
   // << "\n" << std::endl;
   // for(int i=0;i<h.size();++i) {
   //   std::cout << std::left
   //   << std::setw(10) << std::setfill(' ') << h[i]
   //   << std::setw(15) << std::setfill(' ') << pl2[i]
   //   << std::setw(15) << std::setfill(' ') << convp[i]
   //   << std::setw(15) << std::setfill(' ') << ul2[i]
   //   << std::setw(15) << std::setfill(' ') << convu[i]
   //   << std::setw(15) << std::setfill(' ') << divl2[i]
   //   // << std::setw(15) << std::setfill(' ') << convdivPr[i]
   //   // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
   //   // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
   //   << std::setw(15) << std::setfill(' ') << divmax[i]
   //   // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
   //   << std::endl;
   // }
}
#endif

#ifdef FICTITIOUS_STOKES
double computeMean(const CutFESpaceT2 &Vh, R (*f)(const R2)) {
   typedef typename Mesh2::BorderElement BorderElement;
   typedef typename FESpace2::FElement FElement;
   typedef typename Mesh2::Element Element;
   // typedef typename FElement::QFB QFB;
   typedef typename FElement::QF QF;
   typedef typename QF::QuadraturePoint QuadraturePoint;

   double valMean = 0;
   What_d Fop     = Fwhatd(1);
   const QF &qf(*QF_Simplex<R2>(4));
   // const QFB &qfb(*QF_Simplex<R1>(5));

   for (int k = Vh.Th.first_element(); k < Vh.Th.last_element(); k++) {

      // const BorderElement & BE(Vh.Th.be(ifac)); // The face
      // int ifaceK; // index of face of triangle corresp to edge (0,1,2)
      // const int kb = Vh.Th.BoundaryElement(ifac, ifaceK); // index of element
      // (triangle), ifaceK gets modified inside const int k =
      // Vh.idxElementFromBackMesh(kb, 0);

      const FElement &FK((Vh)[k]);
      // R2 normal = FK.T.N(ifaceK);
      const R meas = FK.T.mesure();

      // int nb_face_onB = 0;
      // for(int i=0;i<Element::nea;++i){
      //   int ib = i;
      //   if(Vh.Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
      // }
      // assert(nb_face_onB > 0);
      // double measOnB = nb_face_onB*meas;
      // const int kv = Vh.idxElementFromBackMesh(kb, 0);

      for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

         QuadraturePoint ip(qf[ipq]); // integration point
         const R2 mip = FK.map(ip);   // mip is in the global edge
         const R Cint = meas * ip.getWeight();

         // FK.BF(Fop,FK.T.toKref(mip), fv);
         double val_fh = f(mip);

         for (int d = 0; d < 2; ++d) {
            for (int i = FK.dfcbegin(d); i < FK.dfcend(d); ++i) {
               valMean += Cint * val_fh; // * fv(i,d,op_id)*normal[d];
            }
         }
      }
   }
   return valMean;
}

namespace Erik_Data_StokesDiv {

// Erik
// R shift = 0.5;
// R interfaceRad = 0.25+1e-8;//2./3; // not exactly 1/4 to avoid interface
// cutting exaclty a vertex R fun_levelSet(const R2 P, const int i) {
//   return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) -
//   interfaceRad;
// }
// R fun_rhs(const R2 P, int i, int dom) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1)
//   return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// OLshanskii
R shift        = 0.5;
R interfaceRad = sqrt(0.25); // 2./3; // not exactly 1/4 to avoid interface
                             // cutting exaclty a vertex
R fun_levelSet(const R2 P, const int i) {
   return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) -
          interfaceRad;
}

R fun_div(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   return 0;
}

R fun_rhs(const R2 P, int i, int dom) {
   // R mu=1;
   R x = P.x;
   R y = P.y;
   if (i == 0)
      return 40 * x * (x * x - y * y) - (16 * y - 8);
   else if (i == 1)
      return (16 * x - 8) - 40 * y * (x * x - y * y);
   // if(i==0)      return  40*x*(x*x - y*y);
   // else if(i==1) return - 40*y*(x*x - y*y);
   else
      return 0;
}
R fun_exact_u(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   if (i == 0)
      return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
   else
      return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}
R fun_exact_p(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   return (10 * (x * x - y * y) *
           (x * x - y * y)); // - 670.171/(pi*interfaceRad*interfaceRad);
}
} // namespace Erik_Data_StokesDiv
using namespace Erik_Data_StokesDiv;

int main(int argc, char **argv) {
   typedef TestFunction<2> FunTest;
   typedef FunFEM<Mesh2> Fun_h;
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;

   const double cpubegin = CPUtime();
   MPIcf cfMPI(argc, argv);

   int nx = 10;
   int ny = 10;
   // int d = 2;

   vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

   for (int i = 0; i < 4; ++i) { // i<3

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Kh(nx, ny, 0., 0., 1., 1.);
      // Kh.info();
      const R hi           = 1. / (nx - 1);
      const R penaltyParam = 4e3 / hi;
      const R sigma        = 1e-1;
      double uPenParam     = 1e0;
      double pPenParam     = 1e0;
      double jumpParam     = 1e0;
      double kappa1        = 0.5; // mu(1)/(mu(0)+mu(1));
      double kappa2        = 0.5; // mu(0)/(mu(0)+mu(1));

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);

      Lagrange2 FEvelocity(3);
      Space UUh(Kh, FEvelocity);
      LagrangeDC2 FErho(1);
      Space SSh(Kh, FErho);
      // Space RRh(Kh, DataFE<Mesh>::P2);

      Space Wh(Kh, DataFE<Mesh2>::RT1); // REMOVE KhESE TWO LATER
      // Space Wh(Kh, FEvelocity);
      Space Qh(Kh, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE
      Space PP2h(Kh, DataFE<Mesh2>::P2); // FOR MIXEDSPACE
      Space Sh(Kh, SSh);                 // FOR MIXEDSPACE

      ActiveMesh<Mesh> Khi(Kh); //, interface);
      Khi.truncate(interface, 1);
      Khi.info();

      CutSpace Vh(Khi, Wh);
      CutSpace Ph(Khi, Qh);
      CutSpace P2h(Khi, PP2h);

      CutSpace Uh(Khi, UUh);
      CutSpace Rh(Khi, Sh);

      Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h gh(Uh, fun_exact_u);
      Fun_h exactph(P2h, fun_exact_p);
      ExpressionFunFEM<Mesh> exactp(exactph, 0, op_id);
      // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

      // ProblemOption optionProblem;
      // optionProblem.order_space_element_quadrature_ = 9;
      // optionProblem.order_space_bord_quadrature_ = 9;
      CutFEM<Mesh2> stokes(Vh);
      stokes.add(Ph);
      stokes.add(Rh);

      Normal n;
      Tangent t;
      /* Syntax:
      FunTest (fem space, #components, place in space)
      */
      FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), r(Rh, 2, 0),
          w(Rh, 2, 0);
      R mu = 0.5;
      // a_h(u,v)
      stokes.addBilinear(
          // contractProduct(mu*grad(u),grad(v))
          contractProduct(2 * mu * Eps(u), Eps(v))
          // - innerProduct(p,div(v))
          // + innerProduct(div(u),q)
          ,
          Khi);

      stokes.addBilinear(-innerProduct(p, div(v)) + innerProduct(div(u), q),
                         Khi, INTEGRAL_EXTENSION, 1.);

      stokes.addLinear(innerProduct(fh.expression(2), v), Khi);

      // [NECESSARY FOR RT / BDM ELEMENTS]
      // stokes.addBilinear(
      //   - innerProduct(average(grad(u*t)*n,0.5,0.5), jump(v*t))
      //   + innerProduct(jump(u*t), average(grad(v*t)*n,0.5,0.5))
      //   + innerProduct(1./hi*(sigma*jump(u*t)), jump(v*t))
      //   , Khi
      //   , INTEGRAL_INNER_EDGE_2D
      // );
      stokes.addBilinear(
          -innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2),
                        jump(v * t)) +
              innerProduct(jump(u * t),
                           average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
              innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
          Khi, INTEGRAL_INNER_EDGE_2D);

      stokes.addBilinear(
          // - innerProduct(grad(u)*n, v)  //natural
          // + innerProduct(p, v*n)  // natural
          +innerProduct(r, v) - innerProduct(u, w)
          // + innerProduct(u, grad(v)*n)
          // + innerProduct(penaltyParam*u, v)

          // - innerProduct(2*mu*Eps(u)*n, v)  //natural
          // + innerProduct(p, v*n)  // natural
          //
          // + innerProduct(u, 2*mu*Eps(v)*n)
          // + innerProduct(penaltyParam*u, v)
          ,
          interface);
      // stokes.addLinear(
      // //   // + innerProduct(gh.expression(2), grad(v)*n)
      // //   + innerProduct(gh.expression(2), 2*mu*Eps(v)*n)
      //   + innerProduct(gh.expression(2), penaltyParam*v)
      //
      //   , interface
      // );

      FunTest grad2un = grad(grad(u) * n) * n;
      // FunTest grad2pn = grad(grad(p)*n)*n;
      // FunTest grad2divun = grad(grad(div(u))*n)*n;
      stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
          // //  innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v)) // [Method
          // 1: Remove jump in vel]
          // // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n),
          // jump(grad(v)*n))
          // // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
          // // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
          // // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))
          //
          // +innerProduct(uPenParam*hi*jump(u), jump(v)) // [Method 1: Remove
          // jump in vel] +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n),
          // jump(grad(v)*n))
          // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
          //
          // -innerProduct(pPenParam*hi*jump(p), jump(div(v)))
          // +innerProduct(pPenParam*hi*jump(div(u)), jump(q))
          // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)),
          // jump(grad(div(v))))
          // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) ,
          // jump(grad(q)))
          //
          // // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn),
          // jump(grad2divun))
          // // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) ,
          // jump(grad2pn))
          //
          //
          +innerProduct(uPenParam * hi * jump(r),
                        jump(w)) // [Method 1: Remove jump in vel]

          ,
          Khi
          // , macro
      );

      // Sets uniqueness of the pressure
      R meanP = integral(Khi, exactp, 0);
      // std::cout << meanP << std::endl;
      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,p), meanP
      //   , Khi
      // );

      // matlab::Export(stokes.mat_, "mat0.dat");

      // std::cout << stokes.get_nb_dof() << std::endl;
      stokes.removeDofForHansbo(Rh);
      // matlab::Export(stokes.mat_, "mat1.dat");
      // getchar();

      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,p), meanP
      //   , Khi
      // );
      // std::cout << integral(Khi,exactp,0) << std::endl;

      stokes.solve();

      // EXTRACT SOLUTION
      int idx0_s  = Vh.get_nb_dof();
      Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
      Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Ph, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      {
         Fun_h soluErr(Vh, fun_exact_u);
         Fun_h soluh(Vh, fun_exact_u);
         soluErr.v -= uh.v;
         soluErr.v.map(fabs);
         // Fun_h divSolh(Wh, fun_div);
         // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

         Paraview<Mesh> writer(Khi,
                               "stokesFictitious_" + to_string(i) + ".vtk");
         writer.add(uh, "velocity", 0, 2);
         writer.add(ph, "pressure", 0, 1);
         writer.add(femSol_0dx + femSol_1dy, "divergence");
         writer.add(soluh, "velocityExact", 0, 2);
         writer.add(soluErr, "velocityError", 0, 2);
         // writer.add(solh, "velocityError" , 0, 2);

         // writer.add(fabs(femDiv, "divergenceError");
      }

      R errU      = L2normCut(uh, fun_exact_u, 0, 2);
      R errP      = L2normCut(ph, fun_exact_p, 0, 1);
      R errDiv    = L2normCut(femSol_0dx + femSol_1dy, fun_div, Khi);
      R maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy, fun_div, Khi);

      // // solExVec -= stokesDiv.rhs;
      // for(int i=0;i<solExVec.size();++i){
      //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
      // }
      //
      // Fun_h solEx(mixedSpace, solExVec);
      // writerS.add(solEx,"uh_err", 0, 2);

      // writerS.add(solEx,"uh_ex", 0, 2);
      // writerS.add(solEx,"ph_ex", 2, 1);

      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2.push_back(errDiv);
      divmax.push_back(maxErrDiv);
      h.push_back(1. / nx);
      if (i == 0) {
         convu.push_back(0);
         convp.push_back(0);
      } else {
         convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
         convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
      }

      nx *= 2;
      ny *= 2;
   }
   std::cout << "\n"
             << std::left << std::setw(10) << std::setfill(' ') << "h"
             << std::setw(15) << std::setfill(' ')
             << "err_p"
             // << std::setw(15) << std::setfill(' ') << "conv p"
             << std::setw(15) << std::setfill(' ')
             << "err u"
             // << std::setw(15) << std::setfill(' ') << "conv u"
             << std::setw(15) << std::setfill(' ')
             << "err divu"
             // << std::setw(15) << std::setfill(' ') << "conv divu"
             // << std::setw(15) << std::setfill(' ') << "err_new divu"
             // << std::setw(15) << std::setfill(' ') << "convLoc divu"
             << std::setw(15) << std::setfill(' ')
             << "err maxdivu"
             // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
             << "\n"
             << std::endl;
   for (int i = 0; i < h.size(); ++i) {
      std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i]
                << std::setw(15) << std::setfill(' ')
                << pl2[i]
                // << std::setw(15) << std::setfill(' ') << convpPr[i]
                << std::setw(15) << std::setfill(' ')
                << ul2[i]
                // << std::setw(15) << std::setfill(' ') << convuPr[i]
                << std::setw(15) << std::setfill(' ')
                << divl2[i]
                // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                << std::setw(15) << std::setfill(' ')
                << divmax[i]
                // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                << std::endl;
   }
}
#endif

#ifdef POISSEUILLE_EXAMPLE

namespace Erik_Data_StokesDiv {
double mu1 = 1;
double mu2 = 1;
R2 shift(0., 0.);
double sigma = 1;
R fun_levelSet(const R2 P, int i) {
   return sqrt((P.x - shift.x) * (P.x - shift.x) +
               (P.y - shift.y) * (P.y - shift.y)) -
          1.00001;
}

R fun_exact_p(const R2 P, int i, int dom) {
   R x = P.x;
   R y = P.y;
   return 0.;
}
} // namespace Erik_Data_StokesDiv
using namespace Erik_Data_StokesDiv;

int main(int argc, char **argv) {
   typedef TestFunction<2> FunTest;
   typedef FunFEM<Mesh2> Fun_h;
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;

   const double cpubegin = CPUtime();
   MPIcf cfMPI(argc, argv);

   int nx = 20;
   int ny = 20;
   // int d = 2;

   ProblemOption optionStokes;
   optionStokes.order_space_element_quadrature_ = 5;

   vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

   for (int i = 0; i < 3; ++i) { // i<3

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Kh(nx, ny, -3., -2., 6., 4.);
      const R hi           = 2. / (nx - 1);
      const R penaltyParam = 1e1 / hi;
      // const R sigma = 1;
      double uPenParam     = 1e0;
      double pPenParam     = 1e0;
      double jumpParam     = 1e0;

      // CutFEMParameter mu(1e-1,1e-3);
      // CutFEMParameter invmu(1e1,1e3);
      CutFEMParameter mu(mu1, mu2);
      // double sigma = 0;//24.5;//700;
      double kappa1 = 0.5; // mu(1)/(mu(0)+mu(1));
      double kappa2 = 0.5; // mu(0)/(mu(0)+mu(1));

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);

      Lagrange2 FEvelocity(2);
      Space UUh(Kh, FEvelocity);

      // Space Qh(Kh, DataFE<Mesh>::P1);

      Space Wh(Kh, DataFE<Mesh2>::RT1); // REMOVE KhESE TWO LATER
      // Space Wh(Kh, FEvelocity);
      Space Qh(Kh, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE

      ActiveMesh<Mesh> Khi(Kh, interface);

      ActiveMesh<Mesh> Kh0(Kh);
      Kh0.createSurfaceMesh(interface);
      ActiveMesh<Mesh> Kh1(Kh);
      Kh1.truncate(interface, -1);
      ActiveMesh<Mesh> Kh2(Kh);
      Kh2.truncate(interface, 1);

      // Khi.truncate(interface, 1);
      CutSpace Vh(Khi, Wh);
      CutSpace Ph(Khi, Qh);
      CutSpace Uh(Khi, UUh);

      // Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      // Fun_h gh(Uh, fun_exact_u);
      // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

      CutFEM<Mesh2> stokes(Vh, optionStokes);
      stokes.add(Ph);

      Normal n;
      Tangent t;
      FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0),
          p1(Ph, 1, 0, 0);

      stokes.addBilinear(
          // contractProduct(mu*grad(u),grad(v))
          contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) +
              innerProduct(div(u), q),
          Khi);

      stokes.addBilinear(
          -innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2),
                        jump(v * t)) +
              innerProduct(jump(u * t),
                           average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
              innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
          Khi, INTEGRAL_INNER_EDGE_2D);
      stokes.addBilinear(
          innerProduct(jump(u), -2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
              innerProduct(-2 * mu * average(Eps(u) * n, kappa1, kappa2),
                           jump(v)) +
              innerProduct(1. / hi / hi * jump(u), jump(v)) +
              innerProduct(average(p, kappa1, kappa2), jump(v * n))
          // - innerProduct(jump(u*n), average(q,0.5,0.5))
          ,
          interface);
      // stokes.addBilinear(
      //   innerProduct(jump(u), -2*mu*average(grad(v)*n, 0.5,0.5))
      //   + innerProduct(-2*mu*average(grad(u)*n,0.5,0.5), jump(v))
      //   + innerProduct(10*jump(u), jump(v))
      //   + innerProduct(average(p,0.5,0.5), jump(v*n))
      //   - innerProduct(jump(u*n), average(q,0.5,0.5))
      //   , interface
      // );
      stokes.addLinear(innerProduct(3. / 2, average(v * n, kappa2, kappa1)) *
                           sigma,
                       interface);

      stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                             + innerProduct(u, 2 * mu * Eps(v) * n) +
                             innerProduct(p, v * n) // natural
                             + innerProduct(penaltyParam * u, v)
                         // - innerProduct(u*n, q)  // essential
                         ,
                         Khi, INTEGRAL_BOUNDARY, {1, 3});
      stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                         // + innerProduct(u, 2*mu*Eps(v)*n)
                         // + innerProduct(p, v*n)  // natural
                         // + innerProduct(penaltyParam*u, v)
                         // - innerProduct(u*n, q)  // essential
                         ,
                         Khi, INTEGRAL_BOUNDARY, {4});
      // stokes.addBilinear(
      //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
      //   // + innerProduct(u, 2*mu*Eps(v)*n)
      //   // + innerProduct(p, v*n)  // natural
      //   // + innerProduct(penaltyParam*u, v)
      //   // - innerProduct(u*n, q)  // essential
      //   , Khi
      //   , INTEGRAL_BOUNDARY
      //   , {2}
      // );
      stokes.addLinear(innerProduct(1, v * n) // natural
                       ,
                       Khi, INTEGRAL_BOUNDARY, {4});
      // stokes.addLinear(
      //   innerProduct(-1, v*n)  // natural
      //   , Khi
      //   , INTEGRAL_BOUNDARY
      //   ,{2}
      // );

      // l(v)_Omega
      // stokes.addLinear(
      //   innerProduct(fh.expression(2),v)
      //   , Khi
      // );

      FunTest grad2un = grad(grad(u) * n) * n;
      stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
          innerProduct(uPenParam * pow(hi, 1) * jump(u),
                       jump(v)) // [Method 1: Remove jump in vel]
              + innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n),
                             jump(grad(v) * n)) +
              innerProduct(uPenParam * pow(hi, 5) * jump(grad2un),
                           jump(grad2un))
              // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
              // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

              // +innerProduct(uPenParam*pow(hi,1)*jump(u), mu*jump(v)) //
              // [Method 1: Remove jump in vel]
              // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n),
              // mu*jump(grad(v)*n))
              // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))

              - innerProduct(pPenParam * hi * jump(p), jump(div(v))) +
              innerProduct(pPenParam * hi * jump(div(u)), jump(q))
          // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)),
          // (1./mu)*jump(grad(div(v))))
          // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) ,
          // (1./mu)*jump(grad(q)))

          // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)), mu*jump(div(v))) //
          // [Method 1: Remove jump in vel]
          // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))),
          // mu*jump(grad(div(v))))
          ,
          Khi
          // , macro
      );

      // Sets uniqueness of the pressure
      // double tt0 = MPIcf::Wtime();
      // int N = stokes.get_nb_dof();
      // std::map<int, double> df2rm;
      // R2 P = Qh[0].Pt(0);
      // double val = fun_exact_p(P, 0, 0);
      // df2rm.insert({Vh.get_nb_dof(), val});
      //
      // eraseRow(N, stokes.mat_, stokes.rhs_, df2rm);
      // std::cout << "Time setting p val " << MPIcf::Wtime() - tt0 <<
      // std::endl; stokes.addLagrangeMultiplier(
      //   innerProduct(1.,p), 0., Khi
      // );

      stokes.solve();

      // EXTRACT SOLUTION
      int idx0_s  = Vh.get_nb_dof();
      Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
      Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Ph, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      {
         // Fun_h soluErr(Vh, fun_exact_u);
         // Fun_h soluh(Vh, fun_exact_u);
         // soluErr.v -= uh.v;
         // soluErr.v.map(fabs);
         // Fun_h divSolh(Wh, fun_div);
         // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

         Paraview<Mesh> writer(Khi, "poisseuilleSp_" + to_string(i) + ".vtk");
         writer.add(uh, "velocity", 0, 2);
         writer.add(ph, "pressure", 0, 1);
         writer.add(fabs(femSol_0dx + femSol_1dy), "divergence");
         // writer.add(soluh, "velocityExact" , 0, 2);
         // writer.add(soluErr, "velocityError" , 0, 2);

         // writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");

         // writer.add(solh, "velocityError" , 0, 2);

         // writer.add(fabs(femDiv, "divergenceError");
      }

      R errU      = 0.; // L2normCut(uh,fun_exact_u,0,2);
      R errP      = 0.; // L2normCut(ph,fun_exact_p,0,1);
      R errDiv1   = L2normCut(femSol_0dx + femSol_1dy, Khi, 1);
      R errDiv0   = L2normCut(femSol_0dx + femSol_1dy, Khi, 0);
      R maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy, Khi);

      // // solExVec -= stokesDiv.rhs;
      // for(int i=0;i<solExVec.size();++i){
      //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
      // }
      //
      // Fun_h solEx(mixedSpace, solExVec);
      // writerS.add(solEx,"uh_err", 0, 2);

      // writerS.add(solEx,"uh_ex", 0, 2);
      // writerS.add(solEx,"ph_ex", 2, 1);

      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2.push_back(errDiv0);
      divmax.push_back(errDiv1);

      // divmax.push_back(maxErrDiv);
      h.push_back(1. / nx);
      if (i == 0) {
         convu.push_back(0);
         convp.push_back(0);
      } else {
         convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
         convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
      }

      nx *= 2;
      ny *= 2;
   }
   std::cout << "\n"
             << std::left << std::setw(10) << std::setfill(' ') << "h"
             << std::setw(15) << std::setfill(' ') << "err_p u" << std::setw(15)
             << std::setfill(' ') << "conv p" << std::setw(15)
             << std::setfill(' ') << "err u" << std::setw(15)
             << std::setfill(' ') << "conv u" << std::setw(15)
             << std::setfill(' ') << "err divu0" << std::setw(15)
             << std::setfill(' ')
             << "err divu1"

             // << std::setw(15) << std::setfill(' ') << "conv divu"
             // << std::setw(15) << std::setfill(' ') << "err_new divu"
             // << std::setw(15) << std::setfill(' ') << "convLoc divu"
             // << std::setw(15) << std::setfill(' ') << "err maxdivu"
             // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
             << "\n"
             << std::endl;
   for (int i = 0; i < h.size(); ++i) {
      std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i]
                << std::setw(15) << std::setfill(' ') << pl2[i] << std::setw(15)
                << std::setfill(' ') << convp[i] << std::setw(15)
                << std::setfill(' ') << ul2[i] << std::setw(15)
                << std::setfill(' ') << convu[i] << std::setw(15)
                << std::setfill(' ')
                << divl2[i]
                // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                << std::setw(15) << std::setfill(' ')
                << divmax[i]
                // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                << std::endl;
   }
}
#endif
