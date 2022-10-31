#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseCutProblem.hpp"

#include "../num/gnuplot.hpp"

#define PROBLEM_POISSON_MIXED_BDM1
// #define PROBLEM_POISSON_MIXED_ARTIFICIAL_BDM1
// #define PROBLEM_POISSON_MIXED_CUTFEM_BDM1





#ifdef PROBLEM_POISSON_MIXED_BDM1

/*
TODO:
+ Konvergensordning 1 korrekt för RT0/P0
+ Konvergensordning 2 korrekt för RT1/P1dc

*/

namespace Data_PoissonMixed {
// The FE is composed of three components, RT0={comp0,comp1} and P0={comp2}
  R fun_boundary(const R2 P, int elementComp) {
    if (elementComp == 2) return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    else return 0;
  }
  R fun_rhs(const R2 P, int elementComp) {
    if (elementComp == 2) return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
    else return 0;
  }
  R fun_exact(const R2 P, int elementComp) { // First RT0 flux, then P0 potential
    if (elementComp == 0)
      return 4*pi*cos(2*pi*P[0])*sin(4*pi*P[1]);
    else if (elementComp == 1)
      return 8*pi*sin(2*pi*P[0])*cos(4*pi*P[1]);
    else
      return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
}
using namespace Data_PoissonMixed;

int main(int argc, char** argv )
{
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  int nx = 20; // used to be 10
  int ny = 20;
  vector<double> ul2,pl2,h, convu, convp;

  for(int i=0;i<4;++i) { // i<6 works too, but takes too much time for debugging
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::BDM1;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);

    const CutFEM_Parameter& lambdaB(Parameter::lambdaB); // penalty param

    FESpace2 Qh(Th, DataFE<Mesh2>::P0);

    Fun_h fh(Vh, fun_rhs);// interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Vh, fun_boundary);

    FEM<Mesh2> poissonMixed(Vh);

    int d = 2;
    Normal n;
    FunTest p(Vh,d), u(Vh,1,d), q(Vh,d), v(Vh,1,d);
    FunTest q_n = q.t()*n;

    // a(u,v)_Domain
    //     matlab::Export(poissonMixed.mat, "mat.dat");
    poissonMixed.addBilinear(
      innerProduct(p, q)
      - innerProduct(div(p), v)
      + innerProduct(u, div(q))
      // + innerProduct(1e-10*u, v)
    );
    // a(u,v)_dDomain
    // poissonMixed.addBilinearFormBorder(
    //   innerProduct(lambdaB*u, v)
    // );

    // l(v)_Domain
    ExpressionFunFEM<Mesh2> fx(fh, d, op_id);
    poissonMixed.addLinear(
      innerProduct(fx, v)
    );

    // l(v)_dDomain
    ExpressionFunFEM<Mesh2> gx(gh, d, op_id);
    poissonMixed.addLinearFormBorder(
      innerProduct(gx, q_n)
    );


    std::cout << "Check: Assembly is passed." << std::endl;
    poissonMixed.solve();

    Fun_h sol(Vh, poissonMixed.rhs);
    // Paraview2 writerS(Vh, "poissonMixed_"+to_string(i)+".vtk");
    // writerS.add(sol,"flux", 0, 2);
    // writerS.add(sol,"potential", 2, 1);

    Rn solExVec(Vh.NbDoF());
    interpolate(Vh, solExVec, fun_exact);


    R errU = L2norm(sol, fun_exact,0,2);//poissonMixed.L2norm(solExVec,0,2);
    R errP = L2norm(sol, fun_exact,2,1);

    // uex -= poissonMixed.rhs;
    // Fun_h solEx(Vh, solExVec);
    // Paraview2 writerS2(Vh, "poissonMixedEx"+to_string(i)+".vtk");
    // writerS2.add(solEx,"flux", 0, 2);
    // writerS2.add(solEx,"potential", 2, 1);

    pl2.push_back(errP);
    ul2.push_back(errU);
    h.push_back(1./nx);

    if(i==0) {convp.push_back(0);convu.push_back(0);}
    else {
      convp.push_back( log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    }


    nx *= 2;
    ny *= 2;
  }

  std::cout << "\n\n h \t err u \t\t convU \t\t err p \t\t convP" << std::endl;
  for(int i=0;i<ul2.size();++i) {
    std::cout << h[i] << "\t" << ul2[i] << "\t" << convu[i] << "\t\t"
              << pl2[i] << "\t" << convp[i] << std::endl;
  }


  std::cout << "Time computation " << MPIcf::Wtime() - cpubegin << std::endl;
}
#endif


#ifdef PROBLEM_POISSON_MIXED_ARTIFICIAL_BDM1



#endif

#ifdef PROBLEM_POISSON_MIXED_CUTFEM_BDM1
namespace Data_CutPoissonMixed {

  R fun_boundary(const R2 P, int elementComp) {
    if (elementComp == 2)
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    // return 0;
    else
    return 0;
  }
  R fun_rhs(const R2 P, int elementComp) {
    if (elementComp == 2)
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
    // return 1;
    else
    return 0;
  }
  R fun_exact(const R2 P, int elementComp, int dom) {
    if (elementComp < 2)
    return 0;
    else
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;
  }
}

using namespace Data_CutPoissonMixed;

int main(int argc, char** argv )
{
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;
  vector<double> ul2,pl2,h, convu, convp;

  for(int i=0;i<4;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);

    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    Interface2 interface(Th, levelSet.v);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT0;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);

    // CutFEM stuff
    SubDomain2 Vh1(Vh, levelSet, 1);
    SubDomain2 Vh2(Vh, levelSet,-1);

    KN<SubDomain2*> subDomainsVh(2);
    subDomainsVh(0) = &Vh1;  // domain index 0
    subDomainsVh(1) = &Vh2;  // domain index 1
    CutFESpace2 Wh(subDomainsVh, interface);

    CutFEM<Mesh2,Interface2> poissonMixed(Wh, interface);

    CutFEM_Parameter mu("mu",100.,1.);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    // const R h = Th[0].lenEdge(0);

    // We define fh on the cutSpace
    Fun_h fh(Wh, fun_rhs);
    Fun_h gh(Wh, fun_boundary);

    int d = 2;
    Normal n;
    Tangent t;
    FunTest u(Wh,1,d), v(Wh,1,d), sig(Wh,d), tau(Wh,d);
    FunTest sig_n = sig.t()*n;

    // a(u,v)_Omega
    poissonMixed.addBilinearFormDomain(
      innerProduct(sig, tau)
      + innerProduct(mu*u, div(tau))
      - innerProduct(div(sig), mu*v)
    );
    poissonMixed.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
    );

     poissonMixed.addBilinearFormInterface(
       innerProduct(average1(u), jump(mu*tau.t()*n))
       + innerProduct(jump(mu*sig.t()*n), average1(v))
       + innerProduct(lambdaG*jump(u), jump(v))
     );

     // poissonMixed.addFaceStabilization(
     //   innerProduct(1e-2*h*jump(grad(sig).t()*t), mu*jump(grad(tau).t()*t))
     // );
     // poissonMixed.addFaceStabilization(
     //   innerProduct(1e-2*h*h*h*jump(dx(sig)), mu*jump(dx(tau)))
     //   + innerProduct(1e-2*h*h*h*jump(dy(sig)), mu*jump(dy(tau)))
     // );

     ExpressionFunFEM<Mesh2> fx(fh, d, op_id);
     // l(v)_Omega
     poissonMixed.addLinearFormDomain(
       innerProduct(fx, v)
     );
     ExpressionFunFEM<Mesh2> gx(gh, d, op_id);
     poissonMixed.addLinearFormBorder(
       innerProduct(gx, sig_n)
     );

     std::cout << "Check: Assembly is passed." << std::endl;
     poissonMixed.solve();

     Fun_h sol(Wh, poissonMixed.rhs);

     R errU = L2normCut(sol, fun_exact,0,2);//poissonMixed.L2norm(solExVec,0,2);
     R errP = L2normCut(sol, fun_exact,2,1);

     pl2.push_back(errP);
     ul2.push_back(errU);
     h.push_back(1./nx);

     if(i==0) {convp.push_back(0);convu.push_back(0);}
     else {
       convp.push_back( log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
       convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
     }


     // Paraview2 writer(Wh, levelSet, "poissonCutMixed_"+to_string(i)+".vtk");
     // writer.add(sol,"flux", 0, 2);
     // writer.add(sol,"potential", 2, 1);

     nx *= 2;
     ny *= 2;
   }
   std::cout << "\n\n h \t err u \t\t convU \t\t err p \t\t convP" << std::endl;
   for(int i=0;i<ul2.size();++i) {
     std::cout << h[i] << "\t" << ul2[i] << "\t" << convu[i] << "\t\t"
               << pl2[i] << "\t" << convp[i] << std::endl;
   }

 }


#endif

/*
namespace Data_MixedDarcy {
  bool artificialInterface = false; // [Change to true/false depending on artificial interface used/not]

  R interfaceRad = 0.25;
  R shift = 0.5;

  R fun_radius2(const R2 P){
    return (P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift);
  }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - interfaceRad;
  }

  R fun_dirichlet(const R2 P, int compInd) {
    if (compInd == 0)
    // return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    return 0;
    else
    return 0;
  }
  R fun_neumann(const R2 P, int compInd, int dom) {
    if (compInd == 2) {
      R r2 = fun_radius2(P);
      R radius2 = interfaceRad*interfaceRad;
      return r2/(2*radius2)+3./2.;
    }
    else
      return 0;
  }
  R fun_force(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {
    if (compInd == 2) {
      R radius2 = interfaceRad*interfaceRad;
      return -2./radius2;
    }
    else
      return 0;
  }
  R fun_exact(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-0.5)*(X-0.5) + (Y-0.5)*(Y-0.5);
    R radius2 = interfaceRad*interfaceRad;
    Diff<R, 2> val = r2/(2*radius2) + 3./2;
    if (compInd==2) {
      return val.val;
    }
    else {
      return -val.d[compInd];
    }
  }

  R fun_test(const R2 P, int i) {
    return P[i]*P[i];
  }

}

using namespace Data_MixedDarcy;

int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 3;
  int ny = 3;

  vector<double> ul2,pl2,h, convu, convp;

  for(int i=0;i<4;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT0;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);


    FEM<Mesh2> darcyMixed(Vh);
    // if (artificialInterface) {
    //   FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    //   Fun_h levelSet(Lh, fun_levelSet);
    //   Interface2 interface(Th, levelSet.v);
    //   CutFEM<Mesh2,Interface2> darcyMixed(Vh, interface); // [Overwrites darcyMixed]
    // }

    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;
    R mu = 1;// CutFEM_Parameter mu("mu",100.,1.);
    R mu_G = 2./3*interfaceRad;
    // const CutFEM_Parameter& invh(Parameter::invh);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R hei = Th[0].lenEdge(0);
    const R invh = 1./hei;

    Fun_h fv(Vh, fun_force);
    Fun_h fq(Vh, fun_div);
    Fun_h u0(Vh, fun_dirichlet);
    Fun_h p0(Vh, fun_neumann);

    int d = 2;
    Normal n;
    Tangent t;
    FunTest p(Vh,1,d), q(Vh,1,d), u(Vh,d), v(Vh,d);

    // a(u,v)_Omega
    darcyMixed.addBilinearFormDomain(
      innerProduct(mu*u, v)
      - innerProduct(p, div(v)) + innerProduct(div(u), q) // b(p,v)-b(q,u)
    );
    // Only on Gamma_D (vel normal comp)
    // darcyMixed.addBilinearFormBorder(
    //   innerProduct(lambdaB*u.t()*n, invh*v.t()*n)
    //   + innerProduct(p, v.t()*n) - innerProduct(u.t()*n, q)
    // );

    // picks out values and puts in place component 0 (which is tested for in linear form)
    ExpressionFunFEM<Mesh2> ffq(fq, 2, op_id);

    // l(v)_Omega
    darcyMixed.addLinearFormDomain(
      innerProduct(fv, v)
      + innerProduct(ffq, q)
    );

    // ExpressionFunFEM<Mesh2> u0(gh, d, op_id);
    ExpressionFunFEM<Mesh2> pp0(p0, 2, op_id);
    darcyMixed.addLinearFormBorder(
      - innerProduct(pp0, v.t()*n) // Only on Gamma_N (pressure)
      // + innerProduct(u0, invh*(v.t()*n)*lambdaB) // Only on Gamma_D (vel normal comp)
      // - innerProduct(u0, q) // Only on Gamma_D (vel normal comp)
    );

    std::cout << "Check: Assembly is passed." << std::endl;
    darcyMixed.solve();



    // Fun2 sol(Vh, darcyMixed.rhs);
    // VTKwriter2 writer(sol, "darcyMixed_"+to_string(i)+".vtk");
    // writer.add("velocity", 0, 2);
    // writer.add("pressure", 2, 1);
    //
    // Rn solExVec(Vh.NbDoF());
    // interpolate(Vh, solExVec, fun_exact);
    // R errU = darcyMixed.L2norm(solExVec,0,2);
    // R errP = darcyMixed.L2norm(solExVec,2,1);
    //


    Fun_h solEx(Vh, fun_exact);
    Paraview2 writer(Vh, "testBDM1_"+to_string(i)+".vtk");
    writer.add(solEx, "velocity", 0, 2);
    writer.add(solEx, "pressure", 2, 1);
    //
    // pl2.push_back(errP);
    // ul2.push_back(errU);
    // h.push_back(1./nx);
    //
    // if(i==0) {convp.push_back(0);convu.push_back(0);}
    // else {
    //  convp.push_back( log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
    //  convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    // }

    nx *= 2;
    ny *= 2;
   }

   std::cout << "\n\n h \t err p \t\t conv p \t\t err u \t\t conv u" << std::endl;
   std::cout << std::setw(6);
   for(int i=0;i<ul2.size();++i) {
     std::cout << h[i] << "\t\t" << pl2[i] << "\t\t" << convp[i] << "\t\t"
   	      << ul2[i] << "\t\t" << convu[i] << std::endl;
   }

 }

 */
