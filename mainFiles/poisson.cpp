#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <omp.h>
#include "../util/cputime.h"

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

// #include "finiteElement.hpp"
#include "baseProblem.hpp"

// #include "../num/gnuplot.hpp"

#define PROBLEM_POISSON_LAGRANGE
//#define PROBLEM_POISSON_MIXED
// #define PROBLEM_CUT_POISSON_LAGRANGE
// #define PROBLEM_ARTIFICIALCUT_POISSON_LAGRANGE
//#define PROBLEM_CUT_POISSON_MIXED

/*   Different Poisson Problem using different element,
            FEM and CutFEM
*/

#ifdef PROBLEM_POISSON_LAGRANGE

namespace Data_PoissonLagrange {
// The FE is composed of three components, RT0={comp0,comp1} and P0={comp2}
R fun_boundary(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
}
R fun_rhs(const R2 P, int elementComp) {
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
}
R fun_exact(const R2 P, int elementComp) {
  return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
}

using namespace Data_PoissonLagrange;
int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionFunFEM<Mesh2> Expression;

  MPIcf cfMPI(argc,argv);

  int nx = 10; // used to be 10
  int ny = 10;
  vector<double> ul2,pl2,h, convu, convp, tid;

  const bool writeVTKFiles = true;

  ProblemOption optionPoisson;
  optionPoisson.order_space_element_quadrature_ = 9;
  int order_element = 4;
  for(int i=0;i<4;++i) {
    const double cpubegin = CPUtime();

    Mesh Th(nx, ny, 0., 0., 1., 1.);
    double hi = 1./(nx-1);
    FESpace2 Vh(Th, DataFE<Mesh>::P4);

    // const CutFEMParameter& lambdaB(Parameter::lambdaB);
    const double lambdaB = 1./pow(hi, order_element+1);

    const R meshSize = Th[0].lenEdge(0);
    const R penaltyParam = 1e1/(meshSize);
    Fun_h fh(Vh, fun_rhs);
    Fun_h gh(Vh, fun_boundary);


    FEM<Mesh> poisson(Vh, optionPoisson);

    int d = 2;
    Normal n;
    FunTest u(Vh,1), v(Vh,1);

    double t0 = MPIcf::Wtime();

    // a(u,v)_Domain
    poisson.addBilinear(
      innerProduct(grad(u), grad(v))
      , Th
    );
    // l(v)_Domain
    poisson.addLinear(
      innerProduct(fh.expression(), v)
      , Th
    );

    // weak boundary condition
    poisson.addBilinear(
      innerProduct(lambdaB*u, v)
      , Th
      , boundary
    );

    poisson.addLinear(
      innerProduct(gh.expression(), lambdaB*v) //+ innerProduct(gh, v)*penaltyParam
      , Th
      , boundary
    );
    poisson.solve();


    Fun_h uh(Vh,  poisson.rhs_);
    R errU = L2norm(uh, fun_exact, 0,1);


    // if(MPIcf::IamMaster() && writeVTKFiles) {
    //   Fun_h sol(Vh, poisson.rhs);
    //   Paraview2 writerS(Vh, pathOutpuFigure+"poissonLagrange_"+to_string(i)+".vtk");
    //   writerS.add(sol, "poisson", 0, 1);
    // }

    ul2.push_back(errU);
    h.push_back(1./nx);

    if(i==0) {convu.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;

    tid.push_back(CPUtime() - cpubegin);
  }

  std::cout << "\n \n " << std::endl;
  std::cout << "h \t err u \t\t convU \t\t  time" << std::endl;
  for(int i=0;i<ul2.size();++i) {
    std::cout << h[i] << "\t"
  	      << ul2[i] << "\t" << convu[i] << "\t \t" << tid[i] << std::endl;
  }
}
#endif





#ifdef PROBLEM_POISSON_MIXED

/*
TODO:
+ Konvergensordning 1 korrekt för RT0/P0
+ Konvergensordning 2 korrekt för RT1/P1dc

*/

namespace Data_PoissonMixed {
// The FE is composed of three components, RT0={comp0,comp1} and P0={comp2}
  R fun_boundary(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_rhs(const R2 P, int elementComp) {
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
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

  const double cpubegin = CPUtime();

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

    FESpace2 Qh(Th, DataFE<Mesh2>::P2);

    Fun_h fh(Qh, fun_rhs);// interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Qh, fun_boundary);

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
    // ExpressionFunFEM<Mesh2> fx(fh, d, op_id);
    poissonMixed.addLinear(
      innerProduct(fh.expression(), v)
    );

    // l(v)_dDomain
    // ExpressionFunFEM<Mesh2> gx(gh, d, op_id);
    poissonMixed.addLinearFormBorder(
      innerProduct(gh.expression(), q_n)
    );


    std::cout << "Check: Assembly is passed." << std::endl;
    poissonMixed.solve();

    Fun_h sol(Vh, poissonMixed.rhs);
    Paraview2 writerS(Vh, "poissonMixed_"+to_string(i)+".vtk");
    writerS.add(sol,"flux", 0, 2);
    writerS.add(sol,"potential", 2, 1);

    Fun_h femSolh(Vh,  poissonMixed.rhs);
    R errU = L2norm(femSolh, fun_exact, 0,2);
    R errP = L2norm(femSolh, fun_exact, 2,1);

    // uex -= poissonMixed.rhs;
    // Fun2 solEx(Vh, solExVec);
    // VTKwriter2 writerS2(solEx, "poissonMixedEx"+to_string(i)+".vtk");
    // writerS2.add("flux", 0, 2);
    // writerS2.add("potential", 2, 1);

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
}
#endif


#ifdef PROBLEM_ARTIFICIALCUT_POISSON_LAGRANGE

namespace Erik_Data_CutPoissonLagrange {
  R fun_boundary(const R2 P, int elementComp) {
    // return 0;
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_rhs(const R2 P, int elementComp) {
    // return 1;
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
  }
  R fun_exact(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;
    // return sqrt(P.x*P.x + P.y*P.y) - 2./3;
  }
}

using namespace Erik_Data_CutPoissonLagrange;

int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionFunFEM<Mesh2> Expression;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10; // used to be 10
  int ny = 10;
  vector<double> ul2,h, convu;

  std::string pathOutpuFolder = "../../outputFiles/poisson/";
  std::string pathOutpuFigure = "../../outputFiles/poisson/paraview/";
  const bool writeVTKFiles = true;

  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }


  for(int i=0;i<4;++i) { // i<6 works too, but takes too much time for debugging
    // Mesh Th(nx, ny, -1., -1., 2., 2.);   // [-1,1]*[-1,1]
    Mesh Th(nx, ny, 0., 0., 1., 1.);   // [-1,1]*[-1,1]
    FESpace2 Vh(Th, DataFE<Mesh2>::P1);

    // We define here the levelSet array. In this case we just interpolate the fun_levelSet defined above
    // No need to move it or do anything else with it
    // In order to create the cut we always use a P1 levelSet. Higher order is obtained by using a mapping of higher order,
    // but ALWAYS levelSet P1 to define the interface and the cuts
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);


    // Now we create the interface between the domains. Array of nodes, edges, normals etc
    Interface2 interface(Th, levelSet.v);
    // You can look at it using GNUPLOT
    // command : plot "gammaErik.dat"  w l, "meshErik.dat" w l
    // gnuplot::save(Th, "meshErik.dat");
    // gnuplot::save(interface, "gammaErik.dat");
    //

    // // Now we create the subDomains.
    // We give the background space, the P1 levelSet and the sign consider
    // We always consider Omega1 the outter subDomain and Omega2 the inner Subdomain
    // SubDomain2 Vh1(Vh, levelSet, 1);
    // SubDomain2 Vh2(Vh, levelSet,-1);
    //
    //
    // // Create the cutFEM Space
    // // We put the subDomain adresses in an array and create the CutFEM Space
    // KN<SubDomain2*> subDomainsVh(2);
    // subDomainsVh(0) = &Vh1;  // domain index 0
    // subDomainsVh(1) = &Vh2;  // domain index 1
    // CutFESpace2 Wh(subDomainsVh, interface);

    // We can take a look at them
    // gnuplot::save(Wh, 0, "Vh1Erik.dat");
    // gnuplot::save(Wh, 1, "Vh2Erik.dat");
    // Try first, still in gnuplot
    // 1) plot "gammaErik.dat"  w l, "Vh1Erik.dat" w l
    // 2) plot "gammaErik.dat"  w l, "Vh2Erik.dat" w l
    // You will notice the elements define in both subDomain. Those are the cut elements
    // and there dof are defined twice (one for each subdomain)

    // TEHN WE DO THE CLASSIC THINGS
    CutFEM<Mesh,Interface2> poisson(Vh, interface);
    // Here we need to give the interface with it.
    // This will change soon. It is actually useless but no time to change Now

    // Need to defined mu using this cutFem_parameter because it takes different values
    // Doing so, we can just define one integral and not one on each subDomain (we could do that as well )
    // CutFEM_Parameter mu("mu",1.,1.);
    // double mu = 1;
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R meshSize = Th[0].lenEdge(0);

    // We define fh on the cutSpace
    Fun_h fh(Vh, fun_rhs);
    Fun_h gh(Vh, fun_boundary);

    Normal n;
    FunTest u(Vh,1), v(Vh,1);
    FunTest Dun = (grad(u).t()*n), Dvn = (grad(v).t()*n);

    // a(u,v)_Omega
    poisson.addBilinearFormDomain(
      innerProduct(grad(u), grad(v))
    );
    // l(v)_Omega
    Expression fx(fh, 0, op_id);
    poisson.addLinearFormDomain(
      innerProduct({fx}, v)
    );


    poisson.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
    );

    Expression gx(gh, 0, op_id);
    poisson.addLinearFormBorder(
      innerProduct(gx, lambdaB*v) //+ innerProduct(gh, v)*penaltyParam
    );

    poisson.solve();


    if(MPIcf::IamMaster() && writeVTKFiles) {
      Fun_h sol(Vh, poisson.rhs);
      Paraview2 writerS(Vh, levelSet, pathOutpuFigure+"poissonCutLagrange_"+to_string(i)+".vtk");
      writerS.add(sol, "poisson", 0, 1);
    }

    // ul2.push_back(errU);
    // h.push_back(1./nx);
    //
    // if(i==0) {convu.push_back(0);}
    // else {
    //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    // }

    nx *= 2;
    ny *= 2;
  }

  // std::cout << "\n\n h \t err u \t\t convU " << std::endl;
  // for(int i=0;i<ul2.size();++i) {
  //   std::cout << h[i] << "\t"
  //         << ul2[i] << "\t" << convu[i] << std::endl;
  // }

}
#endif


#ifdef PROBLEM_CUT_POISSON_LAGRANGE

namespace Erik_Data_CutPoissonLagrange {
  R fun_boundary(const R2 P, int elementComp, int dom) {
    // return 0;
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_rhs(const R2 P, int elementComp, int dom) {
    // return 1;
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
  }
  R fun_exact(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;
    // return sqrt(P.x*P.x + P.y*P.y) - 2./3;
  }
  R2 fparam(double t){ return R2(0.5+0.25*cos(t+1./3), 0.5+0.25*sin(t+1./3));}


}

using namespace Erik_Data_CutPoissonLagrange;
#define _USE_MARKER

int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionFunFEM<Mesh2> Expression;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 50; // used to be 10
  int ny = 50;
  vector<double> ul2,h, convu;

  std::string pathOutpuFolder = "../../outputFiles/poisson/";
  std::string pathOutpuFigure = "../../outputFiles/poisson/paraview/";
  const bool writeVTKFiles = true;

  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }

  int nmarker = 50;
  for(int i=0;i<1;++i) {
    Mesh Th(nx, ny, 0., 0., 1., 1.);   // [-1,1]*[-1,1]

    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    // Interface2 tempI(Th, levelSet.v);

    // Marker interface(Th, fparam, 0, 2*M_PI, nmarker);
    // interface.make_interface();

    // // gnuplot::save(, "marker.dat");
    // gnuplot::save(Th, "Th.dat");
    // gnuplot::saveMarker(interface, "interface.dat");
    // gnuplot::save(tempI, "levelSetInterface.dat");
// return 0;
    Interface2 interface(Th, levelSet.v);

    FESpace2 Vh(Th, DataFE<Mesh2>::P1);

    CutFESpace2 Wh(Vh, interface, {1,-1});

    // McDonald macro(Wh, 0.1);
    // macro.make_S();
    //
    // return 0;

    // FESpace2 Wh(Th, interface, DataFE<Mesh2>::P1);

    CutFEM<Mesh2> poisson(Wh);

    CutFEM_Parameter mu("mu",1.,1.);

    // double mu = 1;
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R meshSize = Th[0].lenEdge(0);

    // We define fh on the cutSpace
    Fun_h fh(Wh, fun_rhs);
    Fun_h gh(Wh, fun_boundary);

    Normal n;
    FunTest u(Wh,1), v(Wh,1);
    FunTest Dun = (grad(u).t()*n), Dvn = (grad(v).t()*n);

    // a(u,v)_Omega
    poisson.addBilinear(
      innerProduct(mu*grad(u), grad(v))
    );
    // l(v)_Omega
    poisson.addLinear(
      innerProduct(fh.expression(1), v)
    );

    poisson.addBilinear(
      - innerProduct(mu*average1(Dun), jump(v))
      - innerProduct(mu*jump(Dun), average2(v))
      - innerProduct(mu*jump(u), average1(Dvn))
      + innerProduct(lambdaG*jump(u), jump(v))
      , interface
    );

    poisson.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
    );
    Expression gx(gh, 0, op_id);
    poisson.addLinearFormBorder(
      innerProduct(gx, lambdaB*v) //+ innerProduct(gh, v)*penaltyParam
    );

    poisson.solve();

    if(MPIcf::IamMaster() && writeVTKFiles) {
      Fun_h sol(Wh, poisson.rhs);
      Paraview2 writerS(Wh, levelSet, pathOutpuFigure+"poissonCutLagrange_"+to_string(i)+".vtk");
      writerS.add(sol, "poisson", 0, 1);
    }

    Fun_h femSolh(Wh,  poisson.rhs);
    // R errU = L2normCut(femSolh, fun_exact, 0,1);
    R errU = L2norm(femSolh, fun_exact, 0,1);

    ul2.push_back(errU);
    h.push_back(1./nx);

    if(i==0) {convu.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
    nmarker *= 2;
  }

  std::cout << "\n\n h \t err u \t\t convU " << std::endl;
  for(int i=0;i<ul2.size();++i) {
    std::cout << h[i]   << "\t"
              << ul2[i] << "\t" << convu[i] << std::endl;
  }

}
#endif

#ifdef PROBLEM_CUT_POISSON_MIXED

namespace Erik_Data_CutPoissonMixed {
  // R fun_boundary(const R2 P, int elementComp) {
  //   return 0;
  // }
  // R fun_rhs(const R2 P, int elementComp) {
  //   return (elementComp == 2);
  // }
  // R fun_levelSet(const R2 P, const int i) {
  //   return sqrt(P.x*P.x + P.y*P.y) - 2./3;
  // }
  R fun_boundary(const R2 P, int elementComp) {
    if (elementComp == 2)
    // return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    return 0;
    else
    return 0;
  }
  R fun_rhs(const R2 P, int elementComp) {
    if (elementComp == 2)
    // return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
    return 1;
    else
    return 0;
  }
  R fun_exact(const R2 P, int elementComp) {
    if (elementComp < 2)
    return 0;
    else
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;
  }
}

using namespace Erik_Data_CutPoissonMixed;

int main(int argc, char** argv )
{
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;

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
    const R h = Th[0].lenEdge(0);

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

     Fun2 sol(Wh, poissonMixed.rhs);


     VTKcutWriter2 writer(sol, levelSet, "poissonCutMixed_"+to_string(i)+".vtk");
     writer.add("flux", 0, 2);
     writer.add("potential", 2, 1);

     nx *= 2;
     ny *= 2;
   }

 }
#endif
