#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#  include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "baseCutProblem.hpp"
#include "../num/gnuplot.hpp"

double fun_boundary(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
}
double fun_rhs(const R2 P, int elementComp) {
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
}
double fun_exact(const R2 P, int elementComp) {
  return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
}
double fun_levelSet(const R2 P, const int i) {
  return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;
}


int main(int argc, char** argv )
{

  // 2 DIMENSIONAL PROBLEM
  // =====================================================
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  // INITIALIZE MPI
  // =====================================================
  MPIcf cfMPI(argc,argv);


  // MESH AND PROBLEM PARAMETERS
  // =====================================================
  int nx = 20, ny = 20;
  int lx = 1., ly = 1.;
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& h(Parameter::h);
  CutFEM_Parameter mu("mu",1.,1.);
  const R meshSize =  lx / nx;

  // CONSTRUCTION OF THE MESH
  // =====================================================
  Mesh Th(nx, ny, 0., 0., lx, ly);


  // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
  // =====================================================
  FESpace Lh(Th, DataFE<Mesh>::P1);
  Fun_h levelSet(Lh, fun_levelSet);
  Interface2 interface(Th, levelSet.v);


  // CONSTRUCTION OF THE FE SPACE AND THE CUT SPACE
  // =====================================================
  FESpace Vh(Th, DataFE<Mesh>::P1dc);
  CutFESpace2 Wh(Vh, interface, {1,-1});


  Mesh2 cutTh(interface);
  FESpace Sh(cutTh, DataFE<Mesh>::P1dc);
  Sh.backSpace = &Vh;
  // Mesh2 cutTh(interface);
  // gnuplot::save(Th, "Th.dat");
  // gnuplot::save(cutTh, "cutTh.dat");
  // gnuplot::save(interface, "interface.dat");
  // gnuplot::save(Wh,0,"Vh1.dat");
  // gnuplot::save(Wh,1,"Vh2.dat");

  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh> poisson(Wh);
  poisson.add(Sh);

  Normal n;
  Fun_h fh(Wh, fun_rhs);
  Fun_h gh(Wh, fun_boundary);
  FunTest u(Wh,1), v(Wh,1);

  FunTest u1(Wh,1,0,0), u2(Wh,1,0,1); // to specify subdomain. Given in the fourth input
  FunTest u0(Sh,1);
  FunTest Dun = (grad(u)*n), Dvn = (grad(v)*n);


  std::cout << jump(u1,u0) << std::endl;
  std::cout << jump(u2,u0) << std::endl;
  poisson.addBilinear(
        innerProduct(jump(u1,u0), u0)
      + innerProduct(jump(u2,u0), u0)
    , interface
  );
  // you could do later when you know what are the kappas
  // innerProduct(jump(K1*u1,K0*u0), u0)
  // innerProduct(jump(K2*u2,K0*u0), u0)
  getchar();

  // ASSEMBLY OF THE LINEAR SYSTEM
  // =====================================================
  // poisson.addBilinear(
  //   innerProduct(mu*grad(u), grad(v))
  // );
  // poisson.addLinear(
  //   innerProduct(fh.expression(1), v)
  // );
  // poisson.addBilinear(
  //   - innerProduct(mu*average1(Dun), jump(v))
  //   - innerProduct(mu*jump(u), average1(Dvn))
  //   + innerProduct(lambdaG*jump(u), jump(v))
  //   , interface
  // );
  // poisson.addBilinearFormBorder(
  //   innerProduct(lambdaB*u, v)
  // );
  // poisson.addLinearFormBorder(
  //   innerProduct(gh.expression(1), lambdaB*v)
  // );

  // poisson.addFaceStabilization(
  //   innerProduct(1e-2*h*jump(grad(u)*n), jump(grad(v)*n))
  // );

  // poisson.addEdgeIntegral(
  //   // - innerProduct(average(beta*u*n), jump(v))
  //   // - innerProduct(0.5*lambdaE*jump(u)  , jump(v))
  // );

  // RESOLUTION OF THE LINEAR SYSTEM
  // =====================================================
  poisson.solve();


  // COMPUTE THE L2 ERROR
  // =====================================================
  Fun_h femSolh(Wh,  poisson.rhs);
  R errU = L2norm(femSolh, fun_exact, 0,1);


  // PRINT THE SOLUTION TO PARAVIEW
  // =====================================================
  Paraview2 writer(Wh, levelSet, "poissonCut.vtk");
  writer.add(femSolh, "poisson", 0, 1);

}
