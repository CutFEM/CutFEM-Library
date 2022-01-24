#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#  include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "baseCutProblem.hpp"


double fun_boundary(const R2 P, int elementComp) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
}
double fun_rhs(const R2 P, int elementComp) {
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
}
double fun_exact(const R2 P, int elementComp) {
  return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
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
  int nx = 100, ny = 100;
  int lx = 1., ly = 1.;
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);


  // CONSTRUCTION OF THE MESH AND FINITE ELEMENT SPACE
  // =====================================================
  Mesh Th(nx, ny, 0., 0., lx, ly);
  FESpace Vh(Th, DataFE<Mesh>::P1);


  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  FEM<Mesh> poisson(Vh);
  Fun_h fh(Vh, fun_rhs);
  Fun_h gh(Vh, fun_boundary);
  FunTest u(Vh,1), v(Vh,1);
  Normal n;


  // ASSEMBLY OF THE LINEAR SYSTEM
  // =====================================================
  poisson.addBilinear(
    innerProduct(grad(u), grad(v))
  );
  poisson.addLinear(
    innerProduct(fh.expression(1), v)
  );
  poisson.addBilinearFormBorder(
    innerProduct(lambdaB*u, v)
  );
  poisson.addLinearFormBorder(
    innerProduct(gh.expression(1), lambdaB*v)
  );



  // RESOLUTION OF THE LINEAR SYSTEM
  // =====================================================
  poisson.solve();


  // COMPUTE THE L2 ERROR
  // =====================================================
  Fun_h femSolh(Vh,  poisson.rhs);
  R errU = L2norm(femSolh, fun_exact, 0,1);


  // PRINT THE SOLUTION TO PARAVIEW
  // =====================================================
  Paraview2 writer(Vh, "poisson.vtk");
  writer.add(femSolh, "poisson", 0, 1);

}
