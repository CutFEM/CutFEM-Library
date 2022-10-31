#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseCutProblem.hpp"

#include "../num/gnuplot.hpp"
#include "../util/redirectOutput.hpp"

#define PROBLEM_POISSON_LAGRANGE
// #define PROBLEM_POISSON_MIXED
// #define PROBLEM_CUT_POISSON_LAGRANGE
//#define PROBLEM_CUT_POISSON_MIXED

/*   Different Poisson Problem using different element,
            FEM and CutFEM
*/

#ifdef PROBLEM_POISSON_LAGRANGE

namespace Data_PoissonLagrange {
// The FE is composed of three components, RT0={comp0,comp1} and P0={comp2}
R fun_boundary(const R3 P, int elementComp) {
    return 2*sin(2*pi*P.x)*sin(4*pi*P.y)*P.z;//sin(pi*P.z);
}
R fun_rhs(const R3 P, int elementComp) {
    return 40*pi*pi*sin(2*pi*P.x)*sin(4*pi*P.y)*P.z;//sin(pi*P.z);
}
R fun_exact(const R3 P, int elementComp) {
  return 2*sin(2*pi*P[0])*sin(4*pi*P[1])*P.z;//sin(pi*P[2]);
  }
}

using namespace Data_PoissonLagrange;

int main(int argc, char** argv )
{
  typedef Mesh3 Mesh;
  typedef FESpace3 FESpace;
  typedef TestFunction<3> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionFunFEM<Mesh> Expression;

  MPIcf cfMPI(argc,argv);

  const double cpubegin = CPUtime();

  int nx = 10; // used to be 10
  int ny = 10;
  int nz = 10;
  vector<double> ul2,pl2,h, convu, convp, ndof_v;

  std::string pathOutpuFolder = "../../outputFiles/poisson3";
  std::string pathOutpuFigure = "../../outputFiles/poisson3/paraview/";

  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");

  const bool writeVTKFiles = false;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }

  for(int i=0;i<4;++i) { // i<6 works too, but takes too much time for debugging
    Mesh Th(nx, ny, nz,0., 0., 0., 1., 1., 1.);

    FESpace Vh(Th, DataFE<Mesh>::P1);

    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);

    const R meshSize = Th[0].lenEdge(0);
    const R penaltyParam = 1e1/(meshSize);
    Fun_h fh(Vh, fun_rhs);
    Fun_h gh(Vh, fun_boundary);


    FEM<Mesh> poisson(Vh);

    int d = 3;
    Normal n;
    FunTest u(Vh,1), v(Vh,1);

    // a(u,v)_Domain
    poisson.addBilinear(
      innerProduct(grad(u), grad(v))
    );
    // l(v)_Domain
    Expression fx(fh, 0, op_id);
    poisson.addLinear(
      innerProduct(fx, v)
    );


    Expression gx(gh, 0, op_id);
    // std::cout << gh.v << std::endl;
    // weak boundary condition
    // poisson.addBilinearFormBorder(
    //   innerProduct(lambdaB*u, v)
    // );
    // poisson.addLinearFormBorder(
    //   innerProduct(gx, lambdaB*v) //+ innerProduct(gh, v)*penaltyParam
    // );

    // impose strong boundary condition
    poisson.addStrongBC(
      gx
    );

    std::cout << "Check: Assembly is passed." << std::endl;
    poisson.solve();


    Fun_h femSolh(Vh,  poisson.rhs);
    R errU = L2norm(femSolh, fun_exact, 0,1);

    // Fun_h uex(Vh, fun_exact);
    // uex.v -= poisson.rhs;
    // R errU = L2norm(uex, 0,1);


    if(MPIcf::IamMaster() && writeVTKFiles) {
      Fun_h sol(Vh, poisson.rhs);
      Paraview3 writerS(Vh, pathOutpuFigure+"poissonLagrange_"+to_string(i)+".vtk");
      writerS.add(sol, "poisson", 0, 1);
      writerS.add(gh, "poisson_ex", 0, 1);

    }

    ul2.push_back(errU);
    h.push_back(1./nx);
    ndof_v.push_back(poisson.nDoF);

    if(i==0) {convu.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
  }

  std::cout << "\n \n " << std::endl;
  myCout << "ndof \t h \t err u \t\t convU " << std::endl;
  for(int i=0;i<ul2.size();++i) {
    myCout << ndof_v[i] << "\t" << h[i] << "\t"
  	      << ul2[i] << "\t" << convu[i] << std::endl;
  }
}
#endif
