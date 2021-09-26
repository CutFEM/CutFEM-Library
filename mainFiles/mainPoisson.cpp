#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "baseProblem.hpp"
#include "levelSet.hpp"

namespace Poisson_Data {

  R fun_boundary(const R2 P, int i=0) {
    return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
  }
  R fun_rhs(const R2 P, int i=0) {
    return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
  }

  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}
  R fun_velocity(const R2 P, const int i) {
    if(i==0) return 0.;
    else return P.x*(1-P.x);
  }

}

using namespace Poisson_Data;



int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;


  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  int nx = 10;
  int ny = 10;
  vector<double> ul2,uh1,hv, convul2, convuh1;

  for(int i=0;i<4;++i) {
    Mesh Th(nx, ny, 0., 0., 1., 1.);
    FESpace Vh(Th, DataFE<Mesh2>::P1);

    const R h = Th[0].lenEdge(0);
    const R penaltyParam = 1e1/(h);
    // Rn fhv, ghv(Vh.NbDoF());
    // interpolate(Vh, fhv, fun_rhs);
    // interpolate(Vh, ghv, fun_boundary);
    // Fun_h gh(Vh, ghv);

    Fun_h fh(Vh, fun_rhs);
    Fun_h gh(Vh, fun_boundary);


    Normal n;
    FunTest u(Vh,1), v(Vh,1);
    FunTest dun = grad(u).t()*n, dvn = dun;

    FEM<Mesh> poisson(Vh);
    poisson.addBilinearFormDomain(
      innerProduct(grad(u), grad(v))
    );
    poisson.addBilinearFormBorder(
      - innerProduct(dun,v) + innerProduct(u,dvn)
      + innerProduct(u,v)*penaltyParam
    ) ;
    poisson.addLinearFormDomain(
      innerProduct(fh,v)
    );
    poisson.addLinearFormBorder(
      innerProduct(gh,v)*penaltyParam
      + innerProduct(gh,dvn)
    );
    poisson.solve();

    // initialization
    // Poisson<Mesh2> poisson(Vh);
    // poisson.fun_boundary = fun_boundary;
    // poisson.fun_source = fun_rhs;
    // poisson.fun_exact = fun_boundary;
    // poisson.solve();

    Rn uex(Vh.NbDoF());
    interpolate(Vh, uex, fun_boundary);
    R errU = poisson.L2norm(uex);
    ul2.push_back(errU);
    errU = poisson.H1norm(uex);
    uh1.push_back(errU);


    // Fun2 sol(Vh, poisson.rhs);
    // VTKwriter2 writerS(sol, "poissonErik"+to_string(i)+".vtk");
    // writerS.add("velocity", 0, 1);

    hv.push_back(1./nx);
    if(i==0) {convul2.push_back(0);convuh1.push_back(0);}
    else {
      convul2.push_back( log(ul2[i]/ul2[i-1])/log(hv[i]/hv[i-1]));
      convuh1.push_back( log(uh1[i]/uh1[i-1])/log(hv[i]/hv[i-1]));
    }


    nx *= 2;
    ny *= 2;
  }
  //
  std::cout << "\n\n h \t errL2 u \t\t convL2 U\t errH1 u \t\t convH1 U" << std::endl;
  for(int i=0;i<hv.size();++i) {
    std::cout << hv[i] << "\t" << ul2[i] << "\t" << convul2[i]
              << "\t" << uh1[i] << "\t" << convuh1[i]
              << std::endl;
  }
}
