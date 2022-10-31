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
#include "curvature.hpp"
#include "GenericMapping.hpp"
#include "DA.hpp"
#include "projection.hpp"


namespace Curvature_Data {
  template<typename A>
  A fun_levelSet(const A x, const A y) {
    return  sqrt(x*x + y*y/0.25) - 1.0;}
    // return sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) - 0.25;}

  R fun_levelSet(const R2 P, const int i) {
    return sqrt(P.x*P.x + P.y*P.y/0.25) - 1.0;
  }
  // R fun_levelSet(const R2 P, const int i) {
  //   return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}


  R2 fun_normal(const R2 P) {
    Diff<R, 2> X(P.x,0), Y(P.y,1);
    Diff<R, 2> rho = fun_levelSet(X,Y);
    R2 normal(rho.d[0], rho.d[1]);
    return normal / (normal.norm());
  }


  R fun_meanCurvature(const R2 P, int i, double t) {
    R eps = 0.0;
    R2 normal = fun_normal(P);
    Diff<Diff<R,2>, 2> X(P.x,0), Y(P.y+eps,1);
    X.val.d[0]=1;                                      // init val of dx
    Y.val.d[1]=1;
    Diff<Diff<R,2>, 2> phi = fun_levelSet(X,Y);
    Diff<R,2> den = sqrt(phi.d[0] * phi.d[0]  +  phi.d[1] * phi.d[1] );
    Diff<R,2> dphiX, dphiY;
    dphiX = phi.d[0] / den;
    dphiY = phi.d[1] / den;
    R lap = dphiX.d[0] + dphiY.d[1];

    return -lap*normal[i];
  }

}

using namespace Curvature_Data;


int main(int argc, char** argv )
{
  const int d =2;
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GCurvature<Mesh2> Curvature;

  MPIcf cfMPI(argc,argv);

  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;
  vector<double> ul2,uh1,h, convul2, convuh1;

  int order_ls=1, order_curv=1;
  // std::cout << "choose order element for curvature \t ";
  // std::cin >> order_curv;
  // std::cout << "choose order element for levelSet \t ";
  // std::cin >> order_ls;
  Lagrange2 FEcurvature(order_curv);
  Lagrange2 FElevelSetk(order_ls);

  // GTypeOfFESum<Mesh2> FEcurv(FEcurvature);

  // std::string pathOutpuFolder = "../../outputFiles/curvature/ellipse/";
  // std::string pathOutpuFigure = "../../outputFiles/curvature/ellipse/paraview/";
  // assert(IsPathExist(pathOutpuFigure));
  // std::ofstream outputData(pathOutpuFolder+"dataP"+to_string(order_curv)+"P"+to_string(order_ls)+".dat", std::ofstream::out);


  // outputData << " phi(x) = sqrt(x*x + y*y/0.25) - 1.0" << std::endl;
  // outputData << " Curvature " << order_curv << std::endl;
  // outputData << " Mapping & LevelSet "<< order_ls << std::endl;
  // outputData << " Curvature constant \t 1e-2" << std::endl;
  // outputData << " Curvature h^(2k-2)" << std::endl;

  for(int i=0;i<6;++i) {
    const double t0 = CPUtime();

    Mesh2 Th(nx, ny, -1, -1, 2, 2);  // Ellipse

    // Build LevelSet order k
    FESpace Lh_k(Th, DataFE<Mesh>::P1);
    Fun_h levelSetPk(Lh_k, fun_levelSet);

    // Build LevelSet order 1 for interface
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    projection(levelSetPk, levelSet);

    Interface2 interface(Th, levelSet.v);

    // Build cut Space
    Mesh2 cutTh(interface);
    FESpace2 Vh(cutTh, interface, FEcurvature);

    // Build mapping
    Mapping2 mapping(Vh, levelSetPk);
    mapping.errorInfty(interface, fun_levelSet);

    // ------------  MEAN CURVATURE  -------------
    Curvature curvature(Vh, interface);//, mapping);

    // Fun_h uex(Vh, fun_meanCurvature);
    Fun_h usol(Vh, curvature.rhs);
    R errU = L2normSurf(usol,fun_meanCurvature, 0,0, 2);
    ul2.push_back(errU);
    errU = 0;//curvature.H1norm(uex.v, 0, 2);
    uh1.push_back(errU);

    // Fun_h solMC(Vh, curvature.rhs);
    // Paraview2 writerLS(Vh, "curvature_"+to_string(i)+".vtk");
    //
    // // Paraview2 writerLS(Vh, pathOutpuFigure+"curvature_"+to_string(i)+".vtk");
    // writerLS.add(solMC,"meanCurvature", 0, 2);
    // writerLS.add(levelSet, "levelSet", 0, 1);


    h.push_back(1./nx);
    if(i==0) {convul2.push_back(0);convuh1.push_back(0);}
    else {
      convul2.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      convuh1.push_back( log(uh1[i]/uh1[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
    std::cout << " Time iteration \t" << CPUtime() - t0 << std::endl;
  }

  std::cout << std::setprecision(2);
  std::cout << "\n\n h \t\t errL2 u \t convL2 U \t errH1 u \t convH1 U" << std::endl;
  for(int i=0;i<h.size();++i) {
    std::cout << h[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << uh1[i] << "\t\t" << convuh1[i]
    << std::endl;
  }

  // outputData << std::setprecision(2);
  // outputData << "\n\n h \t\t errL2 u \t convL2 U \t errH1 u \t convH1 U" << std::endl;
  // for(int i=0;i<h.size();++i) {
  //   outputData << h[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
  //   << "\t\t" << uh1[i] << "\t\t" << convuh1[i]
  //   << std::endl;
  // }
}
