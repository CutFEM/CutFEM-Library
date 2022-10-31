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


namespace Curvature_Donut{
  template<typename A>
  A fun_levelSet(const A x, const A y, const A z) {
    return sqrt(z*z + (sqrt(x*x+y*y)-1.0) * (sqrt(x*x+y*y)-1.0) ) - 0.6;
  }


  R fun_levelSet(const R3 P, const int i) {
    return sqrt(P.z*P.z + (sqrt(P.x*P.x+P.y*P.y)-1.0)*(sqrt(P.x*P.x+P.y*P.y)-1.0) ) - 0.6;
  }
  // R fun_levelSet(const R2 P, const int i) {
  //   return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}

  R3 fun_normal(const R3 P) {
     Diff<R, 3> X(P.x,0), Y(P.y,1), Z(P.z,2) ;
     Diff<R, 3> rho = fun_levelSet(X,Y,Z);
     R3 normal(rho.d[0], rho.d[1], rho.d[2]);
     return normal / (normal.norm());
   }

   R3 fun_meanCurvature(const R3 P) {
     R eps = 0.0;
     R3 normal = fun_normal(P);
     Diff<Diff<R,3>, 3> X(P.x,0), Y(P.y,1), Z(P.z,2);
     X.val.d[0]=1;                                      // init val of dx
     Y.val.d[1]=1;
     Z.val.d[2]=1;
     Diff<Diff<R,3>, 3> phi = fun_levelSet(X,Y,Z);
     Diff<R,3> den = sqrt(phi.d[0]*phi.d[0]  +  phi.d[1]*phi.d[1] +  phi.d[2]*phi.d[2]);
     Diff<R,3> dphiX, dphiY, dphiZ;
     dphiX = phi.d[0] / den;
     dphiY = phi.d[1] / den;
     dphiZ = phi.d[2] / den;
     R lap = dphiX.d[0] + dphiY.d[1] + dphiZ.d[2];

     return -lap*normal;
   }

}

using namespace Curvature_Donut;


int main(int argc, char** argv )
{
  const int d =3;
  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface3 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GCurvature<Mesh3> Curvature;

  MPIcf cfMPI(argc,argv);

  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;
  int nz = 10;

  vector<double> ul2,uh1,h, convul2, convuh1;

  int order_ls=1, order_curv=1;
  // std::cout << "choose order element for curvature \t ";
  // std::cin >> order_curv;
  // std::cout << "choose order element for levelSet \t ";
  // std::cin >> order_ls;
  //
  // order_ls = MPIcf::Bcast(order_ls, MPIcf::Master(), 1);
  // order_curv = MPIcf::Bcast(order_curv, MPIcf::Master(), 1);


  Lagrange3 FEcurvature(order_curv);
  Lagrange3 FElevelSetk(order_ls);

  // GTypeOfFESum<Mesh2> FEcurv(FEcurvature);

  std::string pathOutpuFolder = "../../outputFiles/curvature3/donut/";
  std::string pathOutpuFigure = "../../outputFiles/curvature3/donut/paraview/";
  assert(IsPathExist(pathOutpuFigure));
  std::ofstream outputData(pathOutpuFolder+"dataP"+to_string(order_curv)+"P"+to_string(order_ls)+".dat", std::ofstream::out);


  outputData << " phi(x) = sqrt(P.z*P.z + (sqrt(P.x*P.x+P.y*P.y)-1.0)*(sqrt(P.x*P.x+P.y*P.y)-1.0) ) - 0.6" << std::endl;
  outputData << " Curvature " << order_curv << std::endl;
  outputData << " Mapping & LevelSet "<< order_ls << std::endl;
  outputData << " Curvature constant \t 1e-2" << std::endl;
  outputData << " Curvature h^(2k-2)" << std::endl;

  for(int i=0;i<2;++i) {
    const double t0 = CPUtime();

    Mesh Th(nx,ny,nz,-1.8,-1.8,-0.8,3.6,3.6,1.6);

    // Build LevelSet order k
    FESpace Lh_k(Th, FElevelSetk);
    Fun_h levelSetPk(Lh_k, fun_levelSet);

    // Build LevelSet order 1 for interface
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    // projection(levelSetPk, levelSet);

    Interface interface(Th, levelSet.v);

    // Build cut Space
    Mesh cutTh(interface);
    FESpace Vh(cutTh, FEcurvature);

    // Build mapping
    // Mapping3 mapping(Vh, levelSetPk);
    // mapping.errorInfty(interface, fun_levelSet);
    //
    //
    // ------------  MEAN CURVATURE  -------------
    Curvature curvature(Vh, interface);//, mapping);

    // Fun_h uex(Vh, fun_meanCurvature);
    // R errU = curvature.L2norm(uex.v, 0, 2);
    // ul2.push_back(errU);
    // errU = curvature.H1norm(uex.v, 0, 2);
    // uh1.push_back(errU);
    //
    Fun_h solMC(Vh, curvature.rhs);
    // VTKwriter2 writerLS(solMC, pathOutpuFigure+"curvature_"+to_string(i)+".vtk");
    // writerLS.add("meanCurvature", 0, 2);
    // writerLS.add(levelSet, "levelSet", 0, 1);

    if(cfMPI.IamMaster()){
      Paraview3 writer(Vh, pathOutpuFigure+"curvature"+to_string(i)+".vtk");
      writer.add(solMC, "curvature", 0, 2);
      writer.add(levelSet, "levelSet", 0, 1);
    }
    //
    // h.push_back(1./nx);
    // if(i==0) {convul2.push_back(0);convuh1.push_back(0);}
    // else {
    //   convul2.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
    //   convuh1.push_back( log(uh1[i]/uh1[i-1])/log(h[i]/h[i-1]));
    // }

    nx *= 2;
    ny *= 2;
    nz *= 2;

    std::cout << " Time iteration \t" << CPUtime() - t0 << std::endl;
  }

  // std::cout << std::setprecision(2);
  // std::cout << "\n\n h \t\t errL2 u \t convL2 U \t errH1 u \t convH1 U" << std::endl;
  // for(int i=0;i<h.size();++i) {
  //   std::cout << h[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
  //   << "\t\t" << uh1[i] << "\t\t" << convuh1[i]
  //   << std::endl;
  // }
  //
  // outputData << std::setprecision(2);
  // outputData << "\n\n h \t\t errL2 u \t convL2 U \t errH1 u \t convH1 U" << std::endl;
  // for(int i=0;i<h.size();++i) {
  //   outputData << h[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
  //   << "\t\t" << uh1[i] << "\t\t" << convuh1[i]
  //   << std::endl;
  // }
}
