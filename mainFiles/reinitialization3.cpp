#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "baseCutProblem.hpp"
#include "curvature.hpp"
#include "reinitialization.hpp"
#include "../time_stuff.hpp"
#include "../num/gnuplot.hpp"
#include "projection.hpp"


namespace TestReinitilization3D {

  // R fun_levelSet(const R3 P, const int i) {
  //   return sqrt(P.x*P.x + P.y*P.y+ P.z*P.z) - 1;
  // }
  R fun_velocity(const R3 P, int ci) {
    R3 R(0.,P.z, 0.); return R[ci];
  }

  R fun_levelSet1(const R3 P, const int i) { return sqrt(P.x*P.x + (P.y-1.4)*(P.y-1.4) + (P.z+0.8)*(P.z+0.8)) - 1.;}
  R fun_levelSet2(const R3 P, const int i) { return sqrt(P.x*P.x + (P.y+1.4)*(P.y+1.4) + (P.z-0.8)*(P.z-0.8)) - 1.;}

  R fun_levelSet(const R3 P, const int i) {
    if(2.8*P.y - 1.6*P.z > 0) return fun_levelSet1(P,i);
    else return fun_levelSet2(P,i);
  }

}
using namespace TestReinitilization3D;


int main(int argc, char** argv )
{
  const int d = 3;
  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface3 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef LevelSet3 LevelSet;
  typedef GCurvature<Mesh> Curvature;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  std::string pathOutpuFolder = "../../outputFiles/reinitialization3/";
  std::string pathOutpuFigure = "../../outputFiles/reinitialization3/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"data.dat", std::ofstream::out);


  int nx = 20;
  double dt = 0.02;
  // Mesh Th(nx, 3*nx, nx, -2., -6., -2, 4.,12., 4.);
  Mesh Th("../mesh/doubleDropMesh3D_350000T");

  Lagrange3 FEvelocity(2);
  FESpace Uh(Th, FEvelocity);
  Fun_h velocity(Uh, fun_velocity);

  FESpace Lh(Th, DataFE<Mesh>::P1);
  Fun_h ls(Lh, fun_levelSet);

  CReinitialization<Mesh> reinitialization;
  reinitialization.number_iteration = 1;
  reinitialization.epsilon_diffusion = 1e-1;
  reinitialization.dt = dt/2;
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "ON";
  reinitialization.precision_correction = 1e-6;
  reinitialization.max_iteration_correction = 20;
  // reinitialization.info();
  // reinitialization.perform(ls);

  Uh.info();

  if(MPIcf::IamMaster()){
    Paraview3 writer(Lh, pathOutpuFigure+"testLevelSet3D0.vtk");
    writer.add(ls, "levelSet", 0, 1);
  }
  for(int i=0;i<0;++i){
    LevelSet3 levelSet(ls, velocity, velocity, dt);
    ls.init(levelSet.rhs);
    // reinitialization.perform(ls);
    Interface3 interface2(Th, ls.v);
    double area = reinitialization.computeArea(Lh, interface2);
    std::cout << i << "\t volume = " << area << std::endl;
    if(MPIcf::IamMaster()){
      Paraview3 writer(Lh, pathOutpuFigure+"testLevelSet3D"+to_string(i+1)+".vtk");
      writer.add(ls, "levelSet", 0, 1);
    }
  }

  // if(MPIcf::IamMaster()){
  //   Paraview3 writer(Lh, pathOutpuFigure+"levelSet3DTest.vtk");
  //   writer.add(ls, "levelSet", 0, 1);
  // }
  return 0;

  // Interface3 interface(Th, ls.v);
  // double areaInitial = reinitialization.computeArea(Lh, interface);
  // std::cout << " area 1 ->\t" << areaInitial  << "\t" << 4./3*M_PI*0.5*0.5*0.5 << std::endl;
  // ls.v += 0.01;
  // Interface3 interface2(Th, ls.v);
  // double areaWrong = reinitialization.computeArea(Lh, interface2);
  // std::cout << " area 2 ->\t" << areaWrong  << "\t" << 4./3*M_PI*0.5*0.5*0.5 << std::endl;
  //
  // double t0 = MPIcf::Wtime();
  // R delta = reinitialization.correction(ls, areaInitial);
  // std::cout << "Time correction area \t" << MPIcf::Wtime() << std::endl;
  //
  // Interface3 interface3(Th, ls.v);
  // double areaCorrected = reinitialization.computeArea(Lh, interface3);
  // std::cout << " area 3 ->\t" << areaCorrected << "\t" << 4./3*M_PI*0.5*0.5*0.5 << std::endl;
  //
  //
  // reinitialization.perform(ls);

}
