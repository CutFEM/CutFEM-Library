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


namespace TestReinitilization2D {

  // R fun_levelSet(const R2 P, const int i) {
  //   return (sqrt(sqrt(P.x*P.x + P.y*P.y)) - 0.75)*0.01;
  //   // return pow(P.x*P.x*P.x*P.x + P.y*P.y*P.y*P.y, 1./4) - 0.25;
  // }
  // R fun_velocity(const R2 P, int ci) {
  //   R2 R(P.y,0.); return R[ci];
  // }
  R fun_levelSet1(const R2 P, const int i) { return sqrt((P.x+1)*(P.x+1) + (P.y-0.25)*(P.y-0.25)) - 0.5;}
  R fun_levelSet2(const R2 P, const int i) { return sqrt((P.x-1)*(P.x-1) + (P.y+0.25)*(P.y+0.25)) - 0.5;}

  R fun_levelSet(const R2 P, const int i) {
    if(2*P.x - 0.5*P.y < 0) return fun_levelSet1(P,i)*0.01;
    else return fun_levelSet2(P,i)*0.01;
  }
  R fun_velocity(const R2 P, int ci) {
    R2 R(P.y,0.); return R[ci];
  }
}
using namespace TestReinitilization2D;

void load(string path, const Mesh2& Th, FunFEM<Mesh2>& v) {

  int nt;
  ifstream f(path.c_str());
  if(!f) {cerr << "Load a file to KN<double> " << path << endl; exit(1);}
  cout << " Read On file \"" << path <<"\""<<  endl;
  f >> nt;
  std::cout << nt << "\t" << Th.nt << std::endl;
  assert(nt == Th.nt);
  KN<double> ue(3*nt);
  for (int i=0;i<3*nt;i++) {
    f >> ue[i] ;
    assert(f.good());
  }
  for(int k=0;k<nt;++k){
    for(int e=0;e<3;++e){
      int idx = Th(k,e);
      v.v[idx] = ue(3*k+e);
    }
  }
}


int main(int argc, char** argv )
{
  const int d = 2;
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GLevelSet<Mesh> LevelSet;
  typedef GCurvature<Mesh> Curvature;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  std::string pathOutpuFolder = "/NOBACKUP/frachon/outputFiles/reinitialization2/testShearFlow/h1/";
  std::string pathOutpuFigure = "/NOBACKUP/frachon/outputFiles/reinitialization2/testShearFlow/h1/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  // assert(IsPathExist(pathOutpuFigure));
  std::ofstream outputData(pathOutpuFolder+"data.dat", std::ofstream::out);



  int nx = 120;
  int ny = 40;
  double dt = 1./40;

  // Mesh2 Th(nx, ny, -3., -1., 6., 2.);
  Mesh2 Th("../mesh/doubleDropMesh.msh");

  FESpace2 LLh(Th, DataFE<Mesh2>::P1);
  Fun_h lls(LLh);
  load("../../checkPoint/levelSet.txt", Th, lls);
  // Paraview2 writerr(LLh, pathOutpuFigure+"checkPoint.vtk");
  // writerr.add(lls, "levelSet", 0, 1);

  // return 0;
  //
  //
  // Lagrange2 FEvelocity(2);
  // FESpace Uh(Th, FEvelocity);
  // Fun_h velocity(Uh, fun_velocity);


  // FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  // Fun_h ls(Lh, fun_levelSet);
  // FESpace2 Lh_k(Th, DataFE<Mesh2>::P2);
  // Fun_h ls_k(Lh_k, fun_levelSet);
  // projection(ls_k, ls);

  // Interface2 interface(Th, ls.v);

  // Paraview2 writer(Lh, pathOutpuFigure+"levelSetSQRT_4040.vtk");
  // writer.add(ls_k, "noReinit", 0, 1);

  CReinitialization<Mesh> reinitialization;
  reinitialization.number_iteration = 1;
  reinitialization.epsilon_diffusion = 1e-2;
  reinitialization.dt = dt/8;
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "ON";
  reinitialization.precision_correction = 1e-8;
  reinitialization.max_iteration_correction = 20;
  // reinitialization.info();
  // reinitialization.perform(ls_k,ls);//, interface);
  //
  // writer.add(ls_k, "Reinit", 0, 1);

  // double areaInitial = reinitialization.computeArea(Lh, interface);
  // std::cout << " area 1 ->\t" << areaInitial  << "\t" << M_PI*0.75*0.75 << std::endl;
  // double area = 0;
  // for(int i=1;i<=200;++i){
  //   LevelSet2 levelSet(ls, velocity, velocity, 0.02);
  //   ls.init(levelSet.rhs);
  //   // if(i%2 == 0 && i > 1) {
  //   reinitialization.perform(ls);
  //   //   Interface2 interface2(Th, ls.v);
  //   //    area = reinitialization.computeArea(Lh, interface2);
  //   //
  //   // }
  //   std::cout << " iter \t" << i << "\t area = " << area << std::endl;
  //   // if(MPIcf::IamMaster() && i%10==0){
  //   //   Paraview2 writer(Lh, pathOutpuFigure+"testLevelSet"+to_string(i)+".vtk");
  //   //   writer.add(ls, "noReinit", 0, 1);
  //   // }
  // }
  // if(MPIcf::IamMaster()){
    Paraview2 writer(LLh, pathOutpuFigure+"LevelSet.vtk");
    writer.add(lls, "old", 0, 1);
  // }
  reinitialization.number_iteration = 100;
  reinitialization.perform(lls);
  // if(MPIcf::IamMaster()){
    // Paraview2 writer(Lh, pathOutpuFigure+"newLevelSet.vtk");
    writer.add(lls, "new", 0, 1);
  // }

}
/*

int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;

  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  int nx = 40;
  int ny = 80;
  Mesh2 Th(nx, ny, 0., 0., 1., 2.);

  // Space for the velocity field
  KN<const GTypeOfFE<Mesh2>* > myFEVel(2);
  myFEVel(0) = &DataFE<Mesh2>::P2;
  myFEVel(1) = &DataFE<Mesh2>::P2;
  GTypeOfFESum<Mesh2> FEvel(myFEVel);
  FESpace2 VelVh(Th, FEvel);
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);

  Fun_h vel(VelVh, fun_velocity);
  Fun_h up(Lh, fun_levelSet);

  // TestFunction<2> u(2), v(2);
  // Normal n;
  // Tangent t;
  // std::cout << (u.t()*t) << std::endl;
  // std::cout << grad(u.t()*t) << std::endl;
  // std::cout << grad(u.t()*t).t()*n << std::endl;
  //
  // return 0;


  for(int i=0;i<10;++i) {

    // LevelSet2 levelSet(Lh);//, fun_levelSet);
    LevelSet2 levelSet(up, vel, vel, 0.1);

    // levelSet.solve(up, vel, vel, 0.1);
    up.init(levelSet.rhs);

    if(MPIcf::IamMaster() ){
      Fun2 sol(Lh, levelSet.rhs);
      VTKwriter2 writerS(sol, "levelSet"+to_string(i)+".vtk");
      writerS.add("levelSet", 0, 1);
    }
  }


  const double cpuend = CPUtime();
  std::cout << cpuend - cpubegin << std::endl;
  return 0;
}
*/
