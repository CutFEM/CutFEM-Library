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
#include "curvature.hpp"
#include "projection.hpp"


#include "DA.hpp"


// R3 shift(0.5,0.5,0.5);
R3 shift(0.,0.,0.);
R fun_levelSet(const R3 P, int i) {
  return sqrt((P.x-shift.x)*(P.x-shift.x) + (P.y-shift.y)*(P.y-shift.y)+ (P.z-shift.z)*(P.z-shift.z)) - 2./3;
}

// R fun_boundary(const R3 P,const int i) {return (i==0)? 0.1 : 0;}
R fun_boundary(const R3 P, int i) { return (i==0)?0.5*P.z : 0;}


int main(int argc, char** argv )
{
  const int d = 3;

  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface3 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GLevelSet<Mesh> LevelSet;
  typedef GCurvature<Mesh> Curvature;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;
  int nz = 10;

  int order_ls=1, order_curv=1;
  // std::cout << "choose order element for curvature \t ";
  // std::cin >> order_curv;
  Lagrange3 FEcurvature(order_curv);

  std::string pathOutpuFolder = "../../outputFiles/cutStokes3/Shear/";
  std::string pathOutpuFigure = "../../outputFiles/cutStokes3/Shear/paraview/";
  std::cout << " path to the output files : "+pathOutpuFolder << std::endl;
  std::cout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataP"+to_string(order_curv)+"P"+to_string(order_ls)+".dat", std::ofstream::out);

  outputData << " Curvature " << order_curv << std::endl;
  outputData << " Mapping & LevelSet "<< order_ls << std::endl;
  outputData << " Curvature constant \t 1e-2" << std::endl;
  outputData << " Curvature h^(2k-2)" << std::endl;


  CutFEM_Parameter mu("mu",10,1);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",1e-1,1e0);
  const R sigma = 24.5;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG3);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB3);
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  outputData << "\n Parameters of the problem : \n"
            << " mu_1 = " << mu.val1 << " and mu_2 = "<< mu.val2 << " \n"
            << " the surface tension sigma = " << sigma
            << std::endl;

  double meshSize(2./nx);
  // Mesh Th(nx, ny, nz, 0, 0, 0, 1., 1., 1.);
  Mesh3 Th(nx,ny,nz, -1., -1., -1., 2., 2., 2.);

  // levelSet stuff
  FESpace Lh_k(Th, DataFE<Mesh3>::P1);
  Fun_h levelSetPk(Lh_k, fun_levelSet);
  FESpace Lh(Th, DataFE<Mesh>::P1);
  Fun_h levelSet(Lh, fun_levelSet);
  projection(levelSetPk, levelSet);


  Interface3 interface(Th, levelSet.v);
  // ------------  MEAN CURVATURE  -------------
  // Build cut Space
  Mesh3 cutTh(interface);
  FESpace3 cutVh(cutTh, FEcurvature);

  // Build mapping
  // Mapping2 mapping(cutVh, levelSetPk);
  // mapping.errorInfty(interface, fun_levelSet);

  // ------------  MEAN CURVATURE  -------------
  Curvature curvature(cutVh, interface);//, mapping);
  Fun_h meanCurvature(cutVh, curvature.rhs);

  // ------------  CUT STOKES PROBLEM  ---------------
  TaylorHood3 FEstokes;
  FESpace Vh(Th, FEstokes);

  CutFESpace3 Wh(Vh, interface, {1,-1});

  Fun_h gh(Wh, fun_boundary);

  Normal n;
  FunTest u(Wh,d), p(Wh,1,d), v(Wh,d), q(Wh,1,d), p1(Wh,1,d,0);
  FunTest Eun = (Eps(u)*n);

  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh3> stokes(Wh);

  double t0 = CPUtime();
  //a(u,v)_Omega
  stokes.addBilinear(
    contractProduct(2*mu*Eps(u),Eps(v))
    - innerProduct(p, div(v)) + innerProduct(div(u), q)
  );
  std::cout << "Time addBilinear \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();
  // a(u,v)_gamma
  stokes.addBilinear(
    innerProduct(jump(u), -2*mu*average1(Eps(v)*n))
    + innerProduct(-2*mu*average1(Eps(u)*n), jump(v))
    + innerProduct(lambdaG*jump(u), jump(v))
    + innerProduct(average1(p), jump(v.t()*n))
    - innerProduct(jump(u.t()*n), average1(q))
    , interface
  );
  std::cout << "Time addBilinear Interface \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();
  // a(u,v)_dOmega
  stokes.addBilinearFormBorder(
    innerProduct(lambdaB*u,v)
    + innerProduct(p, v.t()*n)
    - innerProduct(u.t()*n, q)
    - innerProduct(2.*mu*Eps(u)*n, v)
    - innerProduct(u, 2.*mu*Eps(v)*n)
    // ,{6}
  ) ;
  std::cout << "Time addBilinear Border \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();
  stokes.addLinear(
    innerProduct(meanCurvature.expression(3), average2(v))*sigma
    , interface
  );
  std::cout << "Time addLinear Curvature \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();
  // l(v)_dOmega
  stokes.addLinearFormBorder(
    innerProduct(gh.expression(3),lambdaB*v)
    - innerProduct(gh.expression(3),2.*mu*Eps(v)*n)
    - innerProduct(gh.expression(3), p*n)
    // ,{6}
  );
  std::cout << "Time addLinear Border \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();

  stokes.addFaceStabilization(
    innerProduct(1e-2*h_E*jump(grad(u)*n), mu*jump(grad(v)*n))
  );
  // R cch = h*h*h;
  R cch = pow(meshSize,3);
  stokes.addFaceStabilization(
    innerProduct(1e-2*cch*jump(grad(p).t()*n), invmu*jump(grad(q).t()*n))
  );

  std::cout << "Time addStabilization \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();

  stokes.addLagrangeMultiplier(
    innerProduct(1.,p1), 0.
    ,0
  );
  std::cout << "Time addLagrange \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();


  stokes.solve();

  std::cout << "Time Solver \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();

  if(MPIcf::IamMaster()){
    Fun_h uh(Wh, stokes.rhs);
    Paraview3 writer(Wh, levelSet, pathOutpuFigure+"cavity.vtk");
    writer.add(levelSet,"levelSet", 0 , 1);
    writer.add(uh, "velocity", 0, 3);
    writer.add(uh, "pressure", 3, 1);
  }

  std::cout << " Total Time \t" << CPUtime() - cpubegin << std::endl;

}
