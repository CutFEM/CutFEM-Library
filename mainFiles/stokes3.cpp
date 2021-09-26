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
#include "baseProblem.hpp"

namespace Stokes3_Data {

  R fun_rhs(const R3 P, int i) { return 0;  }
  R fun_boundary(const R3 P, int i) { return (i==0)?0.5*P.z : 0;}

}

using namespace Stokes3_Data;

int main(int argc, char** argv )
{
  typedef TestFunction<3> FunTest;
  typedef FunFEM<Mesh3> Fun_h;
  const int d = 3;
  int nx=10,ny=10,nz=10;
  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  TaylorHood3 FEstokes;

  std::string meshstr = to_string(nx)+to_string(ny)+to_string(nz);
  std::string pathOutpuFolder = "../../outputFiles/FEMstokes3/shear/";
  std::string pathOutpuFigure = "../../outputFiles/FEMstokes3/shear/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"data.dat", std::ofstream::out);

  Mesh3 Th(nx,ny,nz, -1., -1., -1., 2., 2., 2.);
  FESpace3 Vh(Th, FEstokes);

  Lagrange3 FE_velocity(2);
  FESpace3 Uh(Th, FE_velocity);
  FESpace3 Ph(Th, DataFE<Mesh3>::P1);
  FEM<Mesh3> stokes({&Uh, &Ph});
  long idx0_ph = Uh.NbDoF();

  // Stokes CutFEM problem
  Fun_h fh(Uh, fun_rhs);
  Fun_h gh(Uh, fun_boundary);

  Normal n;
  FunTest u(Uh,d), p(Ph,1), v(Uh,d), q(Ph,1);
  R mu = 1;
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);


  double t0 = CPUtime();

  stokes.addBilinear(
    contractProduct(2*mu*Eps(u),Eps(v))
    - innerProduct(p, div(v)) + innerProduct(div(u), q)
  );
  std::cout << "Time addBilinear \t" << CPUtime() - t0 << std::endl;
return 0;


  // a(u,v)_dOmega
  stokes.addBilinearFormBorder(
    innerProduct(lambdaB*u,v)
    + innerProduct(p, v.t()*n)
    - innerProduct(u.t()*n, q)
    - innerProduct(2.*mu*Eps(u)*n, v)
    - innerProduct(u, 2.*mu*Eps(v)*n)
    // , {1,2,3,4,5,6}
  ) ;

  // // l(v)_Omega
  // // stokes.addLinear
  // //   innerProduct(fh.expression(3),u)
  // // );
  //
  // l(v)_dOmega
  stokes.addLinearFormBorder(
    innerProduct(gh.expression(3), lambdaB*v)
    - innerProduct(gh.expression(3), 2.*mu*Eps(v)*n)
    - innerProduct(gh.expression(3), p*n)
    // , {1,2,3,4,5,6}
  );

  stokes.addLagrangeMultiplier(
    innerProduct(1.,q), 0.
  );

  stokes.solve();

  Fun_h uh(Uh, stokes.rhs);
  Rn_ data_ph(stokes.rhs(SubArray(Ph.NbDoF(),idx0_ph)));
  Fun_h ph(Ph, data_ph);
  // // R errU = L2norm(uh, fun_boundary,0,2);
  // // R errP = L2norm(ph, fun_pressure,0,1);
  //
  if(MPIcf::IamMaster()){
    Fun_h sol(Uh, stokes.rhs);
    Paraview3 writerS(Vh, pathOutpuFigure+"poiseuille"+meshstr+".vtk");
    writerS.add(uh,"velocity", 0, 3);
    writerS.add(ph,"pressure", 0, 1);
  }
  //   // ul2.push_back(errU);
  //   // pl2.push_back(errP);
  //   // h.push_back(1./nx);
  //   // if(i==0) {convu.push_back(0); convp.push_back(0);}
  //   // else {
  //   //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
  //   //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
  //   // }
  //
  // //   nx *= 2;
  // //   ny *= 2;
  // // }
  // // std::cout << "\n\n h \t err u \t\t convU \t\t err p \t\t convP" << std::endl;
  // // for(int i=0;i<ul2.size();++i) {
  // //   std::cout << h[i] << "\t" << ul2[i] << "\t" << convu[i] << "\t\t"
  // //   << pl2[i] << "\t" << convp[i] << std::endl;
  // // }

  outputData.close();
  std::cout << " Time computation \t" << CPUtime() - cpubegin << std::endl;

}
