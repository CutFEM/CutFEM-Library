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

namespace Stokes_Data {

  // R fun_rhs(const R2 P, int i) {
  //   R mu=1;
  //   R x = P.x;
  //   R y = P.y;
  //   if(i==0) return -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1)) + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3));
  //   else if(i==1) return 2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2)) + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2));
  //   else return 0;
  // }
  //
  // R fun_boundary(const R2 P, int i) {
  //   R x = P.x;
  //   R y = P.y;
  //   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
  //   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
  //   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
  // }

  R fun_rhs(const R2 P, int i) {
    R mu=2;
    if(i==0) return 120*P.x*P.y*(1-mu);
    else if(i==1) return 60*(P.x*P.x-P.y*P.y)*(1-mu);
    // else return 0;
  }
  R fun_boundary(const R2 P, int i) {
    if(i==0) return 20*P.x*pow(P.y,3);
    else if(i==1) return 5*pow(P.x,4)-5*pow(P.y,4);
    // else return 60*P.x*P.x*P.y-20*P.y*P.y*P.y-5;
  }
  R fun_pressure(const R2 P, int i) {
    return 60*P.x*P.x*P.y-20*P.y*P.y*P.y-5;
  }
}

using namespace Stokes_Data;

// NEW VERSION

int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  int nx = 10;
  int ny = 10;
  int d = 2;

  vector<double> ul2,pl2,h, convu, convp;
  TaylorHood2 FEstokes;

  std::string pathOutpuFolder = "../../outputFiles/FEMstokes2/test/";
  std::string pathOutpuFigure = "../../outputFiles/FEMstokes2/test/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"data.dat", std::ofstream::out);

  for(int i=0;i<3;++i) {
    std::cout << "\n ------------------------------------- " << std::endl;
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);
    FESpace2 Vh(Th, FEstokes);

    Lagrange2 FE_velocity(2);
    FESpace2 Uh(Th, FE_velocity);
    FESpace2 Ph(Th, DataFE<Mesh2>::P1);
    FEM<Mesh2> stokes({&Uh, &Ph});
    long idx0_ph = Uh.NbDoF();

    // Stokes CutFEM problem
    Fun_h fh(Uh, fun_rhs);
    Fun_h gh(Uh, fun_boundary);

    Normal n;
    // FunTest u(Vh,d), p(Vh,1,d), v(Vh,d), q(Vh,1,d);
    FunTest u(Uh,d), p(Ph,1), v(Uh,d), q(Ph,1);
    R mu = 2;
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);

    stokes.addBilinear(
      contractProduct(2*mu*Eps(u),Eps(v))
      - innerProduct(p, div(v))
      + innerProduct(div(u), q)
    );
    // a(u,v)_dOmega
    stokes.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
      + innerProduct(p, v.t()*n)
      - innerProduct(u.t()*n, q)
      - innerProduct(2.*mu*Eps(u)*n, v)
      - innerProduct(u, 2.*mu*Eps(v)*n)
      , {1,2,3,4}
    ) ;

    // l(v)_Omega
    stokes.addLinear(
      innerProduct(fh.expression(2),u)
    );

    // l(v)_dOmega
    stokes.addLinearFormBorder(
        innerProduct(gh.expression(2), lambdaB*v)
      - innerProduct(gh.expression(2), 2.*mu*Eps(v)*n)
      - innerProduct(gh.expression(2), p*n)
      , {1,2,3,4}
    );

    stokes.addLagrangeMultiplier(
      innerProduct(1.,q), 0.
    );
    matlab::Export(stokes.mat, "matStokes.dat");
    return 0;
    stokes.solve();

    Fun_h uh(Uh, stokes.rhs);
    Rn_ data_ph(stokes.rhs(SubArray(Ph.NbDoF(),idx0_ph)));
    Fun_h ph(Ph, data_ph);
    R errU = L2norm(uh, fun_boundary,0,2);
    R errP = L2norm(ph, fun_pressure,0,1);

    // Fun_h sol(Uh, stokes.rhs);
    // Paraview2 writerS(Vh, pathOutpuFigure+"stokesTest"+to_string(i)+".vtk");
    // writerS.add(uh,"velocity", 0, 2);
    // writerS.add(ph,"pressure", 0, 1);

    ul2.push_back(errU);
    pl2.push_back(errP);
    h.push_back(1./nx);
    if(i==0) {convu.push_back(0); convp.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
  }
  std::cout << "\n\n h \t err u \t\t convU \t\t err p \t\t convP" << std::endl;
  for(int i=0;i<ul2.size();++i) {
    std::cout << h[i] << "\t" << ul2[i] << "\t" << convu[i] << "\t\t"
    << pl2[i] << "\t" << convp[i] << std::endl;
  }

  std::cout << " Time computation \t" << CPUtime() - cpubegin << std::endl;

}
