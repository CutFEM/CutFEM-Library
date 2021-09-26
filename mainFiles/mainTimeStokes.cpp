#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "finiteElement.hpp"
#include "baseProblem.hpp"


namespace TimeStokes_Data {

}

using namespace TimeStokes_Data;


int main(int argc, char** argv )
{

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;


  const double cpubegin = CPUtime();

  int d = 2;              // dimension of the problem
  int nx = 10;            // node in x direction
  int ny = 10;            // node in y direction
  double dt = 2./nx/2;;
  GTime::time_step = dt;
  GTime::total_number_iteration = int(3./dt);//40;//int(0.25/dt);

  TaylorHood2 FEstokes;
  const GTypeOfFE<Mesh>& FETime(DataFE<Mesh1>::P1Poly);

  for(int i=0;i<1;++i) {

    // Mesh Th(nx, ny, -1., -1., 2., 2.);
    // Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
    //
    // FESpace2 Vh(Th, FEstokes);
    // FEM<Mesh2> stokes(Vh);

    // // Stokes CutFEM problem
    // Rn fhv, ghv;
    // interpolate(Vh, fhv, fun_rhs);
    // Fun2 fh(Vh, fhv);
    // interpolate(Vh, ghv, fun_boundary);
    // Fun2 gh(Vh, ghv);
    //
    // Normal n;
    // FunTest u(d), p(1,d), v(d), q(1,d);
    // R mu = 2;
    // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    //
    // stokes.addBilinearFormOmega(
    //   contractProduct(2*mu*Eps(u),Eps(v))
    //   - innerProduct(p, div(v)) + innerProduct(div(u), q)
    // );
    // // a(u,v)_dOmega
    // stokes.addBilinearFormBorder(
    //   innerProduct(lambdaB*u,v)
    //   + innerProduct(p, v.t()*n)
    //   - innerProduct(u.t()*n, q)
    //   - innerProduct(2.*mu*Eps(u)*n, v)
    //   - innerProduct(u, 2.*mu*Eps(v)*n)
    // ) ;
    // // l(v)_Omega
    // stokes.addLinearFormOmega(
    //   innerProduct(fh,u)
    // );
    // // l(v)_dOmega
    // stokes.addLinearFormBorder(
    //     innerProduct(gh,lambdaB*v)
    //   - innerProduct(gh,2.*mu*Eps(v)*n)
    //   - innerProduct(gh, p*n)
    // );
    //
    // stokes.addLagrangeMultiplier(
    //   innerProduct(1.,p), 0.
    // );
    //
    // stokes.solve();

    // FunCut2 sol(Wh, stokes.rhs);
    // VTKcutWriter2 writer(sol, levelSet, "solnew_"+to_string(i)+".vtk");
    // writer.add("velocity", 0, 2);
    // writer.add("pressure", 2, 1);

    nx *= 2;
    ny *= 2;
  }
}
