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

  R fun_init(const R2 P, const int i){
    return 16*P.x*(1-P.x)*P.y*(1-P.y);
  }

}

using namespace TimeStokes_Data;


void convectionDiffusionReaction(FEM<Mesh2>& cdr ) {


}

int main(int argc, char** argv )
{

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;
  typedef ListItemVF BilinearForm;

  MPIcf cfMPI(argc,argv);

  const double cpubegin = CPUtime();

  int d = 2;              // dimension of the problem
  int nx = 20;            // node in x direction
  int ny = 20;            // node in y direction
  double dT = 0.01;
  GTime::time_step = dT;
  GTime::total_number_iteration = int(0.25/dT);

  // TaylorHood2 FEstokes;
  const GTypeOfFE<Mesh2>& FEspace(DataFE<Mesh2>::P1);
  const GTypeOfFE<Mesh1>& FEtime(DataFE<Mesh1>::P1Poly);

  const double E = 1;
  const R2 Beta(0,0);
  const double alpha = 0;
  FunTest u(1), v(1);
  // BilinearForm

  Mesh2 Th(nx, ny, 0., 0., 1., 1.);
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());

  FESpace2 Vh(Th, FEspace);
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);

  FEM<Mesh2> cdr(Vh, Ih);

  // // Stokes CutFEM problem
  // Rn fhv, ghv;
  // interpolate(Vh, fhv, fun_rhs);
  // Fun2 fh(Vh, fhv);
  // interpolate(Vh, ghv, fun_boundary);
  // Fun2 gh(Vh, ghv);
  Rn datau0;
  interpolate(Vh, datau0, fun_init);
  Fun2 u0(Vh, datau0);

  // Normal n;
  // FunTest u(1), v(1);

  ListItemVF ah_K = innerProduct(dt(u), v)
  + innerProduct(E*grad(u), grad(v));
  // + innerProduct(Beta*grad(u), v)
  // + innerProduct(alpha* u, v);
  std::cout << ah_K << std::endl;
  int iter = 0;
  while( iter < GTime::total_number_iteration ) {
    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);
  //
  //   // convectionDiffusionReaction()
  //
    cdr.addBilinearFormOmega(
      ah_K, In
    );
    cdr.addBilinearFormOmega(
      innerProduct(u,v)
    );
    cdr.addLinearFormOmega(
      innerProduct(u0,v)
    );

    cdr.solve();


    KN_<double> solv0(cdr.rhs(SubArray(Vh.NbDoF(), 0)));
    KN_<double> solv1(cdr.rhs(SubArray(Vh.NbDoF(), Vh.NbDoF())));

    Rn solv(Vh.NbDoF(), 0.);
    solv += solv0;
    solv += solv1;
    Fun2 sol(Vh, solv);
    VTKwriter2 writer(sol, "solTime"+to_string(iter)+".vtk");
    // writer.add("heat0", 0, 1);
    writer.add("heat", 0, 1);

    datau0 = solv;
    cdr.rhs=0.;

    iter += 1;
  }




  const double cpuend = CPUtime();

  std::cout << "Running time \t" << cpuend - cpubegin << std::endl;

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


}
