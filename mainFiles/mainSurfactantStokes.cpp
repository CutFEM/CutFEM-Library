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
#include "baseCutProblem.hpp"


namespace SurfactantStokes_Data {

}

using namespace SurfactantStokes_Data;


int main(int argc, char** argv )
{

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<2> FunTest;


  const double cpubegin = CPUtime();

  int d = 2;
  int nx = 10;
  int ny = 10;

  for(int i=0;i<1;++i) {

    Mesh Th(nx, ny, -1., -1., 2., 2.);

    // levelSet stuff
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Rn levelSet_arr;
    interpolate(Lh, levelSet_arr, fun_levelSet);

    // Init interface
    Interface2 interface(Th, levelSet_arr);

    // ------------  MEAN CURVATURE  -------------
    Mesh2 cutTh(interface);
    Lagrange2 FEcurv;
    FESpace2 cutVh(Th, FEcurv);
    // Curvature<Mesh> curvature(Vh, interface);
    // curvature.solve();
    Rn curvaturev;
    interpolate(cutVh, curvaturev, fun_meanCurvature);
    Fun2 meanCurvature(cutVh,curvaturev);

    // ------------  CUT STOKES PROBLEM  ---------------
    TaylorHood2 FEstokes;
    FESpace Vh(Th, FEstokes);

    // Init subDomain
    Fun2 levelSet(Lh, levelSet_arr);
    SubDomain2 Vh1(Vh, levelSet, 1);
    SubDomain2 Vh2(Vh, levelSet,-1);

    // Create the cutFEM Space
    KN<SubDomain2*> subDomainsVh(2);
    subDomainsVh(0) = &Vh1;
    subDomainsVh(1) = &Vh2;
    CutFESpace2 Wh(subDomainsVh, interface);

    // Stokes CutFEM problem
    Rn fhv, ghv;
    interpolate(Wh, fhv, fun_rhs);
    Fun2 fh(Wh, fhv);
    interpolate(Wh, ghv, fun_boundary);
    Fun2 gh(Wh, ghv);

    Normal n;
    FunTest u(d), p(1,d), v(d), q(1,d), p1(1,d,0);

    CutFEM<Mesh,Interface2> stokes(Wh, interface, levelSet);
    CutFEM_Parameter mu("mu",10.,1.);
    CutFEM_Parameter invmu("invmu",1./10,1.);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);

    const R h = Th[0].lenEdge(0);
    const R sigma = 1;

    //a(u,v)_Omega
    stokes.addBilinearFormOmega(
      contractProduct(2*mu*Eps(u),Eps(v))
      - innerProduct(p, div(v)) + innerProduct(div(u), q)
    );
    // a(u,v)_gamma
    stokes.addBilinearFormInterface(
      innerProduct(jump(u), -2*mu*average1(Eps(v)*n))
      + innerProduct(-2*mu*average1(Eps(u)*n), jump(v))
      + innerProduct(lambdaG*jump(u), jump(v))
      + innerProduct(average1(p), jump(v.t()*n))
      - innerProduct(jump(u.t()*n), average1(q))
    );
    // a(u,v)_dOmega
    stokes.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
      + innerProduct(p, v.t()*n)
      - innerProduct(u.t()*n, q)
      - innerProduct(2.*mu*Eps(u)*n, v)
      - innerProduct(u, 2.*mu*Eps(v)*n)
    ) ;
    // l(v)_Omega
    stokes.addLinearFormOmega(
      innerProduct(fh,u)
    );
    // l(v)_Gamma
    stokes.addLinearFormInterface(
      innerProduct(meanCurvature, average2(v))*sigma
    );
    // l(v)_dOmega
    stokes.addLinearFormBorder(
      innerProduct(gh,lambdaB*v)
      - innerProduct(gh,2.*mu*Eps(v)*n)
      - innerProduct(gh, p*n)
    );

    // stokes.addIntegralOnFace(
    //   innerProduct(1e-1*h*h*h*invmu*jump(grad(p).t()*n), jump(grad(q).t()*n))
    // + innerProduct(1e-2*h*mu*jump(grad(u)*n), jump(grad(v)*n))
    // );
    //
    stokes.addLagrangeMultiplier(
      innerProduct(1.,p1), 0.
    );
    stokes.solve();

    Rn uex;
    interpolate(Wh, uex, fun_solution);

    R errU = stokes.L2norm(uex,0,2);
    R errP = stokes.L2norm(uex,2,1);

    // TO DO :
    // change back curvature
    // and base 'cut' problem


    ul2.push_back(errU);
    pl2.push_back(errP);
    hv.push_back(1./nx);
    if(i==0) {convul2.push_back(0); convpl2.push_back(0);}
    else {
      convul2.push_back(log(ul2[i]/ul2[i-1])/log(hv[i]/hv[i-1]));
      convpl2.push_back(log(pl2[i]/pl2[i-1])/log(hv[i]/hv[i-1]));
    }

    // FunCut2 sol(Wh, stokes.rhs);
    // VTKcutWriter2 writer(sol, levelSet, "solnew_"+to_string(i)+".vtk");
    // writer.add("velocity", 0, 2);
    // writer.add("pressure", 2, 1);

    nx *= 2;
    ny *= 2;
  }

  std::cout << std::setprecision(2);
  std::cout << "\n\n h \t\t errL2 u \t convL2 U \t errL2 p \t convL2 P" << std::endl;
  for(int i=0;i<hv.size();++i) {
    std::cout << hv[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << pl2[i] << "\t\t" << convpl2[i]
    << std::endl;
  }
}
