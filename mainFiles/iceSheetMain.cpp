#include <iostream>
#include "cutFEMConfig.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
// #include "iceSheet/markerIceSheet.hpp"
// #include "iceSheet/iceSheet.hpp"
#include "../iceSheet/dataSara.hpp"

#include "finiteElement.hpp"
#include "baseCutProblem.hpp"
#include "../num/gnuplot.hpp"

namespace IceSheet_Data {
  R fun_rhs(const R2 P, int i) {
    const double yearinsec = 365.25*24*60*60;
    const R gravity = 9.8*yearinsec*yearinsec;
    return (i==1)? -gravity : 0;
  }
  R Pwater(const R2 P, int i) {
    double yearinsec = 365.25*24*60*60;
    const R rhow = 1000.0/(1.0e6*yearinsec*yearinsec);
    const R gravity = 9.8*yearinsec*yearinsec;
    return (P.y<0)? -rhow*gravity*P.y : 0;
  }

}

using namespace IceSheet_Data;

int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef TestFunction<2> FunTest;


  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  // import the background mesh
  Mesh Th("../../mesh/iceSheet.mesh");
  Th.info();

  FESpace Vh(Th);
  Paraview2 writerS(Vh, "iceSheetMesh.vtk");
  // DataSara data("../../iceSheet/");

  Marker surface(Th, "../../iceSheet/surf.dat");
  Marker bedRock(Th, "../../iceSheet/bed.dat");

  surface.make_interface();

  gnuplot::save(Th, "Th.dat");
  gnuplot::saveMarker(surface, "surface.dat");
  // gnuplot::saveMarker(bedRock, "bed.dat");

  // MarkerIceSheet marker(Th);
  //
  //
  // // gnuplot::save(marker);
  //
  //
  // // Space for Stokes problem
  // TaylorHood2 FEstokes;
  // FESpace Vh(Th, FEstokes);
  // Vh.info();
  //
  //
  // FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  // Rn lll(data.ls_sign.size());
  // lll -= data.ls_sign;
  // Fun2 levelSet(Lh, lll);//data.ls_sign);
  //
  // SubDomain2 Qh1(Lh, levelSet,  1);
  // CutFESpace2 Qh(Qh1, marker);
  //
  // SubDomain2 Vh1(Vh, levelSet,  1);
  // CutFESpace2 Wh(Vh1, marker);
  // Wh.info();
  //
  //
  // CutFEM<Mesh, Marker> iceDynamic(Wh, marker);
  //
  // enum labelBoundary { surface = 1, floating = 2, bedRock  = 3, symmetry = 4, cliff = 5};
  // const double yearinsec = 365.25*24*60*60;
  // const R rhoi = 900.0 /(1.0e6*yearinsec*yearinsec);
  // const R rhow = 1000.0/(1.0e6*yearinsec*yearinsec);
  // const R A = (4.6416e-25)*yearinsec*1.0e18;
  // const R mu_ = 1.0/pow(2.0*A,1.0/3) ;
  // const R C = 7.624e6/(1.0e6*yearinsec) ;
  // const R Beta = 7.624e6/(1.0e6*yearinsec) ;
  // const R gravity = 9.8*yearinsec*yearinsec;
  // const R h = Th[0].lenEdge(0);
  // const R sigma = 1;
  // const int d = 2;
  // // CutFEM_Parameter mu("mu",1.,10.);
  // // CutFEM_Parameter invmu("invmu",1.,1./10);
  // // const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  // // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  // const CutFEM_Parameter& invh(Parameter::invh);
  //
  //
  // Rn fhv, ghv, pwaterv;
  // interpolate(Wh, fhv, fun_rhs);
  // Fun2 fh(Wh, fhv);
  // interpolate(Qh, pwaterv, Pwater);
  // Fun2 pwater(Qh, pwaterv);
  //
  // Normal n;
  // Tangent t;
  // FunTest u(d,0,0), p(1,d,0), v(d,0,0), q(1,d,0);
  //
  // const double tt0 = CPUtime();
  //
  //
  // // a(u,v)_Omega
  // iceDynamic.addBilinearFormOmega(
  //   contractProduct(2*mu_*Eps(u),Eps(v))
  //   - innerProduct(p, div(v)) + innerProduct(div(u), q)
  // );
  // // l(v)_Omega
  // iceDynamic.addLinearFormOmega(
  //   innerProduct(fh,rhoi*v)
  // );
  //
  // FunTest Eun = (Eps(u)*n);
  // FunTest Snn = 2*mu_*Eun.t()*n - p;
  // FunTest Snt = 2*mu_*Eun.t()*t - p;
  //
  // // a(u,v)_gamma
  // iceDynamic.addBilinearFormInterface(
  //     innerProduct( -1.*Snn, v.t()*n)
  //   - innerProduct( u.t()*n, Snn)
  //   - innerProduct( Snt, v.t()*t)
  //   + innerProduct(Beta*u, v)
  //   + innerProduct(1e3*invh*u.t()*n, v.t()*n)
  //     ,
  //     symmetry
  // );
  // iceDynamic.addBilinearFormInterface(
  //   innerProduct( -1.*Snn, v.t()*n)
  //   - innerProduct( u.t()*n, Snn)
  //   + innerProduct(Beta*u, v)
  //   + innerProduct(1e3*invh*u.t()*n, v.t()*n)
  //   ,
  //   bedRock
  // );
  //
  // int li[2] = {floating, cliff};
  // KN<int> ll(2, li);
  // iceDynamic.addLinearFormInterface(
  //   innerProduct(pwater,v.t()*n) , ll
  // );
  //
  // iceDynamic.addLagrangeMultiplier(
  //   innerProduct(1.,p), 0.
  // );
  //
  //
  // const double tt1 = CPUtime();
  // std::cout << " Time assembly \t " << tt1 - tt0 << std::endl;
  // iceDynamic.solve();
  //
  // IceSheet iceSheet(Wh);
  // iceSheet.init(Wh);
  // iceSheet.solve(marker);
  //
  // // Fun2 fsol(Wh, iceDynamic.rhs);
  // Fun2 fsol(Wh, iceSheet.rhs);
  // VTKcutWriter2 writer1(fsol, marker,"iceSheet.vtk");
  // writer1.add("velocity", 0, 2);
  // writer1.add("pressure", 2, 1);

  std::cout << "Time computation \t" <<  CPUtime() - cpubegin << std::endl;
}

// need to check normal for interface and marker (should be <=>)
