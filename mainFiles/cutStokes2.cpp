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
#include "../num/gnuplot.hpp"
#include "../util/redirectOutput.hpp"


#include "DA.hpp"


namespace CutStokes_Data {
R2 shift(0,0);

  R fun_levelSet(const R2 P, int i) {
    return sqrt((P.x-shift.x)*(P.x-shift.x) + (P.y-shift.y)*(P.y-shift.y)) - 2./3;
  }

  static R falpha1(const R r) {
    const R MU1 = 1e-1;
    const R MU2 = 1e-3;
    const R r0 = 2./3;
    return 1./MU1 + (1./MU2-1./MU1)*exp(r*r-r0*r0);
  }
  static R falpha2(const R r) { R MU2 = 1e-3;
     return 1./MU2;}

  static R falpha(const R r) {
    const R r0 = 2./3;
    return (r < r0)? falpha2(r) : falpha1(r);
  }

  static R fun_rhs(const R2 P, int i) {
    const R r2 = Norme2_2(P-shift);
    const R s = exp(-r2);
    R2 R(4*s*(r2-2)*P.y + 3*P.x*P.x, -4*s*(r2-2)*P.x );
    return (i<2)?R[i] : 0;
  }
  static R fun_boundary(const R2 P,const int i) {
    const R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha(r)*exp(-r*r)*R;
    return (i<2)?R[i] : 0;
  }

  static R2 fun_velocity1(const R2 P) {
    R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha1(r)*exp(-r*r)*R; return R;
  }
  static R2 fun_velocity2(const R2 P) {
    R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha2(r)*exp(-r*r)*R; return R;
  }
  static R fun_pressure1(const R2 P) { return pow(P.x,3) ;}
  static R fun_pressure2(const R2 P) {
    R sigma = 24.5;//700;
    return pow(P.x,3) + sigma*3./2.;
  }

  static R fun_solution(const R2 P, int ci, int domain) {
    if(domain == 0) return (ci < 2)? fun_velocity1(P)[ci] : fun_pressure1(P);
    else return (ci < 2)? fun_velocity2(P)[ci] : fun_pressure2(P);
  }
  static R fun_velocity(const R2 P, int ci, int domain) {
    if(domain == 0) return fun_velocity1(P)[ci];
    else            return fun_velocity2(P)[ci];
  }
  static R fun_pressure(const R2 P, int ci, int domain) {
    if(domain == 0) return fun_pressure1(P);
    else            return fun_pressure2(P);
  }

  R2 fparam(double t){ return R2(2./3*cos(t+1./3), 2./3*sin(t+1./3));}

}

using namespace CutStokes_Data;
//#define _USE_MARKER
#ifdef _USE_MARKER
int main(int argc, char** argv )
{
  const int d = 2;

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef LevelSet2 LevelSet;
  typedef GCurvature<Mesh2> Curvature;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;
  int nmarker = 20;


  vector<double> ul2,uh1,pl2,ph1,hv,convul2, convuh1,convpl2, convph1;


  int order_ls=2, order_curv=2;
  Lagrange2 FEcurvature(order_curv);

  std::string pathOutpuFolder = "../../outputFiles/cutStokes2/DynamicBubble/";
  std::string pathOutpuFigure = "../../outputFiles/cutStokes2/DynamicBubble/paraview/";
  assert(IsPathExist(pathOutpuFigure));
  std::ofstream outputData(pathOutpuFolder+"dataP"+to_string(order_curv)+"P"+to_string(order_ls)+".dat", std::ofstream::out);


  outputData << " phi(x) = sqrt(x*x + y*y/0.25) - 1.0" << std::endl;
  outputData << " Curvature " << order_curv << std::endl;
  outputData << " Mapping & LevelSet "<< order_ls << std::endl;
  outputData << " Curvature constant \t 1e-2" << std::endl;
  outputData << " Curvature h^(2k-2)" << std::endl;


  CutFEM_Parameter mu("mu",1e-1,1e-3);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",1e1,1e3);
  const R sigma = 24.5;//700;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);

  outputData << "\n Parameters of the problem : \n"
            << " mu_1 = " << mu.val1 << " and mu_2 = "<< mu.val2 << " \n"
            << " the surface tension sigma = " << sigma
            << std::endl;

  for(int i=0;i<3;++i) {

    double meshSize(2./nx);

    Mesh Th(nx, ny, -1., -1., 2., 2.);

    // create marker
    // Marker interface(Th, fparam, 0, 2*M_PI, nmarker);
    // interface.make_interface();

    // gnuplot::save(marker);
    // gnuplot::save(Th);
    // gnuplot::saveNode(marker);
    // // levelSet stuff
    // FESpace Lh_k(Th, DataFE<Mesh2>::P2);
    // Fun_h levelSetPk(Lh_k, fun_levelSet);
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    // projection(levelSetPk, levelSet);

    Interface2 interface(Th, levelSet.v);


    // ------------  MEAN CURVATURE  -------------
    // Build cut Space
    Mesh2 cutTh(interface);
    FESpace2 cutVh(cutTh, FEcurvature);

    // ------------  MEAN CURVATURE  -------------
    Curvature curvature(cutVh, interface);
    Fun_h meanCurvature(cutVh, curvature.rhs);
    {
      // FESpace Lh(Th, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      Paraview2 writer(cutVh, pathOutpuFigure+"curvature.vtk");
      writer.add(meanCurvature, "meanCurvature", 0, 2);
      // writer.add(ls, "levelSet", 0, 1);
    }


    // ------------  CUT STOKES PROBLEM  ---------------

    Lagrange2 FEvelocity(2);
    FESpace Uh(Th, FEvelocity);
    FESpace Ph(Th, DataFE<Mesh>::P1);

    CutFESpace2 WUh(Uh, interface, {1,-1});
    CutFESpace2 WPh(Ph, interface, {1,-1});

    // gnuplot::save(interface, "marker.dat");
    // gnuplot::save(cutTh, "Vh.dat");
    // gnuplot::save(WPh, 0, "vh1.dat");
    // gnuplot::save(WPh, 1, "vh2.dat");
    //
    // return 0;



    long idx0_ph = WUh.NbDoF();

    Fun_h fh(WUh, fun_rhs);
    Fun_h gh(WUh, fun_boundary);

    Normal n;
    FunTest u(WUh,d), p(WPh,1), v(WUh,d), q(WPh,1), p1(WPh,1,0,0);
    FunTest Eun = (Eps(u)*n);


    // Stokes problem
    // ----------------------------------------------
    CutFEM<Mesh2> stokes({&WUh, &WPh});



    //a(u,v)_Omega
    stokes.addBilinear(
      contractProduct(2*mu*Eps(u),Eps(v))
      - innerProduct(p, div(v)) + innerProduct(div(u), q)
    );
    // a(u,v)_gamma
    stokes.addBilinear(
      innerProduct(jump(u), -2*mu*average1(Eps(v)*n))
      + innerProduct(-2*mu*average1(Eps(u)*n), jump(v))
      + innerProduct(lambdaG*jump(u), jump(v))
      + innerProduct(average1(p), jump(v.t()*n))
      - innerProduct(jump(u.t()*n), average1(q))
      , interface
    );

    // a(u,v)_dOmega
    stokes.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
      + innerProduct(p, v.t()*n)
      - innerProduct(u.t()*n, q)
      - innerProduct(2.*mu*Eps(u)*n, v)
      - innerProduct(u, 2.*mu*Eps(v)*n)
    ) ;

    // l(v)_Domain
    stokes.addLinear(
      innerProduct(fh.expression(2),u)
    );

    ExpressionFunFEM<Mesh2> Kx(meanCurvature,0,op_id);
    ExpressionFunFEM<Mesh2> Ky(meanCurvature,1,op_id);
    FunTest vx(WUh,1,0), vy(WUh,1,1);
    stokes.addLinear(
        innerProduct(Kx, average2(vx))*sigma
      + innerProduct(Ky, average2(vy))*sigma
      , interface
    );

    // l(v)_dOmega
    stokes.addLinearFormBorder(
      innerProduct(gh.expression(2),lambdaB*v)
      - innerProduct(gh.expression(2),2.*mu*Eps(v)*n)
      - innerProduct(gh.expression(2), p*n)
    );

    stokes.addFaceStabilization(
      innerProduct(1e-2*h*jump(grad(u)*n), mu*jump(grad(v)*n))
    );
    R cch = pow(meshSize,3);
    stokes.addFaceStabilization(
      innerProduct(1e-2*cch*jump(grad(p).t()*n), invmu*jump(grad(q).t()*n))
    );


    stokes.addLagrangeMultiplier(
      innerProduct(1.,p1), 0.
      ,0
    );

    stokes.solve();


    Fun_h uh(WUh, stokes.rhs);
    Rn_ data_ph(stokes.rhs(SubArray(WPh.NbDoF(),idx0_ph)));
    Fun_h ph(WPh, data_ph);

    R errU = L2normCut(uh, fun_velocity,0,2);
    R errP = L2normCut(ph, fun_pressure,0,1);

    ul2.push_back(errU);
    pl2.push_back(errP);
    hv.push_back(meshSize);
    if(i==0) {convul2.push_back(0); convpl2.push_back(0);}
    else {
      convul2.push_back(log(ul2[i]/ul2[i-1])/log(hv[i]/hv[i-1]));
      convpl2.push_back(log(pl2[i]/pl2[i-1])/log(hv[i]/hv[i-1]));
    }

    // // Fun_h sol(Wh, stokes.rhs);
    // Paraview2 writer(WPh, levelSet, pathOutpuFigure+"dynamicDrop"+to_string(i)+".vtk");
    // writer.add(uh, "velocity", 0, 2);
    // writer.add(ph, "pressure", 0, 1);

    nx *= 2;
    ny *= 2;
    nmarker*=2;
    // shift.y += meshSize / 10;
  }

  std::cout << std::setprecision(2);
  std::cout << "\n\n h \t\t errL2 u \t convL2 U \t errL2 p \t convL2 P" << std::endl;
  for(int i=0;i<hv.size();++i) {
    std::cout << hv[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << pl2[i] << "\t\t" << convpl2[i]
    << std::endl;
  }
  outputData << std::setprecision(2);
  outputData << "\n\n h \t\t errL2 u \t convL2 U \t errL2 p \t convL2 P" << std::endl;
  for(int i=0;i<hv.size();++i) {
    outputData  << hv[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << pl2[i] << "\t\t" << convpl2[i]
    << std::endl;
  }

}

#else

int main(int argc, char** argv )
{
  const int d = 2;

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GLevelSet<Mesh> LevelSet;
  typedef GCurvature<Mesh2> Curvature;


  MPIcf cfMPI(argc,argv);
  // MPIcf::setLoopSplitWork("block");
  const double cpubegin = CPUtime();

  int nx = 50;
  int ny = 50;

  vector<double> ul2,uh1,pl2,ph1,hv,convul2, convuh1,convpl2, convph1;


  int order_ls=2, order_curv=2;
  Lagrange2 FEcurvature(order_curv);

  std::string pathOutpuFolder = "../../outputFiles/cutStokes2/DynamicBubble/";
  std::string pathOutpuFigure = "../../outputFiles/cutStokes2/DynamicBubble/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  CoutFileAndScreen myCout(pathOutpuFolder+"screen.txt");

  std::ofstream outputData(pathOutpuFolder+"dataP"+to_string(order_curv)+"P"+to_string(order_ls)+".dat", std::ofstream::out);


  myCout << "Hey from \t" << MPIcf::my_rank() << std::endl;
  outputData << " phi(x) = sqrt(x*x + y*y/0.25) - 1.0" << std::endl;
  outputData << " Curvature " << order_curv << std::endl;
  outputData << " Mapping & LevelSet "<< order_ls << std::endl;
  outputData << " Curvature constant \t 1e-2" << std::endl;
  outputData << " Curvature h^(2k-2)" << std::endl;


  CutFEM_Parameter mu("mu",1e-1,1e-3);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",1e1,1e3);
  const R sigma = 24.5;//700;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);

  outputData << "\n Parameters of the problem : \n"
            << " mu_1 = " << mu.val1 << " and mu_2 = "<< mu.val2 << " \n"
            << " the surface tension sigma = " << sigma
            << std::endl;

  for(int i=0;i<1;++i) {

    double meshSize(2./nx);

    Mesh Th(nx, ny, -1., -1., 2., 2.);

    // levelSet stuff
    FESpace Lh_k(Th, DataFE<Mesh2>::P2);
    Fun_h levelSetPk(Lh_k, fun_levelSet);
    FESpace Lh(Th, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    projection(levelSetPk, levelSet);

    Interface2 interface(Th, levelSet.v);
    // ------------  MEAN CURVATURE  -------------
    // Build cut Space
    Mesh2 cutTh(interface);
    FESpace2 cutVh(cutTh, FEcurvature);

    // Build mapping
    Mapping2 mapping(cutVh, levelSetPk);
    mapping.errorInfty(interface, fun_levelSet);

    // ------------  MEAN CURVATURE  -------------
    Curvature curvature(cutVh, interface, mapping);
    Fun_h meanCurvature(cutVh, curvature.rhs);

    // ------------  CUT STOKES PROBLEM  ---------------

    Lagrange2 FEvelocity(2);
    FESpace Uh(Th, FEvelocity);
    FESpace Ph(Th, DataFE<Mesh>::P1);


    CutFESpace2 WUh(Uh, interface, {1,-1});
    CutFESpace2 WPh(Ph, interface, {1,-1});

    long idx0_ph = WUh.NbDoF();

    Fun_h fh(WUh, fun_rhs);
    Fun_h gh(WUh, fun_boundary);

    Normal n;
    FunTest u(WUh,d), p(WPh,1), v(WUh,d), q(WPh,1), p1(WPh,1,0,0);
    FunTest Eun = (Eps(u)*n);


    // Stokes problem
    // ----------------------------------------------
    CutFEM<Mesh2> stokes({&WUh, &WPh});

    CutFEM_Parameter myPara("myPara",2,2);


    //a(u,v)_Omega
    stokes.addBilinear(
      contractProduct(2*mu*Eps(u),Eps(v))
      - innerProduct(p, div(v)) + innerProduct(div(u), q)
    );
    // a(u,v)_gamma
    stokes.addBilinear(
      innerProduct(jump(u), -2*mu*average1(Eps(v)*n))
      + innerProduct(-2*mu*average1(Eps(u)*n), jump(v))
      + innerProduct(lambdaG*jump(u), jump(v))
      + innerProduct(average1(p), jump(v.t()*n))
      - innerProduct(jump(u.t()*n), average1(q))
      , interface
    );

    // a(u,v)_dOmega
    stokes.addBilinearFormBorder(
      innerProduct(lambdaB*u,v)
      + innerProduct(p, v.t()*n)
      - innerProduct(u.t()*n, q)
      - innerProduct(2.*mu*Eps(u)*n, v)
      - innerProduct(u, 2.*mu*Eps(v)*n)
    ) ;

    // l(v)_Domain
    stokes.addLinear(
      innerProduct(fh.expression(2),u)
    );

    ExpressionFunFEM<Mesh2> Kx(meanCurvature,0,op_id);
    ExpressionFunFEM<Mesh2> Ky(meanCurvature,1,op_id);
    FunTest vx(WUh,1,0), vy(WUh,1,1);
    stokes.addLinear(
        innerProduct(Kx, average2(vx))*sigma
      + innerProduct(Ky, average2(vy))*sigma
      , interface
    );

    // l(v)_dOmega
    stokes.addLinearFormBorder(
      innerProduct(gh.expression(2),lambdaB*v)
      - innerProduct(gh.expression(2),2.*mu*Eps(v)*n)
      - innerProduct(gh.expression(2), p*n)
    );

    stokes.addFaceStabilization(
      innerProduct(1e-2*h*jump(grad(u)*n), mu*jump(grad(v)*n))
    );
    // R cch = pow(meshSize,3);
    stokes.addFaceStabilization(
      innerProduct(1e-2*pow(h,3)*jump(grad(p).t()*n), invmu*jump(grad(q).t()*n))
    );


    stokes.addLagrangeMultiplier(
      innerProduct(1.,p1), 0.
    );

    stokes.solve();



    Fun_h uh(WUh, stokes.rhs);
    Rn_ data_ph(stokes.rhs(SubArray(WPh.NbDoF(),idx0_ph)));
    Fun_h ph(WPh, data_ph);

    // std::cout << data_ph << std::endl;
    // getchar();

    R errU = L2normCut(uh, fun_velocity,0,2);
    R errP = L2normCut(ph, fun_pressure,0,1);

    ul2.push_back(errU);
    pl2.push_back(errP);
    hv.push_back(meshSize);
    if(i==0) {convul2.push_back(0); convpl2.push_back(0);}
    else {
      convul2.push_back(log(ul2[i]/ul2[i-1])/log(hv[i]/hv[i-1]));
      convpl2.push_back(log(pl2[i]/pl2[i-1])/log(hv[i]/hv[i-1]));
    }

    // Fun_h sol(Wh, stokes.rhs);
    // Paraview2 writer(WPh, levelSet, pathOutpuFigure+"dynamicDrop"+to_string(i)+".vtk");
    // writer.add(uh, "velocity", 0, 2);
    // writer.add(ph, "pressure", 0, 1);

    nx *= 2;
    ny *= 2;
    // shift.y += meshSize / 10;
  }

  myCout<< std::setprecision(2);
  myCout<< "\n\n h \t\t errL2 u \t convL2 U \t errL2 p \t convL2 P" << std::endl;
  for(int i=0;i<hv.size();++i) {
    myCout<< hv[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << pl2[i] << "\t\t" << convpl2[i]
    << std::endl;
  }
  outputData << std::setprecision(2);
  outputData << "\n\n h \t\t errL2 u \t convL2 U \t errL2 p \t convL2 P" << std::endl;
  for(int i=0;i<hv.size();++i) {
    outputData  << hv[i] << "\t\t" << ul2[i] << "\t\t" << convul2[i]
    << "\t\t" << pl2[i] << "\t\t" << convpl2[i]
    << std::endl;
  }

}

#endif
