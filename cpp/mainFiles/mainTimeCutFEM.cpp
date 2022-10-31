#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "baseCutProblem.hpp"
#include "problem/curvature_v2.hpp"
#include "problem/reinitialization.hpp"
#include "time_stuff.hpp"
#include "num/gnuplot.hpp"    // std::iter_swap
#include "projection.hpp"

#define RBUBBLE2D
// #define SURFACTANT

#ifdef SURFACTANT

#endif



#ifdef RBUBBLE2D
namespace TestBenchmark2D {
  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}
  // R fun_levelSet(const R2 P, const int i) { return ((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25*0.25;}

  R fun_rhs(const R2 P,const int i) { R2 R(0,-0.98); return (i<2)?R[i] : 0;}
  R fun_boundary(const R2 P,const int i) { R2 R(0,0);  return (i<2)?R[i] : 0;}
  // R fun_velocity(const R2 P,const int i) { return (i==0)? 0 : P.x*(1-P.x);}
  R fun_test2(const R2 P,const int i) { R2 R(6,-6); return (i<2)?R[i] : 0;}
  R fun_init_surfactant(const R2 P, int i) { return 1.;}
  R sigma(const R w) {return 24.5 - w;}
  R Dsigma(const R w) {return - w;}

}

// using namespace CutStokes_Data;
using namespace TestBenchmark2D;


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

  int nx = 40;
  int ny = 80;

  Mesh2 Th(nx, ny, 0., 0., 1., 2.);
  double meshSize(1./40);
  double dT = meshSize/2;
  GTime::time_step = dT;
  GTime::total_number_iteration = 100;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());

  std::cout << " Mesh size \t" << meshSize << std::endl;
  std::cout << " Time Step \t" << dT << std::endl;
  std::cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  std::cout << " ----------------------------------------------" << std::endl;


  TaylorHood2 FEstokes;
  FESpace2 Vh(Th, FEstokes);//, &periodicBC);
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  // const QuadratureFormular1d& qTime(QF_Euler);
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();

  // Space for the velocity field
  KN<const GTypeOfFE<Mesh2>* > myFEVel(2);
  myFEVel(0) = &DataFE<Mesh2>::P2;
  myFEVel(1) = &DataFE<Mesh2>::P2;
  GTypeOfFESum<Mesh2> FEvel(myFEVel);
  FESpace2 VelVh(Th, FEvel);

  vector<Fun_h> vel(nbTime);
  for(int i=0;i<nbTime;++i) vel[i].init(VelVh);//, fun_velocity);


  // levelSet stuff
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  FESpace2 Lh_k(Th, DataFE<Mesh2>::P2);
  // double dt_levelSet = dT;///(nbTime-1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  for(int i=0;i<nbTime;++i) ls_k[i].init(Lh_k, fun_levelSet);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);
  projection(ls_k[0], ls[nbTime-1]);



  // projection(ls_k[0], ls[0]);
  // Interface2 gamma(Th);
  // gamma.make_patch(ls[0].v);
  //
  // Fun2 solS(Lh, ls[0].data);
  // {
  //   VTKwriter2 writerS(solS, "withoutReinit.vtk");
  //   writerS.add("levelSet", 0, 1);
  // }
  // Reinitialization2(ls_k[0], gamma, dt_levelSet/4);
  // projection(ls_k[0], ls[0]);
  //
  // {
  //   VTKwriter2 writerS(solS, "Reinit.vtk");
  //   writerS.add("reinit", 0, 1);
  // }
  //
  // return 0;


  // Declaration of the interface
  KN<Interface2*> interface(nbTime);

  // FE for curvature problem
  KN<const GTypeOfFE<Mesh2>* > myFECurv(2);
  myFECurv(0) = &DataFE<Mesh2>::P2;
  myFECurv(1) = &DataFE<Mesh2>::P2;
  GTypeOfFESum<Mesh2> FEcurv(myFECurv);


  // Stokes problem
  CutFEM<Mesh2,Interface2> stokes(qTime);
  CutFEM_Parameter mu("mu",10.,1.);
  CutFEM_Parameter rho("rho",1000.,100.);
  CutFEM_Parameter invmu("invmu",1./10,1.);

  const R epsilon_surfactant = 0.1;
  const R sigma = 24.5;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);
  int dataDirichlet[2] = {1,3};
  int dataNeumann[2] = {2,4};
  KN<int> dirichlet(2,dataDirichlet), neumann(2,dataNeumann);


  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;

    // Init subDomain
    SubDomain2 Vh1(Vh,  1);
    SubDomain2 Vh2(Vh, -1);


    // if(ip > 0) {
    //   levelSet.Up = lsv[0];
    //   levelSet.rhs = lsv[0];
    // }

    // ls.begin()->swap(ls[nbTime-1]);
    ls_k.begin()->swap(ls_k[nbTime-1]);
    // computation of the nbTime interfaces
    for(int i=0;i<nbTime;++i) {
      projection(ls_k[i], ls[i]);
      interface[i] = new Interface2(Th);
      interface[i]->make_patch(ls[i].v);

      Vh1.add(Vh, ls[i]);
      Vh2.add(Vh, ls[i]);

      if(i<nbTime-1) {
        LevelSet levelSet(ls_k[i], vel[i], vel[i+1], dt_levelSet);
        // if(iter%5==0) Reinitialization2(ls_k[i], *interface[i], dt_levelSet/4);
        // ls_k[i+1].init(levelSet.rhs);
      }
    }

    Vh1.finalize(Vh);
    Vh2.finalize(Vh);

    KN<SubDomain2*> subDomainsVh(2);
    subDomainsVh(0) = &Vh1;
    subDomainsVh(1) = &Vh2;
    CutFESpace2 Wh(subDomainsVh, interface);

    // Create the Active Mesh
    Mesh2 cutThTime(interface);             // surface mesh
    FESpace2 cutWh(cutThTime, DataFE<Mesh2>::P1);   // FE for surfactant
    cutWh.backSpace = &Lh;  // svae backSpace to save solution

    /*
                    PROBLEM DEFINITION
    */
    Normal n;
    FunTest du(Wh,d), dp(Wh,1,d), v(Wh,d), q(Wh,1,d), dp1(Wh,1,d,0);
    FunTest Eun = (Eps(du)*n);
    FunTest ds(cutWh,1), r(cutWh,1);

    stokes.initSpace(Wh, In);
    stokes.add(cutWh, In);
    stokes.initInterface(interface);

    Rn datau0;
    stokes.initialSolution(datau0, cutWh);
    Rn uh(datau0);
    Fun_h u0(Wh, datau0);
    Fun_h fh(Vh, fun_rhs);

    int n0 = Wh.NbDoF()*In.NbDoF();
    KN_<double> datas0(datau0(SubArray(cutWh.NbDoF(),n0)));
    Fun_h u0s(cutWh, datas0);

    // initial value
    if(iter == 0) {
      interpolate(cutWh, datas0, fun_init_surfactant);
    }

    // loop for Newton
    int iterNewton  = 0;
    while(1) {

      std::cout << " Begin assembly " << std::endl;
      double tt0 = CPUtime();
      stokes.addBilinearFormDomain(
        innerProduct(dt(du), rho*v)
        + contractProduct(2*mu*Eps(du),Eps(v))
        - innerProduct(dp, div(v)) + innerProduct(div(du), q)
        , In
      );
      stokes.addBilinearFormInterface(
        innerProduct(jump(du), -2*mu*average1(Eps(v)*n))
        + innerProduct(-2*mu*average1(Eps(du)*n), jump(v))
        + innerProduct(lambdaG*jump(du), jump(v))
        + innerProduct(average1(dp), jump(v.t()*n))
        - innerProduct(jump(du.t()*n), average1(q))
        , In
      );

      for(int i=0;i<nbTime;++i) {      // computation of the curvature
        Mesh2 cutTh(*interface[i]);
        FESpace2 cutVh(cutTh, FEcurv);
        Mapping2 mapping(cutVh, ls_k[i]);

        Curvature curvature(cutVh, *interface[i], mapping);
        Fun_h meanCurvature(cutVh, curvature.rhs);
        ExpressionFunFEM<Mesh2> Kx(meanCurvature,0,op_id);
        ExpressionFunFEM<Mesh2> Ky(meanCurvature,1,op_id);
        FunTest vx(Wh,1,0), vy(Wh,1,1);

        stokes.addBilinearFormInterface(
           innerProduct(Kx*ds, average2(vx))
          +innerProduct(Ky*ds, average2(vy))
          +innerProduct(grad(ds), average2(v))
          ,i , In
        );
        stokes.addLinearFormInterface(
          -innerProduct(Kx, average2(vx))*sigma
          -innerProduct(Ky, average2(vy))*sigma
          ,i , In
        );
      }

      stokes.addBilinearFormInterface(
          innerProduct(dt(ds), r)
          + innerProduct(epsilon_surfactant*gradS(ds), gradS(r))
        , In
      );
      stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du,1e5*v)
        + innerProduct(dp, v.t()*n)
        - innerProduct(du.t()*n, q)
        - innerProduct(2.*mu*Eps(du)*n, v)
        - innerProduct(du, 2.*mu*Eps(v)*n)
        , In, dirichlet
      );
      stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du.t()*n,1e5*v.t()*n)
        + innerProduct(dp, v.t()*n)
        - innerProduct(du.t()*n, q)
        - innerProduct(2.*mu*Eun.t()*n, v.t()*n)
        - innerProduct(du.t()*n, 2.*mu*Eun.t()*n)
        , In, neumann
      );
      stokes.addFaceStabilization(
        innerProduct(1e-2*h*jump(grad(du)*n), mu*jump(grad(v)*n)),
        In
      );
      R cch = pow(meshSize,3);
      stokes.addFaceStabilization(
        innerProduct(1e-2*cch*jump(grad(dp).t()*n), invmu*jump(grad(q).t()*n)),
        In
      );
      stokes.addEdgeIntegral(
        innerProduct(1e-2*h*jump(grad(ds)*n), mu*jump(grad(r)*n)),
        In
      );
      // impose initial condition
      stokes.addBilinearFormDomain(
        innerProduct(rho*du,v)
      );
      stokes.addBilinearFormInterface(
        innerProduct(ds, r)
      );

      std::cout << " Time Full assembly " << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();
      stokes.addMatMul(uh);
      std::cout << " Time RHS assembly A*u0 " << CPUtime() - tt0 << std::endl;

      tt0 = CPUtime();
      FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      ExpressionFunFEM<Mesh2> u1(u0,0,op_id,0,0), u2(u0,1,op_id,0,0);
      ExpressionFunFEM<Mesh2> s(u0s,0,op_id);
      ExpressionFunFEM<Mesh2> dxu1(u0,0,op_dx,0,0), dxu1nxnx(u0,0,op_dx,0,0), dyu1nxny(u0,0,op_dy,0,0),
                              dxu2nxny(u0,1,op_dx,0,0), dyu2(u0,1,op_dy,0,0), dyu2nyny(u0,1,op_dy,0,0);

      dxu1nxnx.addNormal(0); dxu1nxnx.addNormal(0);
      dyu1nxny.addNormal(0); dyu1nxny.addNormal(1);
      dxu2nxny.addNormal(0); dxu2nxny.addNormal(1);
      dyu2nyny.addNormal(1); dyu2nyny.addNormal(1);

      // stokes.addBilinearFormDomain(
      //     innerProduct(du1*dx(u1) + du2*dy(u1), v1)
      //   + innerProduct(du1*dx(u2) + du2*dy(u2), v2)
      //   + innerProduct(u1*dx(du1) + u2*dy(du1), v1)
      //   + innerProduct(u1*dx(du2) + u2*dy(du2), v2)
      //   ,In
      // );
      // stokes.addLinearFormDomain(
      //     innerProduct(u1*dx(u1) + u2*dy(u1), v1)
      //   + innerProduct(u1*dx(u2) + u2*dy(u2), v2)
      //   , In
      // );

      FunTest dux1(Wh,1,0,0), duy1(Wh,1,1,0);
      // stokes.addBilinearFormInterface(
      //     innerProduct(dux1*dx(s) + duy1*dy(s), r)
      //   + innerProduct(u1*dx(ds)  + u2*dy(ds) , r)
      //   + innerProduct(ds*(dxu1 - dxu1nxnx - dyu1nxny), r)
      //   + innerProduct(ds*(dyu2 - dyu2nyny - dxu2nxny), r)
      //   + innerProduct(s*(dxS(dux1) + dyS(duy1)), r)
      //   , In
      // );
      // stokes.addLinearFormInterface(
      //       innerProduct(u1*dx(s) + u2*dy(s), r)
      //       + innerProduct(s*(dxu1 - dxu1nxnx - dyu1nxny), r)
      //       + innerProduct(s*(dyu2 - dyu2nyny - dxu2nxny), r)
      //   , In
      // );
      std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;


      tt0 = CPUtime();
      stokes.addLinearFormDomain(
        -innerProduct(fh,rho*v)
        , In
      );

      // impose initial condition
      stokes.addLinearFormDomain(
        -innerProduct(u0,rho*v)
      );
      stokes.addLinearFormInterface (
        -innerProduct(u0s, r)
      );


      std::cout << " Time assembly rhs " << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();

      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0., nbTime-1, In
      );
      stokes.solve();

      KN_<double> dw(stokes.rhs(SubArray(stokes.nDoF, 0)));
      uh -= dw;

      KN_<double> dws(dw(SubArray(cutWh.NbDoF()*In.NbDoF(),n0)));
      double dist = dws.linfty();
      std::cout << " dw \t" << dist << std::endl;
      std::cout << " Norm velocity \t" << uh.l2()/Wh.NbDoF() << std::endl;
      stokes.saveSolution(uh, cutWh);

      iterNewton++;
      stokes.rhs.resize(stokes.nDoF); stokes.rhs = 0.0;

      if(iterNewton == 7 || dist < 1e-5) break;

    }

    for(int i=0;i<nbTime;++i) {
      Fun_h sol(Wh,In,uh);
      set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    }


//  -----------------------------------------------------
//                     PLOTTING
//  -----------------------------------------------------
    int ndf = 0;
    int ndfs = n0;
    Rn solv(Wh.NbDoF(), 0.);
    Rn sols(cutWh.NbDoF(), 0.);

    if(iter%1 == 0) {
      for(int i=0;i<Ih[0].NbDoF();++i){

        KN_<double> solvi(uh(SubArray(Wh.NbDoF(), ndf)));
        ndf += Wh.NbDoF();
        solv += solvi;

        KN_<double> solsi(uh(SubArray(cutWh.NbDoF(), ndfs)));
        ndfs += cutWh.NbDoF();
        sols += solsi;
      }

      {
        Fun2 sol(Wh, solv);
        VTKcutWriter2 writerr(sol, ls[nbTime-1], "navierStokesSurfactantNewton40_"+to_string(iterfig)+".vtk");
        writerr.add("velocity", 0, 2);
        writerr.add("pressure", 2, 1);
      }
      // {
      //   Fun2 sol(Lh, ls[nbTime-1].data);
      //   VTKwriter2 writer(sol, "levelSetReinit_"+to_string(iterfig)+".vtk");
      //   writer.add("levelSet", 0, 1);
      // }
      {
        Fun2 sol(cutWh, sols);
        VTKwriter2 writer(sol, "surfactantNewton40_"+to_string(iterfig)+".vtk");
        writer.add("surfactant", 0, 1);
        writer.add(ls[nbTime-1], "levelSet.vtk", 0, 1);
      }

      iterfig++;
    }

    if(nbTime == 1) {
      LevelSet levelSet(ls_k[0], vel[0], vel[0], dt_levelSet);
      // if(iter%5 == 0) {
      // }
      ls_k[0].init(levelSet.rhs);
      // Reinitialization(ls_k[0], *interface[0], dt_levelSet);


    }

    iter += 1;
  }

}
#endif
