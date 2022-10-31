#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "problem/paperII.hpp"
#include "problem/levelSet.hpp"
#include "problem/curvature_v2.hpp"
// #include "num/gnuplot.hpp"


namespace TestStatic2D {
  R fun_levelSet(const R2 P) { return sqrt(P.x*P.x + P.y*P.y) - 2./3;}
  R fun_init_surfactant(const R2 P) { return 1 + P.x*P.y;}
  R2 fun_init_velocity(const R2 P, const R t) { return R2(0,0);}
  R2 fun_rhs(const R2 P,const R t=0) { R2 R(0,0); return R;}
  R2 fun_boundary(const R2 P,const R t=0) { return R2(0.,0.);}
  R sigma(const R w) {return 1 - w;}

  // exact solution
  R2 fun_velocity(const R2 P, const R t, const int r) { return R2(0,0);}
  R fun_pressure(const R2 P , const R t, const int r) { return (r==1)*3./2;}
  R fun_surfactant(const R2 P, const R t) { return 1.;}
}

namespace TestDynamic2D {
  R fun_levelSet(const R2 P) { return sqrt(P.x*P.x + P.y*P.y) - 2./3;}
  R fun_init_surfactant(const R2 P, const R t) { return 1.0;}
  R2 fun_init_velocity(const R2 P, const R t) { return R2(0,0);}
  R2 fun_rhs(const R2 P,const R t=0) {
    const R r2 = Norme2_2(P), s = exp(-r2);
    R2 R(P.y*4*s*(r2-2) + 3*P.x*P.x,
    -P.x*4*s*(r2-2) );
    return R;
  }

  static R falpha(const R r, const int dom) {
    const R MU1 = 1;
    const R MU2 = 10;
    const R r0 = 2./3;
    return (dom==1)? 1./MU2 : 1./MU1 + (1./MU2-1./MU1)*exp(r*r-r0*r0);
  }
  R2 fun_boundary(const R2 P,const R t=0) {
    const R r = Norme2(P);
    R2 R(-P.y,P.x); R = falpha(r, 0)*exp(-r*r)*R;
    return R;
  }

  R sigma(const R w) {return 1.;}// - w;}

  // exact solution
  R2 fun_velocity(const R2 P, const R t, const int dom) {
    R r = Norme2(P);
    R2 R(-P.y,P.x); R = falpha(r, dom)*exp(-r*r)*R; return R;
  }
  R fun_pressure(const R2 P,const R t, const int dom) {
    return  pow(P.x,3) + (dom==1)*3./2;
  }
  R fun_surfactant(const R2 P, const R t) { return 1.;}
}


namespace TestBenchmark2D {
  R fun_levelSet(const R2 P) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}
  R fun_init_surfactant(const R2 P) { return 1.0;}
  R2 fun_init_velocity(const R2 P, const R t) { return R2(0,0);}
  R2 fun_rhs(const R2 P,const R t=0) { R2 R(0,-0.98); return R;}
  R2 fun_boundary(const R2 P,const R t=0) { return R2(0,0);}
  R sigma(const R w) {return 24.5 - w;}
  R Dsigma(const R w) {return - w;}

}

namespace TestShearDrop2D {
  R fun_levelSet(const R2 P) { return sqrt(P.x*P.x + P.y*P.y) - 1;}
  R fun_init_surfactant(const R2 P) { return 1.0;}
  R2 fun_init_velocity(const R2 P, const R t) { return R2(0.,0.);}
  //R Q = 0.2; return R2(Q*P.x,-Q*P.y);}
  R2 fun_rhs(const R2 P,const R t=0) { return R2(0.,0.);}
  R2 fun_boundary(const R2 P,const R t=0) { R Q = 0.2; return R2(Q*P.x,-Q*P.y);}
  R sigma(const R w) {return 1. - 0.5*w;}
  R Dsigma(const R w) {return - 0.5*w;}
}

namespace TestShearDrop2D_2 {
  R fun_levelSet(const R2 P) { return sqrt(P.x*P.x + P.y*P.y) - 1./2;}
  R fun_init_surfactant(const R2 P) { return 0.4;}
  R2 fun_init_velocity(const R2 P, const R t) { return R2(0.,0.);}
  //R Q = 0.2; return R2(Q*P.x,-Q*P.y);}
  R2 fun_rhs(const R2 P,const R t=0) { return R2(0.,0.);}
  R2 fun_boundary(const R2 P,const R t=0) { R ss = (P.y < 0)?-1 : 1;
    if (t <= 1) return R2(ss*(0.5*(1 - cos(M_PI*t))),0);
    else return R2(ss,0);
  }
  R sigma(const R w) {return 1. - 0.5*w;}
  R Dsigma(const R w) {return -w;}
}


int main(int argc, char** argv )
{

  // using namespace MainFunction;
  // using namespace TestDynamic2D;
  // using namespace TestStatic2D;
  // using namespace TestShearDrop2D ;
  // using namespace TestShearDrop2D_2;
  using namespace TestBenchmark2D;
  typedef Mesh2 Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Vertex Vertex;
  typedef typename Mesh::Rd Rd;
  typedef typename Mesh::RdHat RdHat;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;

  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  std::cout << " simplify to get kp and cutK " << std::endl;
  std::cout << " fixe array of interface, not pointer if possible " << std::endl;
  std::cout << " NEED TO FIX INTERPOLATION!! " << std::endl;

  // TEST::TestNStokes(false);
  // return 0;

  int nx = 20;
  int ny = 40;
  Mesh2 Th(nx, ny, 0., 0., 1., 2.);
  // Mesh2 Th(nx, nx, -1., -1., 2., 2.);
  //Mesh2 Th(nx, nx, -2., -2., 4., 4.);
  // Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  // gnuplot::exportMesh(Th);

  // PeriodicBC periodicBC(Th);
  // FESpace2 VVh(Th, DataFE<Mesh2>::P1, periodicBC.nbEqui, periodicBC.periodicBE);
  // for(int k=0;k<VVh.NbElement();++k) {
  //   std::cout << "------------------ " << k << " ---------------- " << std::endl;
  //   for(int i=0;i<VVh[k].NbDoF();++i)
  //     std::cout << (VVh[k](i)) << std::endl;
  // }
  // return 0;


  // initialize time stuff
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const double h = Th[0].lenEdge(0);
  const double dt = h/2;
  GTime::time_step = dt;
  GTime::total_number_iteration = int(3./dt);
  GTime::t0 = 0;
  GTime::current_iteration = 0;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << " mesh size :\t" << h << std::endl;
  std::cout << " dt :\t" << dt << std::endl;
  std::cout << " nx :\t" << nx << std::endl;
  std::cout << " number of time step :\t" << GTime::total_number_iteration << std::endl;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);



  // Space for the velocity field < get from stokes >
  // GTypeOfFERd<Mesh2> FEvel(&DataFE<Mesh2>::P2, 2); // Need same than Stokes
  KN<const GTypeOfFE<Mesh2>* > myFEVel(2);
  myFEVel(0) = &DataFE<Mesh2>::P2;
  myFEVel(1) = &DataFE<Mesh2>::P2;
  GTypeOfFESum<Mesh2> FEvel(myFEVel);
  FESpace2 VelVh(Th, FEvel);

  vector<Rn>  velocityField(nbTime);                 // for levelSet
  vector<Fun2> vel(nbTime);
  for(int i=0;i<nbTime;++i) velocityField[i].init(VelVh.nbDoF);
  for(int i=0;i<nbTime;++i) vel[i].init(VelVh , velocityField[i]);


  // levelSet stuff
  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  LevelSet2 levelSet(Lh, fun_levelSet);
  levelSet.dt = dt/(nbTime-1);
  // levelSet.solver = "MUMPS";

  // Fun2 solS(Lh, levelSet.rhs);
  // VTKwriter2 writerSS(solS, "levelSetNEw.vtk");
  // writerSS.add("levelSet11", 0, 1);


  vector<Rn> lsv(nbTime);
  vector<Fun2> ls(nbTime);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, lsv[i]);


  // Declaration of the interface
  KN<Interface2*> interface(nbTime);

  // FE for curvature problem
  // GTypeOfFERd<Mesh2> FEcurv(&DataFE<Mesh2>::P1, 2);
  KN<const GTypeOfFE<Mesh2>* > myFECurv(2);
  myFECurv(0) = &DataFE<Mesh2>::P1;
  myFECurv(1) = &DataFE<Mesh2>::P1;
  GTypeOfFESum<Mesh2> FEcurv(myFECurv);
  vector<Mesh2*> cutTh(nbTime);
  vector<FESpace2*> cutVh(nbTime);
  vector<Rn> curvatureData(nbTime);
  vector<Fun2> curvature(nbTime);


  // Space for Stokes problem
  KN<const GTypeOfFE<Mesh2>* > myFETH(3);
  myFETH(0) = &DataFE<Mesh2>::P2;
  myFETH(1) = &DataFE<Mesh2>::P2;
  myFETH(2) = &DataFE<Mesh2>::P1;
  GTypeOfFESum<Mesh2> TaylorHood(myFETH);
  FESpace2 Vh(Th, TaylorHood);//, &periodicBC);


  // initialization problem paper2
  PaperII2 problem(Ih);
  // problem.solver = "MUMPS";
  problem.qft = &qTime;
  problem.sigma = sigma;
  problem.Dsigma = Dsigma;
  problem.eps = 0.1;
  problem.mu1 = 10;
  problem.mu2 = 1;
  problem.rho1 = 1000;
  problem.rho2 = 100;
  problem.fun_boundary = fun_boundary;
  problem.fun_source = fun_rhs;
  problem.fun_init_surfactant = fun_init_surfactant;

  int iter = 0;
  while( iter < GTime::total_number_iteration ) {
    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[GTime::current_iteration]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;

    int niteration = 1;//(iter==0) + 1;

    for(int ip=0;ip<niteration;++ip) {  //iterate to get good velocity field on tn

      // Init subDomain
      SubDomain2 Vh1(Vh,  1);
      SubDomain2 Vh2(Vh, -1);


      if(ip > 0) {
        levelSet.Up = lsv[0];
        levelSet.rhs = lsv[0];
      }


      // computation of the nbTime interfaces
      for(int i=0;i<nbTime;++i) {
        interface[i] = new Interface2(Th);
        interface[i]->make_patch(levelSet.rhs);

        // computation of the curvature
        cutTh[i] = new Mesh2(*interface[i]);
        cutVh[i] = new FESpace2(*cutTh[i], FEcurv);
        Curvature<Mesh2> Pcurvature(cutVh[i], interface[i]);
        Pcurvature.assembly();
        Pcurvature.solve();


        curvatureData[i].init(Pcurvature.rhs);         // save the data;
        curvature[i].init(cutVh[i], curvatureData[i]);

        lsv[i].init(levelSet.rhs);
        ls[i].init(Lh, lsv[i]);

        Vh1.add(Vh, ls[i]);
        Vh2.add(Vh, ls[i]);

        if(i<nbTime-1) levelSet.solve(vel[i], vel[i+1]);
      }

      Vh1.finalize(Vh);
      Vh2.finalize(Vh);

      KN<SubDomain2*> subDomainsVh(2);
      subDomainsVh(0) = &Vh1;
      subDomainsVh(1) = &Vh2;
      CutFESpace2 Wh(subDomainsVh, interface);

      //if(iter == 0)
      Wh.info();

      // Create the Active Mesh
      Mesh2 cutThTime(interface);             // surface mesh
      FESpace2 cutWh(cutThTime, DataFE<Mesh2>::P1);   // FE for surfactant

      problem.init(Wh, cutWh);
      problem.solve(ls, curvature, 0);

      std::cout << " DEFORMATION \t : \t D = " << problem.compute_deformation() << std::endl;

      // updating velocity
      for(int i=0;i<nbTime;++i ) {
        R tt = GTime::current_time() + i*dt/(nbTime-1);
        problem.set_velocity(velocityField[i], VelVh, ls[i], tt);
      }

      if(MPIcf::IamMaster() && ip == niteration-1) {
        KN_<double> solv0(problem.F(SubArray(Wh.NbDoF(), 0)));
        KN_<double> solv1(problem.F(SubArray(Wh.NbDoF(), Wh.NbDoF())));
        Rn solv(Wh.NbDoF(), 0.);
        solv += solv0;
        Fun2 sol(Wh, solv);
        // VTKcutWriter2 writer0(sol, ls[0], "solShear0_"+to_string(iter)+".vtk");
        // writer0.add("velocity", 0, 2);
        // writer0.add("pressure", 2, 1);
        VTKcutWriter2 writer1(sol, ls[nbTime-1], "solShear_"+to_string(iter)+".vtk");
        // solv += solv1;
        writer1.add("velocity", 0, 2);
        writer1.add("pressure", 2, 1);

        int idxs0 = Wh.NbDoF()*In.NbDoF()+2;
        KN_<double> solvS0(problem.F(SubArray(cutWh.NbDoF(), idxs0)));
        KN_<double> solvS1(problem.F(SubArray(cutWh.NbDoF(), idxs0+cutWh.NbDoF())));
        Rn solvS(cutWh.NbDoF(), 0.);
        solvS += solvS0;
        Fun2 solS(cutWh, solvS);
        VTKwriter2 writerS(solS, "surfactantShear_"+to_string(iter)+".vtk");
        // writerS.add("surfactant0", 0, 1);
        solvS += solvS0;
        writerS.add("surfactant1", 0, 1);

        Rn lls;
        // interpolate(cutWh,lls,fun_levelSet);
        // restriction(Lh, lsv[0], cutWh, lls );
        // solvS = lls;
        // writerS.add("levelSet0", 0, 1);
        restriction(Lh, lsv[nbTime-1], cutWh, lls );
        solvS = lls;
        writerS.add("levelSet1", 0, 1);



        VTKwriter2 writerV(ls[nbTime-1],
          "shearLevelSet_"+to_string(iter)+".vtk");
          writerV.add("levelSet", 0, 1);


        }


      if(ip == niteration-1) problem.setMap();
      if(ip == niteration-1) problem.export_Surfactant(iter);

      std::cout << "\n \n" << std::endl;


    }


    // Now the problem is solved so we save the solution for the next time step
    // problem.setMap();
    // getchar();

    // R errU = 0.;
    // R errP = 0.;
    // R errS = 0.;
    // problem.L2norm(ls[0],
    // 		  errU, fun_velocity,
    // 		  errP, fun_pressure,
    // 		  errS, fun_surfactant);

    // std::cout << "error velocity   \t" << errU << std::endl;
    // std::cout << "error pressure   \t" << errP << std::endl;
    // std::cout << "error surfactant \t" << errS << std::endl;






    // for(int i=0;i<nbTime;++i) {
    //   VTKwriter2 writerV(vel[i],
    // 			 "benchVel"+to_string(iter)+"_"+to_string(i)+".vtk");
    //   writerV.add("velocity", 0, 2);
    // }
    // getchar();
    iter++;

    // delete temporary data
    for(int i=0;i<nbTime;++i) {
      delete interface[i];
      delete cutTh[i];
      //      delete cutVh[i];
    }
  }


  const double cpuend = CPUtime();
  std::cout << " CPU Time \t " << cpuend - cpubegin << std::endl;
  return 0;
}
