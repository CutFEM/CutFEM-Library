#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "baseCutProblem.hpp"
#include "curvature.hpp"
#include "reinitialization.hpp"
#include "../time_stuff.hpp"
#include "../num/gnuplot.hpp"    // std::iter_swap
#include "projection.hpp"
#include "../util/redirectOutput.hpp"

/*
Simulation of a rising bubble.
We consider the navier Stokes equations coupled with the insoluble surfactant
We are using the levelSet P2 with a reinitializtion step every 5 time iterations
We are using finite element in time, with the simpson rule in time.
*/


namespace TestBenchmark2D {
  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.1)*(P.x-0.1) + (P.y-0.0)*(P.y-0.0)) - 0.3;}
  R fun_rhs(const R2 P,const int i) { R2 R(0,-0.98); return (i<2)?R[i] : 0;}
  R fun_velocity(const R2 P, const int i){
    if(i == 0) return -0.5*(1+cos(M_PI*P.x))*sin(M_PI*P.y);
    else       return  0.5*(1+cos(M_PI*P.y))*sin(M_PI*P.x);
  }
  R fun_init_surface(const R2 P, int i) { return 0.;}
  R fun_w(double r) {
    return 0.5*(1-cos((r-0.3)*M_PI/(0.5*0.3)));
  }
  R fun_init_bulk(const R2 P, int i) {
    double r = sqrt((P.x-0.1)*(P.x-0.1) + P.y*P.y);
    if(r > 1.5*0.3) return 0.5*(1-P.x*P.x)*(1-P.x*P.x);
    else if(0.3 < r && r <= 1.5*0.3) return 0.5*(1-P.x*P.x)*(1-P.x*P.x)*fun_w(r);
    else return 0.;}

  R sigma(const R w) {return 24.5 - 0.5*w;}
  R Dsigma(const R w) {return - 0.5*w;}

  R fun_x(const R2 P,const int i) {return P[i];}
  R fun_1(const R2 P,const int i) {return 1;}
}


using namespace TestBenchmark2D;


int main(int argc, char** argv )
{

  const int d = 2;
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef LevelSet2 LevelSet;
  typedef GCurvature<Mesh> Curvature;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;

  Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  double meshSize(2./nx);
  // Mesh2 Th("../mesh/RBadaptMesh1_3060.msh");
  // int nx = 30;
  // int ny = 60;
  // double meshSize((1*2)/sqrt(nx*ny));

  int divisionMeshsize = 2;
  double Tend = 2;
  double dT = meshSize/divisionMeshsize;
  GTime::total_number_iteration = Tend/dT;
  dT = Tend / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  // Create files to save results
  std::string meshstr = to_string(nx)+to_string(ny);
  std::string pathOutpuFolder = "../../outputFiles/bulkSurfactantNavierStokes/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/bulkSurfactantNavierStokes/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;


  myCout << " SIMULATION OF A BULK SURFACE SURFACTANT COUPLED PROBLEM. \n"
            << " We consider the navier Stokes equations coupled to the insoluble surfactant problem.\n"
            << " We solve the non linear system using the Newton method \n"
            << " We approximate the surfactant concentration using P1/P1 space time finite elements  \n"
            << " We are using the levelSet P2 with a reinitializtion step every 5 time iterations \n"
            << " We are using P1 finite element in time, with the Simpson rule as quadrature rule "
            << std::endl << std::endl;



  myCout << " Mesh size \t" << meshSize << std::endl;
  myCout << " Number of node \t" << Th.nv << std::endl;
  myCout << " Number of element \t" << Th.nt << std::endl;

  myCout << " Time Step \t" << dT << std::endl;
  myCout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  myCout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  myCout << " number of quadrature points in time : \t" << nbTime << std::endl;
  myCout << " number of dof in time per time slab : \t" << ndfTime << std::endl;



  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveSurfactantVTK = true;
  const int frequencyPlottingTheSolution = 10;


  // Space for the velocity field
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange2 FEvelocity(2);
  FESpace2 Vh(Th, FEvelocity);
  Fun_h vel(Vh, fun_velocity);

  // levelSet stuff
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
  FESpace2 Lh  (Th, DataFE<Mesh2>::P1);
  FESpace2 Lh_k(Th, DataFE<Mesh2>::P2);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  for(int i=0;i<nbTime;++i) ls_k[i].init(Lh_k, fun_levelSet);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);
  projection(ls_k[0], ls[nbTime-1]);

  CReinitialization<Mesh> reinitialization;
  reinitialization.number_iteration = 2;
  reinitialization.epsilon_diffusion = 1e-3;
  reinitialization.dt = dT/4;
  reinitialization.ON_OFF = "OFF";
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "OFF";
  reinitialization.precision_correction = 1e-6;
  reinitialization.max_iteration_correction = 10;
  reinitialization.info();
  const int frequencyReinitialization = 5;


  // Declaration of the interface
  TimeInterface2 interface(nbTime);


  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh2> bulkSurface(qTime);

  CutFEM_Parameter mu("mu",1.,1.);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",1.,1.);
  const R sigma0 = 1.;
  const R beta   = 0.25;
  const R epsilon_surfactant = 1.;
  const R Bb = 1.;
  const R Bs = 1.;
  const R Bbs = 1.;
  const R Kb = 0.01;
  const R Ks = 1.;
  const R Kbs = 0.;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  myCout << " \n Beginning of the time iteration \n"
         << " --------------------------------------- \n " << std::endl;

  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {
    double tbegin_iter = CPUtime();
    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;

    double t_ls0 = CPUtime();
    ls_k.begin()->swap(ls_k[lastQuadTime]);
    for(int i=0;i<nbTime;++i) {

      projection(ls_k[i], ls[i]);
      if(iter%frequencyReinitialization == 0 && i == 1) {
        reinitialization.perform(ls_k[i], ls[i]);
      }
      interface.init(i,Th,ls[i].v);

      if(i<nbTime-1) {
        LevelSet levelSet(ls_k[i], vel, vel, dt_levelSet);
        ls_k[i+1].init(levelSet.rhs);
      }
    }
    // std::cout << " Time moving interface \t" << CPUtime() - t_ls0 << std::endl;

    CutFESpace2 Uh(Vh, interface, {1,-1});
    CutFESpace2 Ph(Lh, interface, {1,-1});
    CutFESpace2 Wh(Lh, interface, {1});

    Mesh2 cutThTime(interface);             // surface mesh
    FESpace2 cutWh(cutThTime, interface, DataFE<Mesh2>::P1);   // FE for surfactant
    cutWh.backSpace = &Lh;  // svae backSpace to save solution

/*
                PROBLEM DEFINITION
*/
    Normal n;
    FunTest u(Uh,d), v(Uh,d);             // velocity functions
    FunTest p(Ph,1), q(Ph,1);             // pressure functions
    FunTest ub(Wh,1), vb(Wh,1);           // bulk functions
    FunTest us(cutWh,1), vs(cutWh,1);     // surface functions

    // Connect spaces to problem
    bulkSurface.initSpace(Uh, In);
    bulkSurface.add(Ph, In);
    bulkSurface.add(Wh, In);
    bulkSurface.add(cutWh, In);

    myCout << " DOF       \t : \t" << bulkSurface.nDoF << std::endl;

    // initialize the solution
    int idx_p0 = Uh.NbDoF()*In.NbDoF();
    int idx_b0 = idx_p0 + Ph.NbDoF()*In.NbDoF();
    int idx_s0 = idx_b0 + Wh.NbDoF()*In.NbDoF();
    Rn data_sol0(bulkSurface.nDoF, 0.);
    KN_<double> data_u0(data_sol0(SubArray(Uh.NbDoF()   ,0)));
    KN_<double> data_b0(data_sol0(SubArray(Wh.NbDoF()   ,idx_b0)));
    KN_<double> data_s0(data_sol0(SubArray(cutWh.NbDoF(),idx_s0)));

    if(iter == 0) {
      interpolate(Uh   , data_u0, fun_velocity);
      interpolate(Wh   , data_b0, fun_init_bulk);
      interpolate(cutWh, data_s0, fun_init_surface);
    } else {
      bulkSurface.initialSolution(data_sol0);
    }
    Fun_h u0h(Uh   , data_u0);
    Fun_h b0h(Wh   , data_b0);
    Fun_h s0h(cutWh, data_s0);

    // set uh for newton
    Rn data_sol(data_sol0);
    KN_<double> data_uh(data_sol(SubArray(   Uh.NbDoF()*In.NbDoF(),0)));
    KN_<double> data_ph(data_sol(SubArray(   Ph.NbDoF()*In.NbDoF(),idx_p0)));
    KN_<double> data_bh(data_sol(SubArray(   Wh.NbDoF()*In.NbDoF(),idx_b0)));
    KN_<double> data_sh(data_sol(SubArray(cutWh.NbDoF()*In.NbDoF(),idx_s0)));
    Fun_h uh(Uh   , In, data_uh);
    Fun_h ph(Ph   , In, data_ph);
    Fun_h bh(Wh   , In, data_bh);
    Fun_h sh(cutWh, In, data_sh);


    int iterNewton  = 0;
    while(1) {
      /*
      CONSTRUCTION BILINEAR PART
      */
      // std::cout << "\n Begin assembly " << std::endl;
      double tt0 = CPUtime();
    //   if(iterNewton == 0){
    //     bulkSurface.addBilinear(
    //       innerProduct(Bb*dt(ub), vb)
    //       +innerProduct((vel.expression(2),grad(ub)), vb)*Bb
    //       +innerProduct(grad(ub), grad(vb))*Bb*Kb
    //       , In
    //     );
    //     bulkSurface.addBilinear(
    //       innerProduct(Bs*dt(us), vs)
    //       +innerProduct((vel.expression(),grad(us)), vs)*Bs
    //       +innerProduct(divS(vel)*us, vs)*Bs
    //       +innerProduct(gradS(us), gradS(vs))*Ks*Bs
    //       +innerProduct(Bb*ub-Bs*us, Bb*vb-Bs*vs)
    //       , interface
    //       , In
    //     );
    //
    //     R h3 = pow(meshSize,3);
    //     bulkSurface.addFaceStabilization(
    //         innerProduct(1e-2*h_E*jump(grad(ub)*n), jump(grad(vb)*n))
    //       ,In
    //     );
    //     bulkSurface.addEdgeIntegral(
    //       innerProduct(1e-2*h_E*jump(grad(us)*n), jump(grad(vs)*n))
    //       ,In
    //     );
    //
    //     // impose initial condition
    //     bulkSurface.addBilinear(
    //       innerProduct(ub,vb)*Bb
    //     );
    //     bulkSurface.addBilinear(
    //       innerProduct(us, vs)*Bs
    //       , *interface[0]
    //     );
    //
    //   }
    //   // std::cout << " Time Full assembly \t" << CPUtime() - tt0 << std::endl;
    //   bulkSurface.addMatMul(data_uh);
    //   // std::cout << " Time  A*u0 \t\t" << CPUtime() - tt0 << std::endl;
    //
    //   /*
    //   CONSTRUCTION LINEAR PART
    //   */
    //   tt0 = CPUtime();
    //   // impose initial condition
    //   bulkSurface.addLinear(
    //     -innerProduct(b0h.expression(), vb)*Bb
    //   );
    //   bulkSurface.addLinear (
    //     -innerProduct(s0h.expression(), vs)*Bs
    //     , *interface[0]
    //   );
    //   // std::cout << " Time assembly rhs \t" << CPUtime() - tt0 << std::endl;
    //
    //   /*
    //   CONSTRUCTION NON LINEAR PART
    //   */
    //   bulkSurface.saveMatrix();
    //   tt0 = CPUtime();
    //   ExpressionFunFEM<Mesh2> wb(bh,0,op_id);
    //   ExpressionFunFEM<Mesh2> ws(sh,0,op_id);
    //
    //   bulkSurface.addBilinear(
    //     - innerProduct(wb*us , Bb*vb-Bs*vs)*Bbs
    //     - innerProduct(ub*ws , Bb*vb-Bs*vs)*Bbs
    //     , interface
    //     , In
    //   );
    //
    //   bulkSurface.addLinear(
    //     - innerProduct(wb*ws , Bb*vb-Bs*vs)*Bbs
    //     , interface
    //     , In
    //   );
    //   // std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;
    // //
    // //   stokes.addLagrangeMultiplier(
    // //     innerProduct(1.,dp1), 0.
    // //     , In
    // //   );
    // //   stokes.addLagrangeMultiplier(
    // //     innerProduct(1.,dp1), 0.
    // //     , 0
    // //   );
    // //   // stokes.addLagrangeMultiplier(
    // //   //   innerProduct(1.,dp1), 0.
    // //   //   , 1
    // //   // );
    // //   // stokes.addLagrangeMultiplier(
    // //   //   innerProduct(1.,dp1), 0.
    // //   //   , 2
    // //   // );
    //
    //   bulkSurface.solve();
    //
    //
      KN_<double> dwb(bulkSurface.rhs(SubArray(Wh.NbDoF()*In.NbDoF()   , 0)));
      KN_<double> dws(bulkSurface.rhs(SubArray(cutWh.NbDoF()*In.NbDoF(),idx_s0)));
      double dist  = dwb.l1()/dwb.size();
      double dists = dws.l1()/dws.size();
      myCout << " [dwb , dws] = [ " << dist << " , " << dists << " ] " << std::endl;

      KN_<double> dw(bulkSurface.rhs(SubArray(bulkSurface.nDoF, 0)));
      data_uh -= dw;

      // Size has been changed by the lagrange multiplier
      bulkSurface.rhs.resize(bulkSurface.nDoF); bulkSurface.rhs = 0.0;

      iterNewton+=1;
      if(iterNewton == 1 || max(dist, dists) < 1e-10 ) {

        bulkSurface.saveSolution(data_uh);
        bulkSurface.cleanMatrix();
        break;
      }
    }
    // // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    // {
    //   Fun_h funX(Wh, fun_x);
    //   Fun_h fun1(Wh, fun_1);
    //   R tt0 = CPUtime();
    //   double areaBubble   = integral(fun1, 0, 1) ;
    //   double centerOfMass = integral(funX, 1, 1) / areaBubble ;
    //
    //   double Pb = integralSurf(fun1, 1);
    //   double ra = sqrt(areaBubble / M_PI);
    //   double Pa = 2*M_PI*ra;
    //   double circularity = Pa / Pb;
    //   double riseVelocity = integral(uh, 1, 1) / areaBubble;
    //
    //   double q = integralSurf(us, 0);
    //
    //   std::cout << "\n Features of the drop " << std::endl;
    //   std::cout << " Time                   ->    " << GTime::current_time()+dT  << std::endl;
    //   std::cout << " Center Of Mass         ->    " << centerOfMass << std::endl;
    //   std::cout << " Circularity            ->    " << circularity << std::endl;
    //   std::cout << " Rise velocity          ->    " << riseVelocity << std::endl;
    //   std::cout << " Surfactant quantity    ->    " << q << std::endl;
    //
    //
    //   outputData << GTime::current_time()+dT << "\t"
    //              << centerOfMass << "\t"
    //              << circularity << "\t"
    //              << riseVelocity << "\t"
    //              << areaBubble <<  "\t"
    //              << q << std::endl;
    //
    // }
    //
    //
    // // std::cout << " Set Velocity " << std::endl;
    // for(int i=0;i<nbTime;++i) {
    //   Fun_h sol(Wh,In,data_uh);
    //   set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    // }
    //
    //
    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      std::cout << " Plotting " << std::endl;

      if(saveSurfactantVTK) {
        Paraview2 writerU(Uh, ls[0], pathOutpuFigure+"navierStokes_"+to_string(iterfig)+".vtk");
        Paraview2 writerB(Wh, ls[0], pathOutpuFigure+"bulk_"+to_string(iterfig)+".vtk");
        // writer.add(bh, "bulk", 0, 1);

        Paraview2 writerS(cutWh, ls[0], pathOutpuFigure+"surface_"+to_string(iterfig)+".vtk");
        // writerS.add(sh, "surface", 0, 1);
        // writerS.add(ls[0], "levelSet", 0, 1);

      }

      iterfig++;
    }
    std::cout << " Timer iteration computation \t" << CPUtime() - tbegin_iter << std::endl;
    iter += 1;
    return 0;
  }
  myCout << " -----------------------------------  " << std::endl;
  myCout << " -----------------------------------  " << std::endl;
  myCout << " END OF COMPUTATION " << std::endl;


  // Closing the opened files.
  outputData.close();
}
