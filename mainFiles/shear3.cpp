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
#include "problem/curvature_v2.hpp"
#include "problem/reinitialization.hpp"
#include "time_stuff.hpp"
#include "num/gnuplot.hpp"


/*
Simulation of a rising bubble.
We consider the navier Stokes equations

We are using the levelSet P2 with a reinitializtion step every 5 time iterations
We are using finite element in time, with a quadrature rule using only one point x0 which
is equivalemnt to an implementation of the Euler method for the time discretization

*/


namespace TestBenchmark3D {
  // R fun_levelSet(const R3 P, const int i) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5) + (P.z-0.5)*(P.z-0.5)) - 0.25;}
  // R fun_rhs(const R3 P,const int i) { R3 R(0,0,-0.98); return 0;}//R[i];}
  // R fun_boundary(const R3 P,const int i) { R3 R(0,0,0);  return R[i];}
  R fun_levelSet(const R3 P, const int i) { return sqrt(P.x*P.x + P.y*P.y + P.z*P.z) - 2./3;}
  R fun_boundary(const R3 P, int i) { return (i==0)?0.5*P.z : 0;}

  R fun_x(const R3 P,const int i) {return P[i];}
  R fun_1(const R3 P,const int i) {return 1;}
}


using namespace TestBenchmark3D;


int main(int argc, char** argv )
{

  const int d = 3;
  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface3 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GLevelSet<Mesh> LevelSet;
  typedef GCurvature<Mesh> Curvature;
  // typedef GReinitialization<Mesh, Interface> Reinitialization;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  std::cout << " SIMULATION OF A RISING BUBBLE In 3D. \n"
            << " We consider the navier Stokes equations \n"
            << " We solve the non linear system using the Newton method \n"
            << " We are using the levelSet P1 with a reinitializtion step every 5 time iterations \n"
            << " We are using P0 finite element in time, with a quadrature rule using only one point x0 which is equivalemnt to an implementation of the Euler method for the time discretization "
            << std::endl << std::endl;


  int nx = 20;
  int ny = 10;
  int nz = 10;

  Mesh Th(nx, ny, nz, -2.,-1, -1., 4., 2., 2.);
  double meshSize(4./nx);
  int divisionMeshsize = 2;
  double dT = meshSize/divisionMeshsize;
  double finalTime = 4;
  GTime::total_number_iteration = finalTime/dT;
  dT = finalTime / GTime::total_number_iteration;
  GTime::time_step = dT;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());

  TaylorHood3 FEstokes;
  FESpace Vh(Th, FEstokes);
  FESpace1 Ih(Qh, DataFE<Mesh1>::P0Poly);
  const QuadratureFormular1d& qTime(QF_Euler);
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();


  // Create files to save results
  std::string meshstr = to_string(nx)+to_string(ny)+to_string(nz);
  std::string pathOutpuFolder = "../../outputFiles/shear3/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/shear3/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  std::cout << " path to the output files : "+pathOutpuFolder << std::endl;
  std::cout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;
  myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;


  myCout << " Mesh size \t" << meshSize << std::endl;
  myCout << " Time Step \t" << dT << std::endl;
  myCout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  myCout << " We are using the Euler method (only one quadrature point)" << std::endl;
  myCout << " number of quadrature points in time : \t" << nbTime << std::endl;
  myCout << " number of dof in time per time slab : \t" << ndfTime << std::endl;



  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveStokesVTK = true;
  const bool curvatureVTKFile = true;
  const bool saveLevelSetVTK = true;
  const int frequencyPlottingTheSolution = 1;

  // Space for the velocity field
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange3 FEvelocity(2);
  FESpace velocitySpace(Th, FEvelocity);

  // The initial velocity is 0, so no function needed
  Fun_h vel(velocitySpace);
  // for(int i=0;i<nbTime;++i) vel[i].init(velocitySpace);

  // levelSet stuff
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " We use a P1 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;

  FESpace Lh_k(Th, DataFE<Mesh3>::P1);
  Fun_h ls_k(Lh_k, fun_levelSet);
  FESpace Lh(Th, DataFE<Mesh>::P1);
  Fun_h ls(Lh, fun_levelSet);
  projection(ls_k, ls);


  LevelSet3 levelSet(Lh_k);
  levelSet.setStrongBC({1,2});


  CReinitialization<Mesh> reinitialization;
  reinitialization.number_iteration = 2;
  reinitialization.epsilon_diffusion = 1e-3;
  reinitialization.dt = dT/4;
  reinitialization.ON_OFF = "ON";
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "ON";
  reinitialization.precision_correction = 1e-4;
  reinitialization.max_iteration_correction = 5;
  reinitialization.info();
  const int frequencyReinitialization = 5;//5*divisionMeshsize;


  TimeInterface3 interface(nbTime);

  // FE for curvature problem
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " The curvature is computed using P1 finite elements " << std::endl;
  Lagrange3 FEcurvature(1);


  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh> stokes(qTime);
  // CutFEM_Parameter mu("mu",35,0.35);
  // CutFEM_Parameter rho("rho",1000.,10.);
  // CutFEM_Parameter invmu("invmu",1./35,1./0.35);
  CutFEM_Parameter mu("mu",10,1);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",1e-1,1);

  const R sigma = 24.5;//1.96;//24.5;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG3);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB3);
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  // list<int> dirichlet = {1,3};
  // list<int> neumann = {2,4};
  list<int> dirichlet = {1,2,3,4,5,6};


  // Fun_h fh(Vh, fun_rhs);
  Fun_h gh(Vh, fun_boundary);

  std::cout << "\n Parameters of the problem : \n"
            << " mu_1 = " << mu.val1 << " and mu_2 = "<< mu.val2 << " \n"
            << " rho_1 =  " << rho.val1 << "  and rho_2 = "<< mu.val2 << " \n"
            << " the surface tension sigma = " << sigma
            << std::endl;

            std::cout << " \n Beginning of the time iteration \n"
            << " --------------------------------------- \n " << std::endl;


  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;

    // // projection of the P2 levelSet to the P1 space to compute all
    // // the sub domain and the linear interface
    // projection(ls_k, ls);


    // computation of the nbTime interfaces
    interface.init(0,Th,ls.v);
    CutFESpace3 Wh(Vh, interface, {1,-1});

    /*
                    PROBLEM DEFINITION
    */
    Normal n;
    FunTest du(Wh,d), dp(Wh,1,d), v(Wh,d), q(Wh,1,d), dp1(Wh,1,d,0);
    FunTest Eun = (Eps(du)*n);
    FunTest D2nu = grad(grad(du)*n)*n, D2nv = grad(grad(v)*n)*n;

    stokes.initSpace(Wh, In);

    std::cout << " DOF       \t : \t" << stokes.nDoF << std::endl;

    // Build the Matrix and rhs
    Rn datau0;
    stokes.initialSolution(datau0);
    if(iter == 0) {
      interpolate(Wh, datau0, fun_boundary);
    }

    Rn data_uh(datau0);     // this makes a copy
    Fun_h u0(Wh, datau0);   // the initial solution
    Fun_h uh(Wh, data_uh);  // the initial guess


    // loop for Newton
    int iterNewton  = 0;
    while(1) {

      std::cout << " Begin assembly " << std::endl;
      double tt0 = CPUtime();
      double tt1 = CPUtime();

      if(iterNewton == 0){

        stokes.addBilinear(
          contractProduct(2*mu*Eps(du),Eps(v))
          - innerProduct(dp, div(v)) + innerProduct(div(du), q)
          , In
        );

        // std::cout << " Time Full assembly " << CPUtime() - tt0 << std::endl;
        // tt0 = CPUtime();

        stokes.addBilinear(
          innerProduct(jump(du), -2*mu*average1(Eps(v)*n))
          + innerProduct(-2*mu*average1(Eps(du)*n), jump(v))
          + innerProduct(lambdaG*jump(du), jump(v))
          + innerProduct(average1(dp), jump(v.t()*n))
          - innerProduct(jump(du.t()*n), average1(q))
          , interface
          , In
        );
        // std::cout << " Time surface assembly " << CPUtime() - tt0 << std::endl;

        tt0 = CPUtime();
        stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du,v)
          + innerProduct(dp, v.t()*n)
          - innerProduct(du.t()*n, q)
          - innerProduct(2.*mu*Eps(du)*n, v)
          - innerProduct(du, 2.*mu*Eps(v)*n)
          , In
        );

        tt0 = CPUtime();
        R h3 = pow(meshSize,3);
        stokes.addFaceStabilization(
          innerProduct(1e-2*h_E*jump(grad(du)*n), rho*mu*jump(grad(v)*n))
          // + innerProduct(1e-2*h3*jump(D2nu), rho*mu*jump(D2nv))
          ,In
        );
        stokes.addFaceStabilization(
          innerProduct(1e-2*h3*jump(grad(dp).t()*n), invmu*jump(grad(q).t()*n))
          ,In
        );

        // impose initial condition
        stokes.addBilinear(
          innerProduct(rho*du,v)
        );
      }

      // to optimize the Newton iteration.
      stokes.saveMatrix();
      tt0 = CPUtime();
      FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      ExpressionFunFEM<Mesh3> u1(uh,0,op_id,0), u2(uh,1,op_id,0);

      stokes.addBilinear(
        innerProduct(du1*dx(u1) + du2*dy(u1), rho*v1)
        + innerProduct(du1*dx(u2) + du2*dy(u2), rho*v2)
        + innerProduct(u1*dx(du1) + u2*dy(du1), rho*v1)
        + innerProduct(u1*dx(du2) + u2*dy(du2), rho*v2)
        ,In
      );
      stokes.addLinear(
          innerProduct(u1*dx(u1) + u2*dy(u1), rho*v1)
        + innerProduct(u1*dx(u2) + u2*dy(u2), rho*v2)
        , In
      );
      std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;


    std::cout << " Time matrix assembly " << CPUtime() - tt1 << std::endl;



    tt0 = CPUtime();
    stokes.addMatMul(data_uh);
    std::cout << " Time RHS assembly A*u0 " << CPUtime() - tt0 << std::endl;

    // tt0 = CPUtime();
    stokes.addLinearFormBorder(
      - innerProduct(gh.expression(3),lambdaB*v)
      + innerProduct(gh.expression(3),2.*mu*Eps(v)*n)
      + innerProduct(gh.expression(3), q*n)
      , In
    );


    // tt0 = CPUtime();
    // stokes.addLinear(
    //   -innerProduct(fh.expression(3),rho*v)
    //   , In
    // );
    // std::cout << " Time frhs assembly " << CPUtime() - tt0 << std::endl;


    // tt0 = CPUtime();

      Mesh cutTh(interface);
      FESpace cutVh(cutTh, FEcurvature);
      Curvature curvature(cutVh, interface[0]);//, mapping);
      Fun_h Kx(cutVh, curvature.rhs);
      stokes.addLinear(
         - innerProduct(Kx.expression(3), average2(v))*sigma
        , interface
        , In
      );
      // std::cout << " Time curvature assembly " << CPUtime() - tt0 << std::endl;
      // if(iter%frequencyPlottingTheSolution == 0 && curvatureVTKFile && MPIcf::IamMaster()){
      //   Fun_h solMC(cutVh, curvature.rhs);
      //   Paraview3 writerLS(cutVh, pathOutpuFigure+"curvature_"+to_string(iterfig)+".vtk");
      //   writerLS.add(solMC,"meanCurvature", 0, 2);
      //   writerLS.add(ls, "levelSet", 0, 1);
      // }
      // impose initial condition
      // stokes.addLinear(
      //   -innerProduct(u0.expression(3),rho*v)
      // );

      std::cout << " Time assembly rhs " << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();

      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , In
      );

      std::cout << " Time Lagrange Mult assembly " << CPUtime() - tt0 << std::endl;


      stokes.solve();
      KN_<double> dw(stokes.rhs(SubArray(stokes.nDoF, 0)));
      double dist  = dw.l1()/dw.size();
      std::cout << " dw \t" << dist << std::endl;
      data_uh -= dw;

      stokes.rhs.resize(stokes.nDoF); stokes.rhs = 0.0;

      iterNewton+=1;
      if(iterNewton == 2 || dist < 1e-7 ) {
        stokes.saveSolution(data_uh);
        stokes.cleanMatrix();
        break;
      }

    }
    //
    // // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    // {
    //   Fun_h funX(Wh, fun_x);
    //   Fun_h fun1(Wh, fun_1);
    //   R tt0 = CPUtime();
    //   double areaBubble   = integralCut(fun1, 0, 1) ;
    //   double centerOfMass = integralCut(funX, 1, 1) / areaBubble ;
    //
    //   double Pb = integralSurf(fun1, 1);
    //   double ra = sqrt(areaBubble / M_PI);
    //   double Pa = 2*M_PI*ra;
    //   double circularity = Pa / Pb;
    //
    //   double riseVelocity = integralCut(u0, 1, 1) / areaBubble;
    //
    //   std::cout << " Center Of Mass    ->    " << centerOfMass << std::endl;
    //   std::cout << " Circularity       ->    " << circularity << std::endl;
    //   std::cout << " Rise velocity     ->    " << riseVelocity << std::endl;
    //
    //   outputData << GTime::current_time() << "\t"
    //              << centerOfMass << "\t"
    //              << circularity << "\t"
    //              << riseVelocity << "\t"
    //              << areaBubble <<  std::endl;
    // }
    //


    //               PLOTTING
    if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      std::cout << " Plotting " << std::endl;

      if(saveStokesVTK) {
        Fun_h sol(Wh, data_uh);
        Paraview3 writerr(Wh, ls, pathOutpuFigure+"NavierStokes_"+to_string(iterfig)+".vtk");
        writerr.add(sol, "velocity", 0, 3);
        writerr.add(sol, "pressure", 3, 1);
        writerr.add(ls,"levelSet", 0 , 1);
      }
      iterfig++;
    }

    {
      std::cout << " Set Velocity " << std::endl;
      Fun_h sol(Wh,In,data_uh);
      set_velocity(sol, vel, ls, In.map(qTime[0]));
    }

    // MOVE THE LEVELSET
    {
      std::cout << " Move LevelSet " << std::endl;
      double dt_levelSet = dT/2;
      levelSet.solve(ls_k, vel, vel, dt_levelSet);
      ls_k.init(levelSet.rhs);
      // if(iter%frequencyReinitialization == 0){
      //   projection(ls_k, ls);
      //   reinitialization.perform(ls_k, ls_k);
      // }
      levelSet.solve(ls_k, vel, vel, dt_levelSet);
      ls_k.init(levelSet.rhs);
      projection(ls_k, ls);
    }





    iter += 1;
  }

  std::cout << " -----------------------------------  " << std::endl;
  std::cout << " -----------------------------------  " << std::endl;
  std::cout << " END OF COMPUTATION " << std::endl;
  std::cout << " Total time  \t" << CPUtime() - cpubegin << std::endl;

  // Closing the opened files.
  outputData.close();

}
