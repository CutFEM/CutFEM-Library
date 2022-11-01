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
  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}
  R fun_rhs(const R2 P,const int i) { R2 R(0,-0.98); return (i<2)?R[i] : 0;}
  R fun_boundary(const R2 P,const int i) { R2 R(0,0);  return (i<2)?R[i] : 0;}

  R fun_init_surfactant(const R2 P, int i) { return 1.;}
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
  int ny = 40;

  Mesh2 Th(nx, ny, 0., 0., 1., 2.);
  double meshSize(1./nx);
  // Mesh2 Th("../mesh/RBadaptMesh1_3060.msh");
  // int nx = 30;
  // int ny = 60;
  // double meshSize((1*2)/sqrt(nx*ny));

  int divisionMeshsize = 2;
  double dT = meshSize/divisionMeshsize;
  GTime::total_number_iteration = 3/dT;
  dT = 3. / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());

  TaylorHood2 FEstokes;
  FESpace2 Vh(Th, FEstokes);
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  // Create files to save results
  std::string meshstr = to_string(nx)+to_string(ny);
  std::string pathOutpuFolder = "../../outputFiles/RBSurfactantSimpson/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/RBSurfactantSimpson/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;


  myCout << " SIMULATION OF A RISING BUBBLE WITH INSOLABLE SURFACTANT. \n"
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
  const bool saveStokesVTK = true;
  const bool saveSurfactantVTK = true;
  const bool curvatureVTKFile = false;
  const bool saveLevelSetVTK = false;
  const int frequencyPlottingTheSolution = 1;


  // Space for the velocity field
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange2 FEvelocity(2);
  FESpace2 VelVh(Th, FEvelocity);

  // The initial velocity is 0, so no function needed
  vector<Fun_h> vel(nbTime);
  for(int i=0;i<nbTime;++i) vel[i].init(VelVh);


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
  reinitialization.ON_OFF = "ON";
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "ON";
  reinitialization.precision_correction = 1e-6;
  reinitialization.max_iteration_correction = 10;
  reinitialization.info();
  const int frequencyReinitialization = 5;


  // Declaration of the interface
  TimeInterface2 interface(nbTime);

  // FE for curvature problem
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " The curvature is computed using second order finite elements " << std::endl;
  Lagrange2 FEcurvature(2);

  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh2> stokes(qTime);
  CutFEM_Parameter mu("mu",10.,1.);
  CutFEM_Parameter rho("rho",1000.,100.);
  CutFEM_Parameter invmu("invmu",1./10,1.);

  const R epsilon_surfactant = 0.1;
  const R sigma0 = 24.5;
  const R beta   = 0.5;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);

  list<int> dirichlet = {1,3};
  list<int> neumann = {2,4};
  myCout << "\n Parameters of the problem : \n"
            << " mu_1 = 10 and mu_2 = 1 \n"
            << " rho_1 = 1000 and rho_2 = 100 \n"
            << " the surface tension sigma = 24,5 \n"
            << " beta = " << beta
            << std::endl;

            myCout << " \n Beginning of the time iteration \n"
            << " --------------------------------------- \n " << std::endl;

  int iter = 0, iterfig = 0;
  int nonConvergingIteration = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;


    // use the last position of the interface during previous iteration
    // as first step in the current iteration
    ls_k.begin()->swap(ls_k[lastQuadTime]);

    // computation of the nbTime interfaces
    for(int i=0;i<nbTime;++i) {

      projection(ls_k[i], ls[i]);

      // here we want to reinitializa the levelSet.
      // !!! WE CANNOT DO ON THE FIRST LEVELSET BECAUSE IT HAS TO MATCH THE LAST STEP
      if(iter%frequencyReinitialization == 0 && i == 1) {
        reinitialization.perform(ls_k[i], ls[i]);
      }
      interface.init(i,Th,ls[i].v);

      if(i<nbTime-1) {
        LevelSet levelSet(ls_k[i], vel[lastQuadTime], vel[lastQuadTime], dt_levelSet);
        ls_k[i+1].init(levelSet.rhs);
      }
    }

    CutFESpace2 Wh(Vh, interface, {1,-1});

    // Create the Active Cut Mesh for insoluble surfactant
    Mesh2 cutThTime(interface);             // surface mesh
    FESpace2 cutWh(cutThTime, interface, DataFE<Mesh2>::P1);   // FE for surfactant
    cutWh.backSpace = &Lh;  // svae backSpace to save solution

/*
                PROBLEM DEFINITION
*/
    Normal n;
    FunTest du(Wh,d), dp(Wh,1,d), v(Wh,d), q(Wh,1,d), dp1(Wh,1,d,0);
    FunTest Eun = (Eps(du)*n);
    FunTest D2nu = grad(grad(du)*n)*n, D2nv = grad(grad(v)*n)*n;
    FunTest ds(cutWh,1), r(cutWh,1);

    // Connect spaces and interfaces to problem
    stokes.initSpace(Wh, In);
    stokes.add(cutWh, In);
    std::cout << " DOF       \t : \t" << stokes.nDoF << std::endl;

    // interpolate rhs
    Fun_h fh(Vh, fun_rhs);

    // initialize the solution
    Rn datau0;
    stokes.initialSolution(datau0);

    // set uh and initial solution for N-S
    Rn data_uh(datau0);
    int n0 = Wh.NbDoF()*In.NbDoF();

    // Initial solution
    Fun_h u0(Wh, datau0);
    KN_<double> datas0(datau0(SubArray(cutWh.NbDoF(),n0)));
    if(iter == 0) {
      interpolate(cutWh, datas0, fun_init_surfactant);
    }
    Fun_h u0s(cutWh, datas0);


    // Set solution
    KN_<double> datas_i(data_uh(SubArray(cutWh.NbDoF()*In.NbDoF(),n0)));
    Fun_h uh(Wh, In, data_uh);
    Fun_h us(cutWh, In, datas_i);

    // u0.print();
    // loop for Newton
    int iterNewton  = 0;
    double previous_delta = 1e300;
    while(1) {
      /*
      CONSTRUCTION BILINEAR PART
      */
      std::cout << "\n Begin assembly " << std::endl;
      double tt0 = CPUtime();
      if(iterNewton == 0){
        stokes.addBilinear(
          innerProduct(dt(du), rho*v)
          + contractProduct(2*mu*Eps(du),Eps(v))
          - innerProduct(dp, div(v)) + innerProduct(div(du), q)
          , In
        );

        stokes.addBilinear(
          innerProduct(jump(du), -2*mu*average1(Eps(v)*n))
          + innerProduct(-2*mu*average1(Eps(du)*n), jump(v))
          + innerProduct(lambdaG*jump(du), jump(v))
          + innerProduct(average1(dp), jump(v.t()*n))
          - innerProduct(jump(du.t()*n), average1(q))
          ,interface
          , In
        );
        stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du,v)
          + innerProduct(dp, v.t()*n)
          - innerProduct(du.t()*n, q)
          - innerProduct(2.*mu*Eps(du)*n, v)
          - innerProduct(du, 2.*mu*Eps(v)*n)
          , In
          , dirichlet
        );
        stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du.t()*n,v.t()*n)
          + innerProduct(dp, v.t()*n)
          - innerProduct(du.t()*n, q)
          - innerProduct(2.*mu*Eun.t()*n, v.t()*n)
          - innerProduct(du.t()*n, 2.*mu*Eun.t()*n)
          , In
          , neumann
        );


        stokes.addBilinear(
            innerProduct(dt(ds), r)
          + innerProduct(epsilon_surfactant*gradS(ds), gradS(r))
          , interface
          , In
        );

        R h3 = pow(meshSize,3);

        stokes.addFaceStabilization(
            innerProduct(1e-2*h *jump(grad(du)*n), rho*mu*jump(grad(v)*n))
          + innerProduct(1e-2*h3*jump(D2nu)      , rho*mu*jump(D2nv))
          ,In
        );
        // // Mass matrix with time term
        // stokes.addFaceStabilization(
        //   innerProduct(1e-2*h3*jump(grad(dt(du))*n), rho*mu*jump(grad(v)*n))
        //   ,In
        // );
        // // Mass matrix
        // stokes.addFaceStabilization(
        //   innerProduct(1e-2*h3*jump(grad(du)*n), rho*jump(grad(v)*n))
        // );
        stokes.addFaceStabilization(
          innerProduct(1e-2*h3*jump(grad(dp).t()*n), invmu*jump(grad(q).t()*n))
          ,In
        );
        stokes.addEdgeIntegral(
          innerProduct(1e-2*h*jump(grad(ds).t()*n), jump(grad(r).t()*n))
          ,In
        );
        // if stab interface, can have h in both stabilization
        stokes.addBilinear(
          innerProduct(1e-2*h*grad(ds).t()*n, grad(r).t()*n)
          , interface
          ,In
        );


        // impose initial condition
        stokes.addBilinear(
          innerProduct(rho*du,v)
        );
        stokes.addBilinear(
          innerProduct(ds, r)
          , *interface[0]
        );

      }
      for(int i=0;i<nbTime;++i) {      // computation of the curvature
        Mesh2 cutTh(*interface[i]);
        FESpace2 cutVh(cutTh, FEcurvature);
        Mapping2 mapping(cutVh, ls_k[i]);

        Curvature curvature(cutVh, *interface[i], mapping);
        Fun_h meanCurvature(cutVh, curvature.rhs);

        ExpressionFunFEM<Mesh2> Kx(meanCurvature,0,op_id);
        ExpressionFunFEM<Mesh2> Ky(meanCurvature,1,op_id);
        FunTest vx(Wh,1,0), vy(Wh,1,1);

        if(iterNewton == 0) {
          // sigma(w) = 24,5(1 - Beta*w)
          // -(sigma(w)k , <v.n>) = -(s0*Bw, <v.n>)
          // +(dSsigma(w) , <v> ) =
          stokes.addBilinear(
              innerProduct(Kx*ds    , average2(vx))*beta*sigma0
            + innerProduct(Ky*ds    , average2(vy))*beta*sigma0
            + innerProduct(gradS(ds), average2(v) )*beta*sigma0
            , interface
            ,i , In
          );
        }
        stokes.addLinear(
            - innerProduct(Kx  , average2(vx))*sigma0
            - innerProduct(Ky  , average2(vy))*sigma0
            , interface
            ,i , In
        );
        if(iter%frequencyPlottingTheSolution == 0 && curvatureVTKFile && i == 0){
          Fun_h solMC(cutVh, curvature.rhs);
          Paraview2 writerLS(cutVh, pathOutpuFigure+"curvature_"+to_string(iterfig)+".vtk");
          writerLS.add(solMC, "meanCurvature", 0, 2);
          writerLS.add(ls[0], "levelSet", 0, 1);
        }
      }
      std::cout << " Time Full assembly \t" << CPUtime() - tt0 << std::endl;

      /*
      CONSTRUCTION LINEAR PART
      */
      tt0 = CPUtime();
      stokes.addLinear(
        -innerProduct(fh.expression(2),rho*v)
        , In
      );
      // impose initial condition
      stokes.addLinear(
        -innerProduct(u0.expression(2),rho*v)
      );
      stokes.addLinear (
        -innerProduct(u0s.expression(), r)
        , *interface[0]
      );
      std::cout << " Time assembly rhs \t" << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();

      stokes.addMatMul(data_uh);
      std::cout << " Time  A*u0 \t\t" << CPUtime() - tt0 << std::endl;



      /*
      CONSTRUCTION NON LINEAR PART
      */

      // to optimize the Newton iteration.
      stokes.saveMatrix();
      tt0 = CPUtime();
      FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      ExpressionFunFEM<Mesh2> u1(uh,0,op_id,0), u2(uh,1,op_id,0);
      ExpressionFunFEM<Mesh2> s(us,0,op_id);

      stokes.addBilinear(
        innerProduct(du1*dx(u1) + du2*dy(u1), rho*v1)
        + innerProduct(du1*dx(u2) + du2*dy(u2), rho*v2)
        +
      	innerProduct(u1*dx(du1) + u2*dy(du1), rho*v1)
        + innerProduct(u1*dx(du2) + u2*dy(du2), rho*v2)
        ,In
      );
      stokes.addLinear(
          innerProduct(u1*dx(u1) + u2*dy(u1), rho*v1)
        + innerProduct(u1*dx(u2) + u2*dy(u2), rho*v2)
        , In
      );
      stokes.addBilinear(
            innerProduct(u1  *dx(ds) + u2  *dy(ds) , r)
          + innerProduct(du1*dx(s)  + du2*dy(s)  , r)
          + innerProduct(ds*divS(uh)   , r)
          + innerProduct(s *dxS(du1) + s *dyS(du2) , r)
          , interface
          , In
      );
      stokes.addLinear(
            innerProduct(u1*dx(s)   + u2*dy(s)  , r)
          + innerProduct(s *dxS(uh) + s *dyS(uh), r)
          , interface
        , In
      );
      std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;

      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , In
      );
      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , 0
      );
      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,dp1), 0.
      //   , 1
      // );
      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,dp1), 0.
      //   , 2
      // );

      stokes.solve();


      KN_<double> dwu(stokes.rhs(SubArray(Wh.NbDoF()*In.NbDoF(), 0)));
      KN_<double> dws(stokes.rhs(SubArray(cutWh.NbDoF()*In.NbDoF(),n0)));
      double dist  = dwu.l1()/dwu.size();
      double dists = dws.l1()/dws.size();
      std::cout << " dwu\t" << dist << std::endl;
      std::cout << " dws\t" << dists << std::endl;


      KN_<double> dw(stokes.rhs(SubArray(stokes.nDoF, 0)));
      data_uh -= dw;

      // Size has been changed by the lagrange multiplier
      stokes.rhs.resize(stokes.nDoF); stokes.rhs = 0.0;

      iterNewton+=1;

      if(iterNewton == 1 || max(dist, dists) < 1e-3 || previous_delta < max(dist, dists)  ) {

        stokes.saveSolution(data_uh);
        stokes.cleanMatrix();
        break;
      }
    }
    // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    {
      Fun_h funX(Wh, fun_x);
      Fun_h fun1(Wh, fun_1);
      R tt0 = CPUtime();
      double areaBubble   = integral(fun1, 0, 1) ;
      double centerOfMass = integral(funX, 1, 1) / areaBubble ;

      double Pb = integralSurf(fun1, 1);
      double ra = sqrt(areaBubble / M_PI);
      double Pa = 2*M_PI*ra;
      double circularity = Pa / Pb;
      double riseVelocity = integral(uh, 1, 1) / areaBubble;

      double q = integralSurf(us, 0);

      std::cout << "\n Features of the drop " << std::endl;
      std::cout << " Time                   ->    " << GTime::current_time()+dT  << std::endl;
      std::cout << " Center Of Mass         ->    " << centerOfMass << std::endl;
      std::cout << " Circularity            ->    " << circularity << std::endl;
      std::cout << " Rise velocity          ->    " << riseVelocity << std::endl;
      std::cout << " Surfactant quantity    ->    " << q << std::endl;


      outputData << GTime::current_time()+dT << "\t"
                 << centerOfMass << "\t"
                 << circularity << "\t"
                 << riseVelocity << "\t"
                 << areaBubble <<  "\t"
                 << q << std::endl;

    }


    // std::cout << " Set Velocity " << std::endl;
    for(int i=0;i<nbTime;++i) {
      Fun_h sol(Wh,In,data_uh);
      set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    }


    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      std::cout << " Plotting " << std::endl;

      if(saveStokesVTK) {
        KN_<double> solv(data_uh(SubArray(Wh.NbDoF(), 0)));
        Fun_h sol(Wh, solv);
        Paraview2 writerr(Wh, ls[0], pathOutpuFigure+"navierStokes_"+to_string(iterfig)+".vtk");
        writerr.add(sol,"velocity", 0, 2);
        writerr.add(sol,"pressure", 2, 1);
      }
      if(saveSurfactantVTK) {

        // Rn sols(cutWh.NbDoF(), 0.);
        //
        // for(int i=0;i<Ih[0].NbDoF();++i){
        //   int ndf = n0;
        //   KN_<double> solsi(uh(SubArray(cutWh.NbDoF(), ndf)));
        //   ndf += cutWh.NbDoF();
        //   sols += solsi;
        // }
        // Fun2 sol(cutWh, sols);

        KN_<double> sols(data_uh(SubArray(cutWh.NbDoF(), n0)));
        Fun_h sol(cutWh, sols);

        Paraview2 writer(cutWh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
        writer.add(sol, "surfactant", 0, 1);
        writer.add(ls[0], "levelSet", 0, 1);
        // writer.add(ls[2], "levelSet1", 0, 1);

      }
      if(saveLevelSetVTK) {
        Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
        writerLS.add(ls[0],"levelSet", 0, 1);
      }

      iterfig++;
    }

    iter += 1;
  }
  myCout << " -----------------------------------  " << std::endl;
  myCout << " -----------------------------------  " << std::endl;
  myCout << " END OF COMPUTATION " << std::endl;
  myCout << " Non converging iteration : \t" << nonConvergingIteration << std::endl;


  // Closing the opened files.
  outputData.close();
}