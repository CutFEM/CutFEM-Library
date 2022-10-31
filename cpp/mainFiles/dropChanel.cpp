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
#define FORMULATION 1;
// #define FORMULATION 2;
namespace ExtensionPairBubble2D {
  R fun_levelSet1(const R2 P, const int i) { return sqrt(P.x*P.x + (P.y-1.201)*(P.y-1.201)) - 1.;}
  R fun_levelSet2(const R2 P, const int i) { return sqrt(P.x*P.x + (P.y+1.201)*(P.y+1.201)) - 1.;}

  R fun_levelSet(const R2 P, const int i) {
    R2 wc(3.5,0.5);
    return sqrt((P.x-wc.x)*(P.x-wc.x) + (P.y-wc.y)*(P.y-wc.y)) - 0.66;
  }
  R fun_rhs(const R2 P,const int i) { return 0;}
  R fun_boundary(const R2 P,const int i) { double Q = 1.;
    // return (i==0)?Q*(2-P.y)*(2+P.y) : 0 ;}
    return (i==0)?Q*(2-P.y)*(1+P.y) : 0 ;}

  R fun_init_surfactant(const R2 P, int i) { return 1.;}

  R fun_x(const R2 P,const int i) {return P[i];}
  R fun_1(const R2 P,const int i) {return 1;}
}


using namespace ExtensionPairBubble2D;


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

  // int nx = 30;
  // int ny = 30;
  // Mesh2 Th(nx, ny, -3., -3., 6., 6.);
  // double meshSize = 6./nx;
  // double meshSize(6./nx);
  Mesh2 Th("../mesh/curveChanel3.msh");
  int nx = 150;
  int ny = 40;
  double meshSize = 6./nx;

  int divisionMeshsize = 2;
  double dT = meshSize/divisionMeshsize;
  double final_time = 1;
  GTime::total_number_iteration = final_time/dT;
  dT = final_time / GTime::total_number_iteration;
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
  std::string pathOutpuFolder = "../../outputFiles/dropInCurveChanel/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/dropInCurveChanel/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;


  myCout << " SIMULATION OF BUBBLE IN CHANEL. \n"
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
  const bool saveSurfactantVTK = false;
  const bool curvatureVTKFile = false;
  const bool saveLevelSetVTK = true;
  const int frequencyPlottingTheSolution = 1;


  // Space for the velocity field
  // ----------------------------------------------
  myCout << "\n --------------------------------------- \n " << std::endl;
  myCout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange2 FEvelocity(2);
  FESpace2 VelVh(Th, FEvelocity);
  Fun_h vel(Vh, fun_boundary);

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

  LevelSet2 levelSet(Lh_k);
  levelSet.setStrongBC({1});


  CReinitialization<Mesh> reinitialization;
  reinitialization.number_iteration = 4;
  reinitialization.epsilon_diffusion = 1e-3;
  reinitialization.dt = dT/8;
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
  CutFEM_Parameter mu("mu",10.,2.);
  CutFEM_Parameter rho("rho",1.,1.);
  CutFEM_Parameter invmu("invmu",0.1,1./2);

  const R epsilon_surfactant = 0.1;
  const R sigma0 = 50.;
  const R beta   = 0.5;
  double initialConcentrationSurfactant;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  list<int> inFlow = {1,3};
  list<int> dirichlet = {1,2,3};

  list<int> neumann = {2};

  // interpolate rhs and boundary condition. Not time dependent and not domain dependent
  Fun_h gh(Vh, fun_boundary);


  myCout << "\n Parameters of the problem : \n"
            << " mu_1 = " << mu.val1 << " and mu_2 = "<< mu.val2 << " \n"
            << " rho_1 = "<< rho.val1 << " and rho_2 = " << rho.val2 << " \n"
            << " the surface tension sigma = " << sigma0 << " \n"
            << " beta = " << beta
            << std::endl;

            myCout << " \n Beginning of the time iteration \n"
            << " --------------------------------------- \n " << std::endl;

  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;

    ls_k.begin()->swap(ls_k[lastQuadTime]);
    for(int i=0;i<nbTime;++i) {

      projection(ls_k[i], ls[i]);
      if(iter%frequencyReinitialization == 0 && i == 1 && iter > 0) {
        reinitialization.perform(ls_k[i], ls[i]);
      }
      interface.init(i,Th,ls[i].v);

      if(i<nbTime-1) {
        levelSet.solve(ls_k[i], vel, vel, dt_levelSet);
        ls_k[i+1].init(levelSet.rhs);
      }
    }

    CutFESpace2 Wh(Vh, interface, {1,-1});
    // Mesh2 cutThTime(interface);             // surface mesh
    // FESpace2 cutWh(cutThTime, interface, DataFE<Mesh2>::P1);   // FE for surfactant
    // cutWh.backSpace = &Lh;  // svae backSpace to save solution

/*
                PROBLEM DEFINITION
*/
    Normal n;
    FunTest du(Wh,d), dp(Wh,1,d), v(Wh,d), q(Wh,1,d), dp1(Wh,1,d,0);
    FunTest Eun = (Eps(du)*n);
    FunTest D2nu = grad(grad(du)*n)*n, D2nv = grad(grad(v)*n)*n;
    // FunTest ds(cutWh,1), r(cutWh,1);

    // Connect spaces and interfaces to problem
    stokes.initSpace(Wh, In);
    // stokes.add(cutWh, In);
    std::cout << " DOF       \t : \t" << stokes.nDoF << std::endl;


    // initialize the solution
    int idx_s0 = Wh.NbDoF()*In.NbDoF();
    Rn data_sol0(stokes.nDoF, 0.);
    KN_<double> data_u0(data_sol0(SubArray(Wh.NbDoF()   ,0)));
    // KN_<double> data_s0(data_sol0(SubArray(cutWh.NbDoF(),idx_s0)));
    if(iter == 0) {
      // interpolate(cutWh, data_s0, fun_init_surfactant);
    } else {
      stokes.initialSolution(data_sol0);
    }
    Fun_h u0h(Wh   , data_u0);
    // Fun_h s0h(cutWh, data_s0);

    // set uh and initial solution for N-S
    Rn data_sol(data_sol0);
    KN_<double> data_uh(data_sol(SubArray(   Wh.NbDoF()*In.NbDoF(),0)));
    // KN_<double> data_sh(data_sol(SubArray(cutWh.NbDoF()*In.NbDoF(),idx_s0)));
    Fun_h uh(   Wh, In, data_uh);
    // Fun_h us(cutWh, In, data_sh);

    // loop for Newton
    int iterNewton  = 0;
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
          , interface
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
        // stokes.addBilinearFormBorder(
        //   innerProduct(lambdaB*du.t()*n,v.t()*n)
        //   + innerProduct(dp, v.t()*n)
        //   - innerProduct(du.t()*n, q)
        //   - innerProduct(2.*mu*Eun.t()*n, v.t()*n)
        //   - innerProduct(du.t()*n, 2.*mu*Eun.t()*n)
        //   , In
        //   , neumann
        // );

        // #ifdef FORMULATION 1
        // stokes.addBilinear(
        //   innerProduct(dt(ds), r)
        // + innerProduct(epsilon_surfactant*gradS(ds), gradS(r))
        //   , interface
        //   , In
        // );
        // double cc = 1./In.T.mesure()*6.;
        // stokes.addBilinear(
        //   innerProduct(cc*ds, r)
        //   , interface
        //   , 0, In
        // );
        // #else
        // stokes.addBilinear(
        //   - innerProduct(ds, dt(r))
        //   + innerProduct(epsilon_surfactant*gradS(ds), gradS(r))
        //   , interface
        //   , In
        // );
        // double cc = 1./In.T.mesure()*6.;
        // stokes.addBilinear(
        //   innerProduct(cc*ds, r)
        //   , interface
        //   , lastQuadTime, In
        // );
        // #endif
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
        // stokes.addEdgeIntegral(
        //   innerProduct(1e-2*jump(grad(ds).t()*n), jump(grad(r).t()*n))
        //   ,In
        // );
        // if stab interface, can have h in both stabilization
        // stokes.addBilinear(
        //   innerProduct(1e-2*h_E*grad(ds).t()*n, grad(r).t()*n)
        //   , interface
        //   , In
        // );


        // impose initial condition
        stokes.addBilinear(
          innerProduct(rho*du,v)
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

        // if(iterNewton == 0) {
        //   // sigma(w) = 24,5(1 - Beta*w)
        //   // -(sigma(w)k , <v.n>) = -(s0*Bw, <v.n>)
        //   // +(dSsigma(w) , <v> ) =
        //   stokes.addBilinear(
        //       innerProduct(Kx*ds    , average2(vx))*beta*sigma0
        //     + innerProduct(Ky*ds    , average2(vy))*beta*sigma0
        //     + innerProduct(gradS(ds), average2(v) )*beta*sigma0
        //     , interface
        //     ,i , In
        //   );
        // }
        stokes.addLinear(
            - innerProduct(Kx  , average2(vx))*sigma0
            - innerProduct(Ky  , average2(vy))*sigma0
            , interface
            , i , In
        );
        // if(iter%frequencyPlottingTheSolution == 0 && curvatureVTKFile && i == 0){
        //   Fun_h solMC(cutVh, curvature.rhs);
        //   Paraview2 writerLS(cutVh, pathOutpuFigure+"curvature_"+to_string(iterfig)+".vtk");
        //   writerLS.add(solMC, "meanCurvature", 0, 2);
        //   writerLS.add(ls[0], "levelSet", 0, 1);
        // }
      }
      std::cout << " Time Full assembly \t" << CPUtime() - tt0 << std::endl;

      /*
      CONSTRUCTION LINEAR PART
      */
      tt0 = CPUtime();
      // impose initial condition
      stokes.addLinear(
        -innerProduct(u0h.expression(2),rho*v)
      );
      // double cc = 1./In.T.mesure()*6.;
      // stokes.addLinear (
      //   -innerProduct(s0h.expression(), cc*r)
      //   , interface
      //   , 0, In
      // );

      stokes.addLinearFormBorder(
        - innerProduct(gh.expression(2),lambdaB*v)
        + innerProduct(gh.expression(2),2.*mu*Eps(v)*n)
        + innerProduct(gh.expression(2), q*n)
        , In
        , inFlow
      );

      std::cout << " Time assembly rhs \t" << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();

      stokes.addMatMul(data_sol);
      std::cout << " Time  A*u0 \t\t" << CPUtime() - tt0 << std::endl;



      /*
      CONSTRUCTION NON LINEAR PART
      */

      // to optimize the Newton iteration.
      // stokes.saveMatrix();
      // tt0 = CPUtime();
      // FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      // ExpressionFunFEM<Mesh2> u1(uh,0,op_id,0), u2(uh,1,op_id,0);
      // // ExpressionFunFEM<Mesh2> s(us,0,op_id);
      //
      // stokes.addBilinear(
      //   innerProduct(du1*dx(u1) + du2*dy(u1), rho*v1)
      //   + innerProduct(du1*dx(u2) + du2*dy(u2), rho*v2)
      //   + innerProduct(u1*dx(du1) + u2*dy(du1), rho*v1)
      //   + innerProduct(u1*dx(du2) + u2*dy(du2), rho*v2)
      //   ,In
      // );
      // stokes.addLinear(
      //     innerProduct(u1*dx(u1) + u2*dy(u1), rho*v1)
      //   + innerProduct(u1*dx(u2) + u2*dy(u2), rho*v2)
      //   , In
      // );
//
// #ifdef FORMULATION2
//       stokes.addBilinear(
//         - innerProduct(ds*u1 , dx(r)) - innerProduct(ds*u2 , dy(r))
//         - innerProduct(s*du , grad(r))
//         , interface
//         , In
//       );
//       stokes.addLinear(
//         - innerProduct(s*u1 , dx(r))
//         - innerProduct(s*u2 , dy(r))
//         , interface
//         , In
//       );
// #else
//       stokes.addBilinear(
//           innerProduct(u1  *dx(ds) + u2  *dy(ds) , r)
//         + innerProduct(du1*dx(s)  + du2*dy(s)  , r)
//         + innerProduct(ds*divS(uh)   , r)
//         + innerProduct(s *dxS(du1) + s *dyS(du2) , r)
//         , interface
//         , In
//       );
//       stokes.addLinear(
//           innerProduct(u1*dx(s)   + u2*dy(s)  , r)
//         + innerProduct(s *dxS(uh) + s *dyS(uh), r)
//         , interface
//         , In
//       );
//
// #endif
      std::cout << " Time assembly NL \t" << CPUtime() - tt0 << std::endl;
      tt0 = CPUtime();

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
      // KN_<double> dws(stokes.rhs(SubArray(cutWh.NbDoF()*In.NbDoF(),idx_s0)));
      double dist  = dwu.l1()/dwu.size();
      double dists = 0;//dws.l1()/dws.size();
      std::cout << " dwu\t" << dist << std::endl;
      // std::cout << " dws\t" << dists << std::endl;

      KN_<double> dw(stokes.rhs(SubArray(stokes.nDoF, 0)));
      data_sol -= dw;

      // Size has been changed by the lagrange multiplier
      stokes.rhs.resize(stokes.nDoF); stokes.rhs = 0.0;

      iterNewton+=1;

      if(iterNewton == 1 || max(dist, dists) < 1e-10  ) {
        stokes.saveSolution(data_sol);
        stokes.cleanMatrix();
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
    //   double q    = integralSurf(us,  0, 0 , In.map(qTime[0]));
    //   double qend = integralSurf(us,  0, lastQuadTime, In.map(qTime[lastQuadTime]));
    //   double q0   = integralSurf(u0s, 0);
    //
    //
    //   if(iter == 0) initialConcentrationSurfactant = qend;
    //   std::cout << "\n Features of the drop " << std::endl;
    //   std::cout << " Time                   ->    " << GTime::current_time()+dT  << std::endl;
    //   std::cout << " Center Of Mass         ->    " << centerOfMass << std::endl;
    //   std::cout << " Circularity            ->    " << circularity << std::endl;
    //   std::cout << " Rise velocity          ->    " << riseVelocity << std::endl;
    //   std::cout << " Surfactant quantity init->    " << q0 << std::endl;
    //   std::cout << " Surfactant quantity    ->    " << q << std::endl;
    //   std::cout << " Surfactant quantity end ->   " << qend << std::endl;
    //   std::cout << " |q_0 - q_end|           ->   " << fabs(q - qend) << std::endl;
    //   std::cout << " Surfactant conservation ->   " << fabs(qend - initialConcentrationSurfactant) << std::endl;
    //
    //
    //   outputData << GTime::current_time()+dT << "\t"
    //              << centerOfMass << "\t"
    //              << circularity << "\t"
    //              << riseVelocity << "\t"
    //              << areaBubble <<  "\t"
    //              << qend << "\t"
    //              << fabs(q - qend) << "\t"
    //              << fabs(qend - initialConcentrationSurfactant) << std::endl;
    //
    // }

    set_velocity(uh, vel, ls[lastQuadTime], In.map(qTime[lastQuadTime]));



    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      std::cout << " Plotting " << std::endl;

      if(saveStokesVTK) {
        Paraview2 writerr(Wh, ls[0], pathOutpuFigure+"navierStokes_"+to_string(iterfig)+".vtk");
        writerr.add(uh,"velocity", 0, 2);
        writerr.add(uh,"pressure", 2, 1);
      }
      // if(saveSurfactantVTK) {
      //   Paraview2 writer(cutWh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
      //   writer.add(us, "surfactant", 0, 1);
      //   writer.add(ls[0], "levelSet", 0, 1);
      // }
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


  // Closing the opened files.
  outputData.close();
}
