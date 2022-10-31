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
We consider the navier Stokes equations

We are using the levelSet P2 with a reinitializtion step every 5 time iterations
We are using finite element in time, with the simpson rule in time.
*/
namespace TestBenchmark2D {
  R fun_levelSet(const R2 P, const int i) { return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.5)*(P.y-0.5)) - 0.25;}
  R fun_rhs(const R2 P,const int i) { R2 R(0,-0.98); return (i<2)?R[i] : 0;}
  R fun_boundary(const R2 P,const int i) { R2 R(0,0);  return (i<2)?R[i] : 0;}

  R fun_x(const R2 P,const int i) {return P[i];}
  R fun_1(const R2 P,const int i) {return 1;}
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
  typedef LevelSet2 LevelSet;
  typedef GCurvature<Mesh> Curvature;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();
  std::cout << " SIMULATION OF A RISING BUBBLE. \n"
            << " We consider the navier Stokes equations \n"
            << " We solve the non linear system using the Newton method \n"
            << " We are using the levelSet P2 with a reinitializtion step every 5 time iterations \n"
            << " We are using P1 finite element in time, with the Simpson rule as quadrature rule "
            << std::endl << std::endl;
  int nx = 40;
  int ny = 80 ;

  Mesh2 Th(nx, ny, 0., 0., 1., 2.);
  double meshSize(1./nx);
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

  std::cout << " Mesh size \t" << meshSize << std::endl;
  std::cout << " Number of node \t" << Th.nv << std::endl;
  std::cout << " Number of element \t" << Th.nt << std::endl;

  std::cout << " Time Step \t" << dT << std::endl;
  std::cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  std::cout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  std::cout << " number of quadrature points in time : \t" << nbTime << std::endl;
  std::cout << " number of dof in time per time slab : \t" << ndfTime << std::endl;



  // Create files to save results
  std::string meshstr = to_string(nx)+to_string(ny);
  std::string pathOutpuFolder = "../../outputFiles/RBsimpson/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/RBsimpson/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  std::cout << " path to the output files : "+pathOutpuFolder << std::endl;
  std::cout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;
  std::cout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;


  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveStokesVTK = true;
  const bool curvatureVTKFile = true;
  const bool saveLevelSetVTK = false;
  const int frequencyPlottingTheSolution = 1;

  // Space for the velocity field
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange2 FEvelocity(2);
  FESpace2 VelVh(Th, FEvelocity);

  // The initial velocity is 0, so no function needed
  vector<Fun_h> vel(nbTime);
  for(int i=0;i<nbTime;++i) vel[i].init(VelVh);


  // levelSet stuff
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
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
  const int frequencyReinitialization = 2 ;


  // Declaration of the interface
  TimeInterface2 interface(nbTime);

  // FE for curvature problem
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " The curvature is computed using second order finite elements " << std::endl;
  Lagrange2 FEcurvature(2);

  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh2> stokes(qTime);
  CutFEM_Parameter mu("mu",10.,1.);
  CutFEM_Parameter rho("rho",1000.,100.);
  CutFEM_Parameter invmu("invmu",1./10,1.);


  const R sigma = 24.5;
  const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
  const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
  const CutFEM_Parameter& h(Parameter::h);

  list<int> dirichlet = {1,3};
  list<int> neumann = {2,4};

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

    // use the last position of the interface during previous iteration
    // as first step in the current iteration
    ls_k.begin()->swap(ls_k[lastQuadTime]);

    // move the levelSet
    for(int i=0;i<nbTime;++i) {
      projection(ls_k[i], ls[i]);

      // here we want to reinitializa the levelSet.
      // !!! WE CANNOT DO ON THE FIRST LEVELSET BECAUSE IT HAS TO MATCH THE LAST STEP
      if(iter%frequencyReinitialization == 0 && i == 1 && iter > 0) {
        std::cout << " ------- reinitialization LevelSet -------" << std::endl;
        reinitialization.perform(ls_k[i], ls[i]);
      }

      if(i<nbTime-1) {
        LevelSet levelSet(ls_k[i], vel[lastQuadTime], vel[lastQuadTime], dt_levelSet);
        ls_k[i+1].init(levelSet.rhs);
      }
    }
    interface.init(Th, ls);
    CutFESpace2 Wh(Vh, interface, {1,-1});
/*
                PROBLEM DEFINITION
*/
    Normal n;
    FunTest du(Wh,d), dp(Wh,1,d), v(Wh,d), q(Wh,1,d), dp1(Wh,1,d,0);
    FunTest Eun = (Eps(du)*n);
    FunTest D2nu = grad(grad(du)*n)*n, D2nv = grad(grad(v)*n)*n;

    stokes.initSpace(Wh, In);
    std::cout << " DOF       \t : \t" << stokes.nDoF << std::endl;

    // interpolate rhs
    Fun_h fh(Vh, fun_rhs);
    Rn datau0;
    stokes.initialSolution(datau0);

    // set uh and initial solution for N-S
    Rn data_uh(datau0);

    Fun_h u0(Wh, datau0);
    Fun_h uh(Wh, In, data_uh);

    // loop for Newton
    int iterNewton  = 0;
    while(1) {

      std::cout << " Begin assembly " << std::endl;
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
        stokes.addBilinearFormBorder(
          innerProduct(lambdaB*du.t()*n,v.t()*n)
          + innerProduct(dp, v.t()*n)
          - innerProduct(du.t()*n, q)
          - innerProduct(2.*mu*Eun.t()*n, v.t()*n)
          - innerProduct(du.t()*n, 2.*mu*Eun.t()*n)
          , In
          , neumann
        );

        R h3 = pow(meshSize,3);
        stokes.addFaceStabilization(
          innerProduct(1e-2*h*jump(grad(du)*n), rho*mu*jump(grad(v)*n))
          + innerProduct(1e-2*h3*jump(D2nu), rho*mu*jump(D2nv))
          ,In
        );
        stokes.addFaceStabilization(
          innerProduct(1e-2*h3*jump(grad(dp).t()*n), invmu*jump(grad(q).t()*n))
          ,In
        );
        // // Mass matrix with time term
        // stokes.addFaceStabilization(
        //   innerProduct(1e-2*cch*jump(grad(dt(du))*n), rho*mu*jump(grad(v)*n))
        //   ,In
        // );
        // // Mass matrix
        // stokes.addFaceStabilization(
        //   innerProduct(1e-2*cch*jump(grad(du)*n), rho*jump(grad(v)*n))
        // );

        // impose initial condition
        stokes.addBilinear(
          innerProduct(rho*du,v)
        );

      }
      std::cout << " Time Full assembly " << CPUtime() - tt0 << std::endl;


      /*
      CONSTRUCTION LINEAR PART
      */
      tt0 = CPUtime();
      stokes.addLinear(
        -innerProduct(fh.expression(2),rho*v)
        , In
      );

      //  COMPUTATION OF THE CURVATURE TERM
      for(int i=0;i<nbTime;++i) {

        Mesh2 cutTh(*interface[i]);
        FESpace2 cutVh(cutTh, FEcurvature);
        Mapping2 mapping(cutVh, ls_k[i]);

        Curvature curvature(cutVh, *interface[i], mapping);
        Fun_h meanCurvature(cutVh, curvature.rhs);

        ExpressionFunFEM<Mesh2> Kx(meanCurvature,0,op_id);
        ExpressionFunFEM<Mesh2> Ky(meanCurvature,1,op_id);
        FunTest vx(Wh,1,0), vy(Wh,1,1);

        stokes.addLinear(
           - innerProduct(Kx, average2(vx))*sigma
           - innerProduct(Ky, average2(vy))*sigma
           , interface
           ,i
          , In
        );

        //plotting the curvsture
        if(iter%frequencyPlottingTheSolution == 0 && curvatureVTKFile && i == 0){
          Fun_h solMC(cutVh, curvature.rhs);
          Paraview2 writerLS(cutVh, pathOutpuFigure+"curvature_"+to_string(iterfig)+".vtk");
          writerLS.add(solMC, "meanCurvature", 0, 2);
          writerLS.add(ls[i], "levelSet", 0, 1);
        }
      }
      // impose initial condition
      stokes.addLinear(
        -innerProduct(u0.expression(2),rho*v)
      );
      std::cout << " Time assembly rhs " << CPUtime() - tt0 << std::endl;


      tt0 = CPUtime();
      stokes.addMatMul(data_uh);
      std::cout << " Time RHS assembly A*u0 " << CPUtime() - tt0 << std::endl;
      // to optimize the Newton iteration.
      stokes.saveMatrix();


      /*
      CONSTRUCTION NON LINEAR PART
      */
      tt0 = CPUtime();
      FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      ExpressionFunFEM<Mesh2> u1(uh,0,op_id,op_id), u2(uh,1,op_id,op_id);

      stokes.addBilinear(
          innerProduct(du1*dx(u1) + du2*dy(u1) , rho*v1)
        + innerProduct(du1*dx(u2) + du2*dy(u2) , rho*v2)
        + innerProduct(u1*dx(du1) + u2 *dy(du1), rho*v1)
        + innerProduct(u1*dx(du2) + u2 *dy(du2), rho*v2)
        ,In
      );

      stokes.addLinear(
          innerProduct(u1*dx(u1) + u2*dy(u1), rho*v1)
        + innerProduct(u1*dx(u2) + u2*dy(u2), rho*v2)
        , In
      );

      std::cout << " Time Full assembly Non linear part " << CPUtime() - tt0 << std::endl;


      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , In
      );
      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,dp1), 0.
      //   , 0
      //   , In
      // );
      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , 0
      );

      // matlab::Export(stokes.DF, "DF"+to_string(iterNewton)+".dat");
      // matlab::Export(stokes.rhs, "rhs"+to_string(iterNewton)+".dat");

      stokes.solve();

      KN_<double> dw(stokes.rhs(SubArray(stokes.nDoF, 0)));
      double dist  = dw.l1()/dw.size();
      std::cout << " dw \t" << dist << std::endl;

      data_uh -= dw;

      stokes.rhs.resize(stokes.nDoF); stokes.rhs = 0.0;

      iterNewton+=1;


      if(iterNewton == 10 || dist < 1e-10) {
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

      double riseVelocity = integral(u0, 1, 1) / areaBubble;

      std::cout << "\n Features of the drop " << std::endl;
      std::cout << " Time              ->    " << GTime::current_time()+dT  << std::endl;
      std::cout << " Center Of Mass    ->    " << centerOfMass << std::endl;
      std::cout << " Circularity       ->    " << circularity << std::endl;
      std::cout << " Rise velocity     ->    " << riseVelocity << std::endl;

      outputData << GTime::current_time() << "\t"
                 << centerOfMass << "\t"
                 << circularity << "\t"
                 << riseVelocity << "\t"
                 << areaBubble <<  std::endl;
    }


    // std::cout << " Set Velocity " << std::endl;
    for(int i=0;i<nbTime;++i) {
      Fun_h sol(Wh,In,data_uh);
      set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    }


//  -----------------------------------------------------
//                     PLOTTING
//  -----------------------------------------------------
  if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles&& MPIcf::IamMaster()){
    std::cout << " Plotting " << std::endl;
    int ndf = 0;

    // Rn solv(Wh.NbDoF(), 0.);
    //
    // for(int i=0;i<Ih[0].NbDoF();++i){
    //   KN_<double> solvi(data_uh(SubArray(Wh.NbDoF(), ndf)));
    //   ndf += Wh.NbDoF();
    //   solv += solvi;
    // }
    if(saveStokesVTK) {
      KN_<double> solv(data_uh(SubArray(Wh.NbDoF(), 0)));
      Fun_h sol(Wh, solv);
      Paraview2 writerr(Wh, ls[0], pathOutpuFigure+"navierStokes_"+to_string(iterfig)+".vtk");
      writerr.add(sol, "velocity", 0, 2);
      writerr.add(sol, "pressure", 2, 1);
    }
    if(saveLevelSetVTK) {
      Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
      writerLS.add(ls[0], "levelSet", 0, 1);
    }
    iterfig++;
  }

    iter += 1;
  }
  std::cout << " -----------------------------------  " << std::endl;
  std::cout << " -----------------------------------  " << std::endl;
  std::cout << " END OF COMPUTATION " << std::endl;


  // Closing the opened files.
  outputData.close();
}
