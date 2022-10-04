/*

This is an implementation of example 2 in Paper "A Cut Finite Element Method for
two-phase flows with insoluble surfactants".
The example is a bubble in a shear flow.

*/

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "global.hpp"
#include "levelSet.hpp"
#include "curvature.hpp"
#include "projection.hpp"
#include "mapping.hpp"
#include "../num/matlab.hpp"
#include "../num/gnuplot.hpp"
// #include "curvature.hpp"
// #include "reinitialization.hpp"
// #include "../time_stuff.hpp"
// #include "projection.hpp"
// #include "../util/redirectOutput.hpp"



R fun_levelSet(const R2 P, const int i) { return sqrt(P.x*P.x + P.y*P.y) - 1.00001;}
R fun_rhs(const R2 P,const int i) { return 0;}
R fun_boundary(const R2 P,const int i) { return (i==0)?0.5*P.y : 0;}
R fun_init_surfactant(const R2 P, int i) { return 1.;}
// R sigma(const R w) {return 24.5 - 0.25*w;}
// R Dsigma(const R w) {return - 0.25*w;}

R fun_x(const R2 P,const int i) {return P[i];}
R fun_1(const R2 P,const int i) {return 1;}


// DEFINE SOME PARAMETER
class LambdaBoundary : public VirtualParameter {
  public :
  const CutFEMParameter& mu_;
  double G_, H_;
  LambdaBoundary(const CutFEMParameter& mu, double G, double H) : mu_(mu) , G_(G) , H_(H) {}
  double evaluate(int domain, double h, double meas, double measK, double measCut) const {
    double gamma = meas / h ;
    double alpha = measK / h / h;
    double val = mu_.evaluate(domain,h,meas,measK,measCut);
    return val/h*(G_+H_*gamma/alpha) ;
  }
};
class LambdaGamma : public VirtualParameter {
  public :
  const double mu1_;
  const double mu2_;
  LambdaGamma(const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
  double evaluate(int domain, double h, double meas, double measK, double measCut) const {
    double gamma = meas / h  ;
    double alphaK = measK / h / h;
    return  (0.5*mu1_ + 0.5*mu2_)*(100 + 10*gamma)/(alphaK);
  }
};
class WeightKappa : public VirtualParameter {
  public :
  const double mu1_;
  const double mu2_;
  WeightKappa (const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
  double mu(int i) const { return (i==0)? mu1_ : mu2_;}
  double evaluate(int dom0, double h, double meas, double measK, double measCut0) const {
    int dom1 = (dom0 == 0);
    double measCut1 = meas - measCut0;
    double alphaK0 = measCut0 / h / h;
    double alphaK1 = measCut1 / h / h;
    return  (mu(dom1)*alphaK0)/(mu(dom0)*alphaK1+mu(dom1)*alphaK0);
  }
};
class Weight2Kappa : public VirtualParameter {
  public :
  const double mu1_;
  const double mu2_;
  Weight2Kappa (const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
  double mu(int i) const { return (i==0)? mu1_ : mu2_;}
  double evaluate(int dom0, double h, double meas, double measK, double measCut0) const {
    int dom1 = (dom0 == 0);
    double measCut1 = meas - measCut0;
    double alphaK0 = measCut0 / h / h;
    double alphaK1 = measCut1 / h / h;
    return  (mu(dom0)*alphaK1)/(mu(dom0)*alphaK1+mu(dom1)*alphaK0);
  }
};


int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;

  MPIcf cfMPI(argc,argv);
  int thread_count_tmp = 1;
  cout << "Threads: ";
  cin >> thread_count_tmp;
  MPIcf::Bcast(thread_count_tmp, MPIcf::Master(), 1);
  const int thread_count = thread_count_tmp;
  omp_set_num_threads(thread_count);
  const double cpubegin = MPIcf::Wtime();

  // MESH DEFINITION
  // ---------------------------------------------
  int nx = 50;
  int ny = 20;
  Mesh Kh(nx, ny, -5., -2., 10., 4.);

  double mesh_size = (10./(nx-1));


  // Mesh Kh("../mesh/shearMesh7530_10040.msh");
  // double mesh_size = Kh.get_mesh_size();

  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Background mesh " << std::endl;
  Kh.info();

  //   //
  //   Mesh2 Th("../mesh/shearMesh7530_10040.msh");
  //   int nx = 100;
  //   int ny = 40;
  //   double meshSize = 10./nx/2;
  //   // Mesh2 Th(nx, ny, -5., -2., 10., 4.);
  //

  ProblemOption optionProblem;
  optionProblem.solver_name_ = "mumps";
  optionProblem.clear_matrix_ = true;
  std::vector<std::map<std::pair<int,int>, double>> mat_NL(thread_count);


  // TIME DEFINITION
  // ---------------------------------------------
  int division_mesh_size = 4;
  double final_time = 1;
  double dT = mesh_size/division_mesh_size;
  GTime::total_number_iteration = final_time/dT;
  dT = final_time / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);

  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nb_quad_time = qTime.n;
  const Uint ndf_time_slab = Ih[0].NbDoF();
  const Uint last_quad_time = nb_quad_time-1;


  // PROBLEM AND PARAMETER DEFINITION
  // ----------------------------------------------
  CutFEM<Mesh2> stokes(qTime,thread_count, optionProblem);
  CutFEMParameter mu(0.1,0.1);
  CutFEMParameter rho(1.,1.);
  CutFEMParameter invmu(10.,10.);
  LambdaBoundary lambdaB(mu, 10, 100);
  LambdaGamma lambdaG(mu(0), mu(1));
  WeightKappa  kappa(mu(0), mu(1));
  Weight2Kappa kappa2(mu(0), mu(1));

  // const MeshParameter& h(Parameter::h);
  double hh = mesh_size;
  const double sigma0 = 0.2;
  const double beta   = 0.;
  const double epsilon_surfactant = 0.1;


  // WRITE THE PROBLEM PARAMETERS
  // ----------------------------------------------
  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Discretization of the problem : \n"
  << " h = " << mesh_size
  << "\n t_begin = " << GTime::t0 << " and t_final = " << final_time
  << "\n dt = " << dT
  << "\n number of iteration = " << GTime::total_number_iteration
  << std::endl;
  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Parameters of the problem : \n"
  << " mu_1 = " << mu(0) << " and mu_2 = "<< mu(1) << " \n"
  << " rho_1 = "<< rho(0) << " and rho_2 = " << rho(1) << " \n"
  << " the surface tension sigma = " << sigma0 << " \n"
  << " beta = " << beta << " \n"
  << " D_Gamma = " << epsilon_surfactant
  << std::endl;
  std::cout << " ------------------------------------ " << std::endl;


  // SPACE DEFINITION
  // ---------------------------------------------
  Lagrange2 FEu(2);
  Space Vh(Kh, FEu);
  Space Lh  (Kh, DataFE<Mesh>::P1);
  Space Lh_k(Kh, DataFE<Mesh>::P2);
  Lagrange2 FEcurv(1);
  Space Vh1(Kh, FEcurv);

  //   // Create files to save results
  //   std::string meshstr = to_string(nx)+to_string(ny);
  //   std::string pathOutpuFolder = "../../outputFiles/shearBubbleConservative/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  //   std::string pathOutpuFigure = "../../outputFiles/shearBubbleConservative/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  //   CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  //   myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  //   myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  //   if(MPIcf::IamMaster()) {
  //     std::experimental::filesystem::create_directories(pathOutpuFigure);
  //   }
  //   std::ofstream outputData(pathOutpuFolder+"dataDrop.dat", std::ofstream::out);
  //   myCout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  //   myCout << " Creating the file \"output.txt\" to save cout" << std::endl;
  //
  //
  //   myCout << " SIMULATION OF A BUBBLE IN SHEAR FLOW WITH INSOLABLE SURFACTANT. \n"
  //             << " We consider the navier Stokes equations coupled to the insoluble surfactant problem.\n"
  //             << " We solve the non linear system using the Newton method \n"
  //             << " We approximate the surfactant concentration using P1/P1 space time finite elements  \n"
  //             << " We are using the levelSet P2 with a reinitializtion step every 5 time iterations \n"
  //             << " We are using P1 finite element in time, with the Simpson rule as quadrature rule "
  //             << std::endl << std::endl;
  //
  //
  //   myCout << " Mesh size \t" << meshSize << std::endl;
  //   myCout << " Number of node \t" << Th.nv << std::endl;
  //   myCout << " Number of element \t" << Th.nt << std::endl;
  //
  //   myCout << " Time Step \t" << dT << std::endl;
  //   myCout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  //   myCout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  //   myCout << " number of quadrature points in time : \t" << nbTime << std::endl;
  //   myCout << " number of dof in time per time slab : \t" << ndfTime << std::endl;
  //   return 0;
  //
  //
  //   // Set parameters for paraview PLOTTING
  //   const bool writeVTKFiles = true;
  //   const bool saveStokesVTK = true;
  //   const bool saveSurfactantVTK = true;
  //   const bool curvatureVTKFile = false;
  //   const bool saveLevelSetVTK = true;
  const int frequency_plotting = 1;

  // INTERPOLATION OF THE VELOCITY (BOUNDARY CONDITION)
  // ----------------------------------------------
  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Interpolation of velocity " << std::endl;
  vector<Fun_h> vel(nb_quad_time);
  for(int i=0;i<nb_quad_time;++i) vel[i].init(Vh, fun_boundary);

  // DECLARATION OF INTERFACE AND LEVELSET
  // ----------------------------------------------
  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Create interface and levelSet " << std::endl;
  TimeInterface<Mesh> interface(qTime);
  double dt_levelSet = dT/(nb_quad_time-1);
  vector<Fun_h> ls_k(nb_quad_time), ls(nb_quad_time);
  for(int i=0;i<nb_quad_time;++i) ls_k[i].init(Lh_k, fun_levelSet);
  for(int i=0;i<nb_quad_time;++i) ls[i].init(Lh, fun_levelSet);


  //   projection(ls_k[0], ls[nbTime-1]);
  //   levelSet.setStrongBC({2,4});

  //   CReinitialization<Mesh> reinitialization;
  //   reinitialization.number_iteration = 4;
  //   reinitialization.epsilon_diffusion = 1e-3;
  //   reinitialization.dt = dT/8;
  //   reinitialization.ON_OFF = "ON";
  //   reinitialization.ODE_method = "Euler";
  //   reinitialization.mass_correction = "ON";
  //   reinitialization.precision_correction = 1e-6;
  //   reinitialization.max_iteration_correction = 10;
  //   reinitialization.info();
  //   const int frequencyReinitialization = 5;



  //   list<int> dirichlet = {1,2,3,4,20,40};


  // INTERPOLATE RHS AND BOUNDARY CONDITION
  // ----------------------------------------------
  std::cout << " ------------------------------------" << std::endl;
  std::cout << " Interpolate boundary and rhs " << std::endl;
  Fun_h fh(Vh, fun_rhs);
  Fun_h gh(Vh, fun_boundary);


  std::cout << " ------------------------------------" << std::endl;
  std::cout << " ------------------------------------" << std::endl;

  int iter = 0, ifig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);
    if(iter == 0){
      std::cout << " ITERATION \t : \t" << iter+1 << " / " << GTime::total_number_iteration << std::endl;
    }
    else {
      std::cout << "\r";
      std::cout << " ITERATION \t : \t" << iter+1 << " / " << GTime::total_number_iteration << " \t Error : ************** " ;
      std::cout.flush();
      globalVariable::verbose = 0;
    }

    // MOVE THE LEVELSET AND CREATE INTERFACES
    // ----------------------------------------------
    {
      if(iter == 0){
        std::cout << " ------------------------------------" << std::endl;
        std::cout << " Move the levelSet and create interfaces " << std::endl;
      }
      double t0 = MPIcf::Wtime();
      globalVariable::verbose = 0;

      ls_k.begin()->swap(ls_k[last_quad_time]);
      // ls.begin()->swap(ls[last_quad_time]);

      for(int i=0;i<nb_quad_time;++i) {

        projection(ls_k[i], ls[i]);

        // here we want to reinitializa the levelSet.
        // !!! WE CANNOT DO ON THE FIRST LEVELSET BECAUSE IT HAS TO MATCH THE LAST STEP
        // if(iter%frequencyReinitialization == 0 && i == 1 && iter > 0) {
        //   std::cout << " ------- reinitialization LevelSet -------" << std::endl;
        //   reinitialization.perform(ls_k[i], ls[i]);
        // }


        interface.init(i,Kh,ls[i]);
        if(i<nb_quad_time-1) {
          LevelSet::move(ls_k[i], vel[last_quad_time], vel[last_quad_time], dt_levelSet, ls_k[i+1]);
        }
      }
      if(iter == 0) {
        globalVariable::verbose = 1;
        std::cout << " Time moving interface : \t" << MPIcf::Wtime() - t0 << std::endl;
      }
    }


    // CREATE ACTIVE MESH AND CUTSPACES
    // ----------------------------------------------
    double t0_cmesh = MPIcf::Wtime();
    CutMesh Kh_i(Kh, interface);
    CutSpace Wh(Kh_i, Vh);
    CutSpace Ph(Kh_i, Lh);
    CutMesh Kh_0(Kh);
    Kh_0.createSurfaceMesh(interface);
    CutSpace Sh(Kh_0, Lh);
    if(iter == 0) {
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Cut Mesh " << std::endl;
      Kh_i.info();
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Cut Space for velocity " << std::endl;
      Wh.info();
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Cut Space for pressure " << std::endl;
      Ph.info();
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Surface Space for surfactant " << std::endl;
      Sh.info();
      std::cout << " TIME BUILDING CUT MESH/SPACE : \t" << MPIcf::Wtime() - t0_cmesh << std::endl;
    }

    // DEFINE TEST FUNCTIONS
    // ----------------------------------------------
    Normal n;
    FunTest du(Wh,2), dp(Ph,1), v(Wh,2), q(Ph,1), dp1(Wh,1,0,0);
    FunTest Eun = (Eps(du)*n);
    FunTest D2nu = grad(grad(du)*n)*n, D2nv = grad(grad(v)*n)*n;
    FunTest ds(Sh,1), r(Sh,1);

    // DEFINE THE PROBLEM ON THE CORRECT SPACES
    // ----------------------------------------------
    stokes.initSpace(Wh, In);
    stokes.add(Ph, In);
    stokes.add(Sh, In);

    if(iter == 0) {
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Problem's DOF : \t" << stokes.get_nb_dof() << std::endl;
    }

    // INITIALIZE VECTORS USED
    // -------------------------------------
    double t0_init_sol = MPIcf::Wtime();
    if(iter == 0) {
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " Initialize vector and solution " << std::endl;
    }
    Rn data_init;
    stokes.initialSolution(data_init);
    Rn data_all(data_init);

    int idxp0 = Wh.NbDoF()*In.NbDoF();
    int idxs0 = idxp0 + Ph.NbDoF()*In.NbDoF();
    KN_<double> data_uh0(data_init(SubArray(Wh.NbDoF(),0)));
    KN_<double> data_sh0(data_init(SubArray(Sh.NbDoF(),idxs0)));
    KN_<double> data_uh(data_all(SubArray(Wh.NbDoF()*In.NbDoF(),0)));
    KN_<double> data_ph(data_all(SubArray(Ph.NbDoF()*In.NbDoF(),idxp0)));
    KN_<double> data_sh(data_all(SubArray(Sh.NbDoF()*In.NbDoF(),idxs0)));

    Fun_h uh(Wh, In, data_uh);
    Fun_h ph(Ph, In, data_ph);
    Fun_h sh(Sh, In, data_sh);
    Fun_h u0(Wh, data_uh0);
    Fun_h s0(Sh, data_sh0);


    if(iter == 0) {
      interpolate(Wh, data_uh0, fun_boundary);
      interpolate(Sh, data_sh0, fun_init_surfactant);
      std::cout << " Time initialization : \t" << MPIcf::Wtime() - t0_init_sol << std::endl;
    }

    // NEWTON ITERATION
    // ----------------------------------------------
    int iterNewton  = 0;
    double t0_newton = MPIcf::Wtime();
    if(iter == 0) {
      std::cout << " ------------------------------------" << std::endl;
      std::cout << " START NEWTON ITERATION " << std::endl;
    }
    while(1) {
      double tt0 = MPIcf::Wtime();
      if(iterNewton == 0){
        stokes.addBilinear(
          innerProduct(dt(du), rho*v)
          + contractProduct(2*mu*Eps(du),Eps(v))
          - innerProduct(dp, div(v)) + innerProduct(div(du), q)
          , Kh_i
          , In
        );
        stokes.addBilinear(
          innerProduct(jump(du), -2*mu*average(Eps(v)*n, kappa))
          + innerProduct(-2*mu*average(Eps(du)*n, kappa), jump(v))
          + innerProduct(lambdaG*jump(du), jump(v))
          + innerProduct(average(dp, kappa), jump(v.t()*n))
          - innerProduct(jump(du.t()*n), average(q, kappa))
          , interface
          , In
        );
        stokes.addBilinear(
          innerProduct(lambdaB*du,v)
          + innerProduct(dp, v.t()*n)
          - innerProduct(du.t()*n, q)
          - innerProduct(2.*mu*Eps(du)*n, v)
          - innerProduct(du, 2.*mu*Eps(v)*n)
          , Kh_i
          , INTEGRAL_BOUNDARY
          , In
          // , dirichlet
        );


        stokes.addBilinear(
          - innerProduct(ds, dt(r))
          + innerProduct(epsilon_surfactant*gradS(ds), gradS(r))
          , interface
          , In
        );



        double h3 = pow(mesh_size,3);
        stokes.addFaceStabilization(
          innerProduct(1e-1*hh*jump(grad(du)*n), mu*jump(grad(v)*n))
          + innerProduct(1e-1*h3*jump(D2nu)      , mu*jump(D2nv))
          + innerProduct(1e-1*h3*jump(grad(dp)*n), invmu*jump(grad(q)*n))
          , Kh_i
          , In
        );

        stokes.addFaceStabilization(
          innerProduct(1e0*jump(grad(ds)*n), jump(grad(r)*n))
          , Kh_0
          , In
        );
        //         // if stab interface, can have h in both stabilization
        // stokes.addBilinear(
        //   innerProduct(1e-2*hh*grad(ds)*n, grad(r)*n)
        //   , interface
        //   , In
        // );

        // impose initial condition
        stokes.addBilinear(
          innerProduct(rho*du,v)
          , Kh_i
        );
        stokes.addBilinear(
          innerProduct(ds, r)
          , *interface(last_quad_time)
          , In
          , last_quad_time
        );


      }


      for(int i=0;i<nb_quad_time;++i) {      // computation of the curvature
        if(i >0 || iterNewton > 0) globalVariable::verbose = 0;

        CutMesh cutTh(Kh);
        cutTh.createSurfaceMesh(*interface[i]);
        CutSpace cutVh(cutTh, Vh);

        Mapping2 mapping(cutVh, ls_k[i], ls[i]);
        Curvature<Mesh> curvature(cutVh, interface[i]);

        auto tq = qTime(i);
        double tid = (double)In.map(tq);
        ExpressionLinearSurfaceTension<Mesh> sigmaW(sh, sigma0, beta, tid);

        Rn data_H(cutVh.get_nb_dof());
        data_H = curvature.solve(sigmaW, mapping);
        // data_H = curvature.solve(mapping);

        Fun_h H(cutVh, data_H);

        ExpressionFunFEM<Mesh2> Hx(H,0,op_id);
        ExpressionFunFEM<Mesh2> Hy(H,1,op_id);
        FunTest vx(Wh,1,0), vy(Wh,1,1);

        stokes.addLinear(
          - innerProduct(Hx , average(vx, kappa2))
          - innerProduct(Hy , average(vy, kappa2))
          , interface
          , In , i
        );
        if(i == 0){
          Paraview<Mesh> writerS(cutTh, "shearExampleCurvature_"+to_string(ifig)+".vtk");
          writerS.add(ls[i], "levelSet", 0, 1);
          writerS.add(H,"meanCurvature", 0, 2);
        }
        // return 0;
      }
      if(iter == 0 && iterNewton == 0) globalVariable::verbose = 1;

      if(iter == 0 && iterNewton == 0) {
        std::cout << " TIME ASSEMBLY MATRIX \t" << MPIcf::Wtime() - tt0 << std::endl;
      }
      /*
      CONSTRUCTION LINEAR PART
      */
      tt0 = MPIcf::Wtime();
      stokes.addLinear(
        -innerProduct(fh.expression(2),rho*v)
        , Kh_i
        , In
      );
      stokes.addLinear(
        -innerProduct(u0.expression(2),rho*v)
        , Kh_i
      );
      stokes.addLinear (
        -innerProduct(s0.expression(), r)
        , *interface(0)
        , In
        , 0
      );

      stokes.addLinear(
        - innerProduct(gh.expression(2),lambdaB*v)
        + innerProduct(gh.expression(2),2.*mu*Eps(v)*n)
        + innerProduct(gh.expression(2), q*n)
        , Kh_i
        , INTEGRAL_BOUNDARY
        , In
      );
      if(iter == 0 && iterNewton == 0) {
        std::cout << " TIME ASSEMBLY RHS \t" << MPIcf::Wtime() - tt0 << std::endl;
      }

      stokes.gather_map();
      stokes.addMatMul(data_all);
      tt0 = MPIcf::Wtime();
      if(iter == 0 && iterNewton == 0) {
        std::cout << " TIME MATRIX MULTIPLICATION \t" << MPIcf::Wtime() - tt0 << std::endl;
      }
      tt0 = MPIcf::Wtime();


      /*
      CONSTRUCTION NON LINEAR PART
      */
      stokes.set_map(mat_NL);
      mat_NL[0] = stokes.mat_[0];

      FunTest du1(Wh,1,0), du2(Wh,1,1), v1(Wh,1,0), v2(Wh,1,1);
      ExpressionFunFEM<Mesh2> ux(uh,0,op_id,op_id), uy(uh,1,op_id,op_id);
      ExpressionFunFEM<Mesh2> s(sh,0,op_id,op_id);

      stokes.addBilinear(
        innerProduct(du1*dx(ux) + du2*dy(ux), rho*v1)
        + innerProduct(du1*dx(uy) + du2*dy(uy), rho*v2)
        + innerProduct(ux*dx(du1) + uy*dy(du1), rho*v1)
        + innerProduct(ux*dx(du2) + uy*dy(du2), rho*v2)
        , Kh_i
        ,In
      );
      stokes.addLinear(
        innerProduct(ux*dx(ux) + uy*dy(ux), rho*v1)
        + innerProduct(ux*dx(uy) + uy*dy(uy), rho*v2)
        , Kh_i
        , In
      );

      stokes.addBilinear(
        - innerProduct(ds*average(ux) , dx(r))
        - innerProduct(ds*average(uy) , dy(r))
        - innerProduct(s*average(du,0.5,0.5) , grad(r))
        , interface
        , In
      );
      stokes.addLinear(
        - innerProduct(s, dx(r)*average(ux))
        - innerProduct(s, dy(r)*average(uy))
        , interface
        , In
      );
      stokes.addLagrangeMultiplier(
        innerProduct(1.,dp1), 0.
        , Kh_i
        , In
      );
      // stokes.addLagrangeMultiplier(
      //   innerProduct(1.,dp1), 0.
      //   , Kh_i
      //   // , 0
      //   // , In
      // );

      if(iter == 0 && iterNewton == 0) {
        std::cout << " ASSEMBLY NON LINEAR \t" << MPIcf::Wtime() - tt0 << std::endl;
      }
      tt0 = MPIcf::Wtime();

      // gather(mat_NL);
      // matlab::Export(mat_NL[0], "mat_1.dat");
      // return 0;
      stokes.solve(mat_NL, stokes.rhs_);



      if(iter == 0 && iterNewton == 0) {
        std::cout << " TIME SOLVER \t" << MPIcf::Wtime() - tt0 << std::endl;
      }

      KN_<double> dwu(stokes.rhs_(SubArray(Wh.NbDoF()*In.NbDoF(), 0)));
      KN_<double> dwp(stokes.rhs_(SubArray(Ph.NbDoF()*In.NbDoF(), idxp0)));
      KN_<double> dws(stokes.rhs_(SubArray(Sh.NbDoF()*In.NbDoF(), idxs0)));
      double dist  = dwu.l1()/dwu.size();
      double dists = dws.l1()/dws.size();
      if(iter == 0){
        std::cout << " Iterartion error || du ||_infty = " << dist << std::endl;
        std::cout << " Iterartion error || ds ||_infty = " << dists << std::endl;
      }
      else {
        std::cout << "\r";
        std::cout << " ITERATION \t : \t" << iter+1 << " / " << GTime::total_number_iteration << " \t Error : " << max(dist, dists);
        std::cout.flush();
      }
      KN_<double> dw(stokes.rhs_(SubArray(stokes.get_nb_dof(), 0)));
      data_all -= dw;

      stokes.rhs_.resize(stokes.get_nb_dof()); stokes.rhs_ = 0.0;

      iterNewton+=1;

      if(iterNewton == 10 || max(dist, dists) < 1e-10) {
        stokes.saveSolution(data_all);
        stokes.cleanBuildInMatrix();
        stokes.set_map();

        if(iter == 0) {
          std::cout << " TIME NEWTON ITERATION : \t" << MPIcf::Wtime()-t0_newton << std::endl;
        }

        break;
      }
    }
    // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
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
    //   << centerOfMass << "\t"
    //   << circularity << "\t"
    //   << riseVelocity << "\t"
    //   << areaBubble <<  "\t"
    //   << qend << "\t"
    //   << fabs(q - qend) << "\t"
    //   << fabs(qend - initialConcentrationSurfactant) << std::endl;
    //
    // }
    //

    //     // std::cout << " Set Velocity " << std::endl;
    //     for(int i=0;i<nbTime;++i) {
    //       Fun_h sol(Wh,In,data_uh);
    //       set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
    //     }
    //
    //
    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(MPIcf::IamMaster() && (iter%frequency_plotting == 0 || iter+1 == GTime::total_number_iteration)) {

      {
        std::string filename = "shearExample_"+to_string(ifig)+".vtk";
        if(iter == 0){ std::cout << " Plotting -> " << filename << std::endl;}
        Paraview<Mesh> writer(Kh_i, filename);
        writer.add(ls[0],"levelSet", 0, 1);
        writer.add(uh,"velocity", 0, 2);
        writer.add(ph,"pressure", 0, 1);
      }
      {
        std::string filename = "shearExampleSurfactant_"+to_string(ifig)+".vtk";
        if(iter == 0){ std::cout << " Plotting -> " << filename << std::endl;}
        Paraview<Mesh> writer(Kh_0, filename);
        writer.add(ls[0],"levelSet_0", 0, 1);
        writer.add(ls[2],"levelSet_2", 0, 1);
        writer.add(sh,"surfactant_0", 0, 1);
      }
      //     if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
      //       std::cout << " Plotting " << std::endl;
      //
      //       if(saveStokesVTK) {
      // KN_<double> solu(data_uh(SubArray(Wh.NbDoF(), 0)));
      //         KN_<double> solp(data_ph(SubArray(Ph.NbDoF(), 0)));

      // Fun_h uh(Wh, solv);
      // Paraview2 writerr(Wh, ls[0], pathOutpuFigure+"navierStokes_"+to_string(iterfig)+".vtk");

      // }
      // if(saveSurfactantVTK) {
      //
      //         // Rn sols(cutWh.NbDoF(), 0.);
      //         //
      //         // for(int i=0;i<Ih[0].NbDoF();++i){
      //         //   int ndf = n0;
      //         //   KN_<double> solsi(uh(SubArray(cutWh.NbDoF(), ndf)));
      //         //   ndf += cutWh.NbDoF();
      //         //   sols += solsi;
      //         // }
      //         // Fun2 sol(cutWh, sols);
      //
      //         KN_<double> sols(data_uh(SubArray(cutWh.NbDoF(), n0)));
      //         Fun_h sol(cutWh, sols);
      //
      //         Paraview2 writer(cutWh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
      //         writer.add(sol, "surfactant", 0, 1);
      //         writer.add(ls[0], "levelSet", 0, 1);
      //         // writer.add(ls[2], "levelSet1", 0, 1);
      //
      //       }
      //       if(saveLevelSetVTK) {
      //         Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
      //         writerLS.add(ls[0],"levelSet", 0, 1);
      //       }
      //
      ifig++;
    }

    iter += 1;
    globalVariable::verbose = 0;
  }
  //   myCout << " -----------------------------------  " << std::endl;
  //   myCout << " -----------------------------------  " << std::endl;
  //   myCout << " END OF COMPUTATION " << std::endl;
  std::cout << "\n\n ------------------------------------" << std::endl;
  std::cout << " -----------------------------------  " << std::endl;
  std::cout << " Computation time : \t " << MPIcf::Wtime() - cpubegin << std::endl;
  //
  //   // Closing the opened files.
  //   outputData.close();
}
