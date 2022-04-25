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
#include "../num/matlab.hpp"
#include "paraview.hpp"

namespace NumericSurfactant2D0 {
  R fun_levelSet(const R2 P, const int i, const double t ) {
    R x = P.x,  y = P.y;
    return sqrt((x-t)*(x-t) + y*y) - 1;
    // return sqrt((x)*(x) + y*y) - 1;
  }
  R fun_init_surfactant(const R2 P, int i) {    // = 8PI
    return P.y/1+2;
  }
  R fun_velocity(const R2 P, int i,const double t) {
    return (i==0)?1:0;
  }
  // R fun_sol_surfactant(const R2 P,  const int i, const R t) {
  //   return P.y/1+2;;
  // }
  R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-0.25*t);}

  R fun_rhs(const R2 P, const int cc, const R t) {
    R x = P.x,  y = P.y;
    R r = y*exp(-t/4) - (x*y*exp(-t/4))/4 - (y*exp(-t/4)*(3*t - 4*x))/(t*t - 2*t*x + x*x + y*y);
    return r;
  }
  // R fun0(const R2 P, const R t) {return P.x*P.y*exp(-4*t);}
  // R fun_sol_surfactant(const R2 P, const R t) {
  //   R x = P.x,  y = P.y;
  //   return exp(-t/4)*y/(sqrt((x-t)*(x-t)+y*y)) + 2;
  // }
  // R2 fun_velocity_field(const R2 P, int i) {
  //   //return R2((P.y+2)*(P.y+2)/3,0);
  //   return (i==0);
  // }
  // R fun_div_surfactant(const R2 P, const R t) {
  //   return 0;
  // }
  // R fun_rhs_surfactant(const R2 P, const R t) {
  //   return 0;
  // }

  // R fun_levelSet(const R2 P) {
  //   R x = P.x,  y = P.y;
  //   return sqrt(x*x + y*y) - 2;
  // }
  // Diff<R,2> fun_levelSet_timeDA(Diff<R,2> P[2], const R t) {
  //   return sqrt((P[0]-t)*(P[0]-t) + P[1]*P[1]) - 2;
  // }

}
namespace NumericSurfactantEllipse2D {
  R fun_init_surfactant(const R2 P, const int i) {return P.x*P.y ;} //2*... = 4*pi
  R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-t/4);}
  R fun_velocity(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    R b = M_PI/4*(cos(2*M_PI*t))/a;
    return (i==0) ? b*P.x : 0;
  }
  R fun_levelSet(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    // a = 1;
    return sqrt(P.x*P.x/a + P.y*P.y) - 1.0000079;
  }
  R fun_rhs(const R2 P, const int cc, const R t) {
    R x = P.x,  y = P.y;
    // R r = (4*x*y*exp(-4*t)*(sin(2*M_PI*t) + 4)*
    //        (2048*x*x*y*y + 336*y*y*y*y*pow(sin(2*M_PI*t),2) +
    //         52*y*y*y*y*pow(sin(2*M_PI*t),3) + 3*y*y*y*y*pow(sin(2*M_PI*t),4)
    //         + 64*x*x*x*x*sin(2*M_PI*t) + 960*y*y*y*y*sin(2*M_PI*t) + 1024*x*x*x*x
    //         + 1024*y*y*y*y + 144*x*x*y*y*pow(sin(2*M_PI*t),2)
    //         + 4*x*x*y*y*pow(sin(2*M_PI*t),3) + 1024*x*x*y*y*sin(2*M_PI*t)))
    //   /pow((y*y*pow(sin(2*t*M_PI),2) + 8*y*y*sin(2*t*M_PI) + 16*x*x + 16*y*y),3)
    //   - 4*x*y*exp(-4*t) + (x*y*M_PI*exp(-4*t)*cos(2*M_PI*t))/(sin(2*M_PI*t) + 4)
    //   + (x*y*y*y*M_PI*exp(-4*t)*cos(2*M_PI*t)*(sin(2*M_PI*t) + 4))/
    //   (y*y*pow(sin(2*M_PI*t),2) + 8*y*y*sin(2*M_PI*t) + 16*x*x + 16*y*y);
    R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t)
    + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi)
    + 16*x*x + 16*y*y, 2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t))/(sin(2*pi*t) + 4)
    + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
    return r;
  }
  // R fun_rhs(const R2 P, const int cc, const R t) {
  //   R x = P.x,  y = P.y;
  //   R r = (4*x*y*exp(-t/4)*(sin(2*pi*t) + 4)*(3*y*y*sin(2*pi*t)*sin(2*pi*t) + 4*x*x*sin(2*pi*t) + 28*y*y*sin(2*pi*t) + 64*x*x + 64*y*y))/pow(y*y*sin(2*t*pi)*sin(2*t*pi) + 8*y*y*sin(2*t*pi) + 16*x*x + 16*y*y,2) - (x*y*exp(-t/4))/4 + (x*y*pi*exp(-t/4)*cos(2*pi*t)*(y*y*sin(2*pi*t) - 4*x*x + 4*y*y))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y) + (x*y*y*y*pi*exp(-t/4)*cos(2*pi*t)*(sin(2*pi*t) + 4))/(y*y*sin(2*pi*t)*sin(2*pi*t) + 8*y*y*sin(2*pi*t) + 16*x*x + 16*y*y);
  //   return r;
  // }
}

#define FORMULATION1;
// #define FORMULATION2;

using namespace NumericSurfactantEllipse2D ;
// using namespace NumericSurfactant2D0 ;
typedef Mesh2 Mesh;
typedef FESpace2 Space;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;

double solveCG_1(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
  const int d = 2;

  const double cpubegin = CPUtime();

  double errL2 = 0.;

  // BACKGROUND MESH
  // Mesh Kh(nx, nx, -2, -2, 4., 4.);
  Space Vh(Kh, DataFE<Mesh2>::P1);
  double meshSize = hi;//(4./(nx-1));

  // TIME PARAMETER
  int divisionMeshsize = 4;
  // double dT = dt;//meshSize/divisionMeshsize;
  GTime::total_number_iteration = int(tfinal/dT);
  dT = tfinal / GTime::total_number_iteration;
  GTime::time_step = dT;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  //CREATE SPACE AND VELOCITY FUNCTION
  Lagrange2 FEvelocity(2);
  Space VelVh  (Kh, FEvelocity);
  vector<Fun_h> vel(nbTime);


  // CREATE THE LEVELSET
  Space Lh  (Kh, DataFE<Mesh2>::P1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls(nbTime);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
  // LevelSet2 levelSet(Lh_k);

  // TIME INTERFACE AT ALL QUADRATURE TIME
  TimeInterface<Mesh> interface(qTime);

  // INITIALIZE THE PROBLEM
  CutFEM<Mesh2> surfactant(qTime);
  const R epsilon_S = 1.;

  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
    double tid = GTime::current_time();
    // INITIALIZE LEVELSET, VELOCITY, INTERFACE
    ls.begin()->swap(ls[nbTime-1]);
    for(int i=0;i<nbTime;++i) {
      R tt = In.Pt(R1(qTime(i).x));
      ls[i].init(Lh, fun_levelSet, tt);
      vel[i].init(VelVh, fun_velocity, tt);
      interface.init(i,Kh,ls[i]);
    }

    // CREATE THE ACTIVE MESH AND CUT SPACE
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.createSurfaceMesh(interface);
    CutSpace Wh(Kh0, Vh);

    // INITIALIZE THE PROBLEM
    Rn datau0;
    surfactant.initSpace(Wh, In);
    surfactant.initialSolution(datau0);
    KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
    if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
    Rn uh(datau0);
    Fun_h u0(Wh, datau0);

    Normal n;
    FunTest s(Wh,1), r(Wh,1);


    std::cout << " Begin assembly " << std::endl;
    double tt0 = MPIcf::Wtime();

    surfactant.addBilinear(
      innerProduct(dt(s), r)
      + innerProduct(epsilon_S*gradS(s), gradS(r))
      , interface
      , In
    );

    for(int i=0;i<nbTime;++i) {      // computation of the curvature
      ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
      surfactant.addBilinear(
        innerProduct(dx(s)*vx + dy(s)*vy, r)
        + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
        , interface
        ,In, i
      );
    }
    surfactant.addBilinear(
      innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
      , Kh0
      , innerFacet
      , In
    );

    //     // surfactant.addBilinearFormInterface(
    //     //   innerProduct(1e-2*grad(s)*n, grad(r)*n)
    //     //   ,In
    //     // );
    //
        surfactant.addBilinear(
          innerProduct(s, r)
          , *interface(0)
          , In
          , 0
        );
        surfactant.addLinear(
          innerProduct(u0.expression(), r)
          , *interface(0)
          , In
          , 0
        );

        Fun_h funrhs(Wh, In, fun_rhs);
        // Fun_h funrhsp(cutVh, In);
        // projection(funrhs, funrhsp, In, interface ,0);
        surfactant.addLinear(
          innerProduct(funrhs.expression(), r)
          , interface
          , In
        );
        // double  intF = integralSurf(funrhsp, In, interface) ;
    //     intFnp = integralSurf(funrhs, In, qTime);
    //
    //
    // matlab::Export(surfactant.mat_, "matA.dat");
    // surfactant.cleanMatrix();
    surfactant.solve();

    KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
    uh = dw;
    surfactant.saveSolution(uh);
    //



    //  -----------------------------------------------------
    //                COMPUTE ERROR
    //  -----------------------------------------------------
    {
      Rn sol(Wh.get_nb_dof(), 0.);
      sol += uh(SubArray(Wh.get_nb_dof(), 0));
      Fun_h funuh(Wh, sol);
      errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(0), tid,0,1);
      std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

      sol  += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
      errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(nbTime-1), tid+dT,0,1);
      std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
    }



    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    // if(MPIcf::IamMaster() ) {
    //     // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
    //     Fun_h sol(Wh, uh);
    //
    //     // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
    //     Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");
    //
    //     writer.add(sol , "surfactant", 0, 1);
    //     writer.add(ls[0], "levelSet",  0, 1);
    //     // writer.add(ls[2], "levelSet2", 0, 1);
    //   }


    // return 0.;
    iter += 1;
  }

  return errL2;
}
double solveCG_2(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
  const int d = 2;
  // typedef Mesh2 Mesh;
  // typedef FESpace2 Space;
  // typedef CutFESpace<Mesh> CutSpace;
  // typedef TestFunction<2> FunTest;
  // typedef FunFEM<Mesh2> Fun_h;

  const double cpubegin = CPUtime();

  double errL2 = 0.;

  // BACKGROUND MESH
  // Mesh Kh(nx, nx, -2, -2, 4., 4.);
  Space Vh(Kh, DataFE<Mesh2>::P1);
  double meshSize = hi;//(4./(nx-1));

  // TIME PARAMETER
  int divisionMeshsize = 4;
  // double dT = 0.01;//meshSize/divisionMeshsize;
  GTime::total_number_iteration = int(tfinal/dT);
  dT = tfinal / GTime::total_number_iteration;
  GTime::time_step = dT;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;



  //  -----------------------------------------------------
  //             Create files to save results
  //  -----------------------------------------------------
  // std::string pathOutpuFolder = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  // std::string pathOutpuFigure = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  // CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  // cout << " path to the output files : "+pathOutpuFolder << std::endl;
  // cout << " path to the vtk files : "+pathOutpuFigure << std::endl;
  std::ofstream outputData("dataDrop"+to_string(i)+".dat", std::ofstream::out);
  // cout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  // cout << " Creating the file \"output.txt\" to save cout" << std::endl;


  //CREATE SPACE AND VELOCITY FUNCTION
  Lagrange2 FEvelocity(2);
  Space VelVh  (Kh, FEvelocity);
  vector<Fun_h> vel(nbTime);


  // CREATE THE LEVELSET
  Space Lh  (Kh, DataFE<Mesh2>::P1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls(nbTime);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
  // LevelSet2 levelSet(Lh_k);

  // TIME INTERFACE AT ALL QUADRATURE TIME
  TimeInterface<Mesh> interface(qTime);

  // INITIALIZE THE PROBLEM
  CutFEM<Mesh2> surfactant(qTime);
  const R epsilon_S = 1.;

  int iter = 0, iterfig = 0;
  double q_init0, q_init1, qp1;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
    double tid = GTime::current_time();
    // INITIALIZE LEVELSET, VELOCITY, INTERFACE
    ls.begin()->swap(ls[nbTime-1]);
    for(int i=0;i<nbTime;++i) {
      R tt = In.Pt(R1(qTime(i).x));
      ls[i].init(Lh, fun_levelSet, tt);
      vel[i].init(VelVh, fun_velocity, tt);
      interface.init(i,Kh,ls[i]);
    }

    // CREATE THE ACTIVE MESH AND CUT SPACE
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.createSurfaceMesh(interface);
    CutSpace Wh(Kh0, Vh);

    // INITIALIZE THE PROBLEM
    Rn datau0;
    surfactant.initSpace(Wh, In);
    surfactant.initialSolution(datau0);
    KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
    if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
    Rn uh(datau0);
    Fun_h u0(Wh, datau0);

    Normal n;
    FunTest s(Wh,1), r(Wh,1);


    std::cout << " Begin assembly " << std::endl;
    double tt0 = MPIcf::Wtime();

    surfactant.addBilinear(
      -innerProduct(s, dt(r))
      + innerProduct(epsilon_S*gradS(s), gradS(r))
      , interface
      , In
    );

    for(int i=0;i<nbTime;++i) {
      ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
      surfactant.addBilinear(
        - innerProduct(s, dx(r)*vx)
        - innerProduct(s, dy(r)*vy)
        , interface
        , In
        , i
      );
    }
    // surfactant.addEdgeIntegral(
    //   innerProduct(1e-2*jump(grad(s).t()*n), jump(grad(r).t()*n)),
    //   In
    // );


    Fun_h funrhs (Wh, In, fun_rhs);
    // Fun_h funrhsp(cutVh, In);
    // projection(funrhs, funrhsp, In, interface ,0);
    surfactant.addLinear (
      innerProduct(funrhs.expression(), r)
      , interface
      , In
    );
    double intF   = integral(funrhs, In, interface, 0) ;
    // intFnp = integralSurf(funrhs , In, qTime);

    surfactant.addBilinear(
      innerProduct(s, r)
      , *interface(lastQuadTime)
      , In
      , lastQuadTime
    );
    surfactant.addLinear(
      innerProduct(u0.expression(), r)
      , *interface(0)
      , In
      , 0
    );

    // FunTest D2un = grad(grad(r)*n)*n;
    surfactant.addBilinear(
      innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
      , Kh0
      , innerFacet
      , In
    );
    // if stab interface, can have h in both stabilization
    // surfactant.addBilinear(
    //   innerProduct(1e-1*grad(s)*n, grad(r)*n)
    //   , interface
    //   , In
    // );


    // matlab::Export(surfactant.mat_, "matA.dat");
    // surfactant.cleanMatrix();
    surfactant.solve();

    KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
    uh = dw;
    surfactant.saveSolution(uh);



    Rn sol(Wh.get_nb_dof(), 0.);
    sol += uh(SubArray(Wh.get_nb_dof(), 0));
    sol += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
    Fun_h funuh_0(Wh, uh);
    Fun_h funuh_1(Wh, sol);

    //  -----------------------------------------------------
    //                COMPUTE ERROR
    //  -----------------------------------------------------
    {
      Fun_h funuh(Wh, sol);
      errL2 =  L2normSurf(funuh_0, fun_sol_surfactant, interface(0), tid,0,1);
      std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

      errL2 =  L2normSurf(funuh_1, fun_sol_surfactant, interface(nbTime-1), tid+dT,0,1);
      std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
    }
    //  -----------------------------------------------------
    // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    //  -----------------------------------------------------
    {
      double q0 = integral(funuh_0, interface(0), 0);
      double q1 = integral(funuh_1, interface(lastQuadTime), 0);
      if(iter==0){
        q_init0 = q0;
        q_init1 = q1;
        qp1     = q1;
        q_init1 = integral(u0, interface(0), 0);
      }

      outputData << setprecision(10);
      outputData << std::setw(10) << std::setfill(' ') << GTime::current_time()
                 << std::setw(20) << std::setfill(' ') << q0
                 << std::setw(20) << std::setfill(' ') << q1
                 << std::setw(20) << std::setfill(' ') << intF
                 << std::setw(20) << std::setfill(' ') << fabs(q1 - qp1 - intF)
                 << std::endl;
                 qp1 = q1;
      }



    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    // if(MPIcf::IamMaster() ) {
    //     // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
    //     Fun_h sol(Wh, uh);
    //
    //     // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
    //     Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");
    //
    //     writer.add(sol , "surfactant", 0, 1);
    //     writer.add(ls[0], "levelSet",  0, 1);
    //     // writer.add(ls[2], "levelSet2", 0, 1);
    //   }


    // return 0.;
    iter += 1;
  }


  outputData.close();
  return errL2;
}
double solveDG_1(int i, const Mesh& Kh, int nx, double tfinal, double hi, double dT) {
  const int d = 2;
  // typedef Mesh2 Mesh;
  // typedef FESpace2 Space;
  // typedef CutFESpace<Mesh> CutSpace;
  // typedef TestFunction<2> FunTest;
  // typedef FunFEM<Mesh2> Fun_h;

  const double cpubegin = CPUtime();

  double errL2 = 0.;

  // BACKGROUND MESH
  // Mesh Kh(2*nx, nx, -2, -2, 8., 4.);
  Space Vh(Kh, DataFE<Mesh2>::P1);
  double meshSize = hi;//(4./(nx-1));
  // double hi = meshSize;

  // TIME PARAMETER
  int divisionMeshsize = 4;
  // double dT = hi/8;//meshSize/divisionMeshsize;
  GTime::total_number_iteration = int(tfinal/dT);
  dT = tfinal / GTime::total_number_iteration;
  GTime::time_step = dT;
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  //CREATE SPACE AND VELOCITY FUNCTION
  Lagrange2 FEvelocity(2);
  Space VelVh  (Kh, FEvelocity);
  vector<Fun_h> vel(nbTime);

  // CREATE THE LEVELSET
  Space Lh  (Kh, DataFE<Mesh2>::P1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls(nbTime);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
  // LevelSet2 levelSet(Lh_k);

  // TIME INTERFACE AT ALL QUADRATURE TIME
  TimeInterface<Mesh> interface(qTime);

  // INITIALIZE THE PROBLEM
  CutFEM<Mesh2> surfactant(qTime);
  const R epsilon_S = 1.;

  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ITERATION \t : \t" << iter << " / " << GTime::total_number_iteration -1 << std::endl;
    double tid = GTime::current_time();
    // INITIALIZE LEVELSET, VELOCITY, INTERFACE
    ls.begin()->swap(ls[nbTime-1]);
    for(int i=0;i<nbTime;++i) {
      R tt = In.Pt(R1(qTime(i).x));
      ls[i].init(Lh, fun_levelSet, tt);
      vel[i].init(VelVh, fun_velocity, tt);
      interface.init(i,Kh,ls[i]);
    }

    // CREATE THE ACTIVE MESH AND CUT SPACE
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.createSurfaceMesh(interface);
    CutSpace Wh(Kh0, Vh);

    // INITIALIZE THE PROBLEM
    Rn datau0;
    surfactant.initSpace(Wh, In);
    surfactant.initialSolution(datau0);
    KN_<double> datas0(datau0(SubArray(Wh.get_nb_dof(),0)));
    if(iter == 0)interpolate(Wh, datas0, fun_init_surfactant);
    Rn uh(datau0);
    Fun_h u0(Wh, datau0);

    Normal n;
    Conormal conormal;
    FunTest s(Wh,1), r(Wh,1);


    std::cout << " Begin assembly " << std::endl;
    double tt0 = MPIcf::Wtime();

    surfactant.addBilinear(
      innerProduct(dt(s), r)
      + innerProduct(epsilon_S*gradS(s), gradS(r))
      , interface
      , In
    );


    for(int i=0;i<nbTime;++i) {      // computation of the curvature
      ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
      // Interior edges surface convective term
      // surfactant.addBilinear(
      //   - innerProduct(s, dx(r)*vx + dy(r)*vy)
      //   //+ innerProduct(dx(s)*vx + dy(s)*vy, r)  // TEST
      //   + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
      //   , interface
      //   , In , i
      // );
      // surfactant.addBilinear(
      //   + innerProduct(average((vel[i] * conormal)*s, 0.5, -0.5), jump(r))
      //   //+ innerProduct(average((vel.at(i) * t)*s, 1, 1), average(r)) // TEST
      //   //- innerProduct(jump(s), average((vel.at(i) * t) * r, 0.5, -0.5))
      //   //+ innerProduct(10 * jump((vel.at(i)*t)*s), jump(r))
      //   + innerProduct(0*average(fabs(vel[i]*conormal)*s, 1, 1), jump(r))
      //   , interface
      //   , innerRidge
      //   , In , i
      // );
      surfactant.addBilinear(
        + innerProduct(dx(s)*vx + dy(s)*vy, r)*0.5
        - innerProduct(s, dx(r)*vx + dy(r)*vy)*0.5
        + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
        , interface
        , In , i
      );

      // "Variant 2"
      surfactant.addBilinear(
          + innerProduct(average((vel[i]*conormal)*s, 0.5, -0.5), jump(r))*0.5
          - innerProduct(jump(s), average((vel[i]*conormal)*r, 0.5, -0.5))*0.5
          + innerProduct(10*jump(s), jump(r))
          //+ innerProduct(lambdaBSdiv*average((vel[i]*t)*s, 1, 1), average((vel[i]*t)*r, 1, 1))                  // (3.16)
          , interface
          , innerRidge
          , In , i
      );

    }
    surfactant.addBilinear(
        // innerProduct(hi*jump(s), jump(s))
      + innerProduct(1*jump(grad(s)*n), jump(grad(r)*n))
      , Kh0
      , innerFacet
      , In
    );
    // surfactant.addBilinear(
    //   innerProduct(hi*grad(s)*n, grad(r)*n)
    //   , interface
    //   , In
    // );
        surfactant.addBilinear(
          innerProduct(s, r)
          , *interface(0)
          , In
          , 0
        );
        surfactant.addLinear(
          innerProduct(u0.expression(), r)
          , *interface(0)
          , In
          , 0
        );

        Fun_h funrhs(Wh, In, fun_rhs);
        // Fun_h funrhsp(cutVh, In);
        // projection(funrhs, funrhsp, In, interface ,0);
        surfactant.addLinear(
          innerProduct(funrhs.expression(), r)
          , interface
          , In
        );
        // double  intF = integralSurf(funrhsp, In, interface) ;
    //     intFnp = integralSurf(funrhs, In, qTime);
    //
    //
    // matlab::Export(surfactant.mat_, "matA.dat");
    // surfactant.cleanMatrix();
    surfactant.solve();

    KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
    uh = dw;
    surfactant.saveSolution(uh);
    //



    //  -----------------------------------------------------
    //                COMPUTE ERROR
    //  -----------------------------------------------------
    {
      Rn sol(Wh.get_nb_dof(), 0.);
      sol += uh(SubArray(Wh.get_nb_dof(), 0));
      Fun_h funuh(Wh, sol);
      errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(0), tid,0,1);
      std::cout << " t_n -> || u-uex||_2 = " << errL2 << std::endl;

      sol  += uh(SubArray(Wh.get_nb_dof(), Wh.get_nb_dof()));
      errL2 =  L2normSurf(funuh, fun_sol_surfactant, *interface(nbTime-1), tid+dT,0,1);
      std::cout << " t_{n+1} -> || u-uex||_2 = " << errL2 << std::endl;
    }



    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(MPIcf::IamMaster() ) {
      // KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
      Fun_h sol(Wh, uh);
      Paraview<Mesh> writer(Kh0, "surfactant_"+to_string(iter)+".vtk");

      writer.add(sol , "surfactant", 0, 1);
      writer.add(ls[0], "levelSet",  0, 1);
      // writer.add(ls[2], "levelSet2", 0, 1);
    }


    // return 0.;
    iter += 1;
  }

  return errL2;
}


int main(int argc, char** argv ) {


  MPIcf cfMPI(argc,argv);
  int nx = 20;
  for(int iter=0; iter<4;++iter) {

    Mesh Kh(nx, nx, -2, -2, 4., 4.);
    double h = 4./(nx-1);
    double dt = h/4;
    double err = solveCG_1(iter, Kh, nx, 0.25, h, dt);

    nx = 2*nx-1;
  }
  return 0;
}


  //
  //
  //
  // cout << " nx        \t" << nx << std::endl;
  // cout << " Mesh size \t" << meshSize << std::endl;
  // cout << " Time Step \t" << dT << std::endl;
  // cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  // cout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  // cout << " number of quadrature points in time : \t" << nbTime << std::endl;
  // cout << " number of dof in time per time slab : \t" << ndfTime << std::endl;
  //
  //
  // // Set parameters for paraview PLOTTING
  // const bool writeVTKFiles = true;
  // const bool saveSurfactantVTK = true;
  // const bool saveVelocityVTK = false;
  // const bool saveLevelSetVTK = false;
  // const int frequencyPlottingTheSolution = 10;




//   const CutFEM_Parameter& h(Parameter::h);
//   const CutFEM_Parameter& h_E(Parameter::meas);
//
//   double q0_0, q0_1, qp_0, qp_1;
//   double intF = 0, intFnp = 0;
//   int iter = 0, iterfig = 0;
//   while( iter < GTime::total_number_iteration ) {
//
//     GTime::current_iteration = iter;
//     const TimeSlab& In(Ih[iter]);
//
//     std::cout << " ------------------------------------------------------------- "<< std::endl;
//     std::cout << " ------------------------------------------------------------- "<< std::endl;
//     std::cout << " ITERATION \t : \t" << iter << std::endl;
//     std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;
//
//
//     ls_k.begin()->swap(ls_k[nbTime-1]);
//     // computation of the nbTime interfaces
//     for(int i=0;i<nbTime;++i) {
//
//       R tt = In.Pt(R1(qTime(i).x));
//       ls_k[i].init(Lh_k, fun_levelSet, tt);
//       ls[i].init(Lh, fun_levelSet, tt);
//
//       mapping[i] = new Mapping2(VelVh, ls_k[i]);
//
//       // projection(vel[i], vel[i]);
//       vel[i].init(VelVh, fun_velocity, tt);
//
//       // projection(ls_k[i], ls[i]);
//       interface.init(i,Th,ls[i].v);
//
//       // gnuplot::save(*interface[i], "interface"+to_string(i)+".dat");
//       // if(i<nbTime-1) {
//       //   LevelSet levelSet(ls_k[i], vel[i], vel[i+1], dt_levelSet);
//       //   ls_k[i+1].init(levelSet.rhs);
//       //   if(iter%frequencyReinitialization == 0 && i == 0)
//       //    reinitialization.perform(ls_k[i], ls_k[i], *interface[i]);
//       // }
//     }
//
//     // vel[i].init(VelVh, fun_velocity, In,tt);
//
//
//     // Create the Active Cut Mesh for insoluble surfactant
//     Mesh2 cutTh(interface);
// 	  FESpace2 cutVh(cutTh, interface, DataFE<Mesh2>::P1);
//     cutVh.backSpace = &Lh_k;  // save backSpace to save solution
//
//     // FESpace2 cutVh2(cutTh, interface, FEvelocity);
//     // cutVh2.backSpace = &VelVh;
//     // FESpace2 VelVh  (Th, FEvelocity);
//     // Fun_h vel2(cutVh2, In, fun_velocity);
//
//
//     surfactant.initSpace(cutVh, In);
//
//     Rn datau0;
//     surfactant.initialSolution(datau0);
//
//     KN_<double> datas0(datau0(SubArray(cutVh.NbDoF(),0)));
//     if(iter == 0)interpolate(cutVh, datas0, fun_init_surfactant);
//     Rn uh(datau0);
//     Fun_h u0(cutVh, datau0);
//     Normal n;
//     FunTest s(cutVh,1), r(cutVh,1);
//     cout << " Number of dof of the problem \t" << surfactant.nDoF << std::endl;
//
//
//
//     std::cout << " Begin assembly " << std::endl;
//     double tt0 = CPUtime();
//
//
//
// #ifdef FORMULATION1
//
//     surfactant.addBilinear(
//           // innerProduct(dt(s), r)
//         + innerProduct(epsilon_S*gradS(s), gradS(r))
//         , interface
//         , In
//     );
//
//     // for(int i=0;i<nbTime;++i) {      // computation of the curvature
//     //   ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
//     //   ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
//     //   // ExpressionFunFEM<Mesh2> vx(vel,0,op_id);
//     //   // ExpressionFunFEM<Mesh2> vy(vel,1,op_id);
//     //   surfactant.addBilinear(
//     //         innerProduct(dx(s)*vx + dy(s)*vy, r)
//     //       + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
//     //       , interface
//     //       ,i , In
//     //   );
//     // }
//
//
//     surfactant.addEdgeIntegral(
//       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n)),
//       In
//     );
//
//     // surfactant.addBilinearFormInterface(
//     //   innerProduct(1e-2*grad(s)*n, grad(r)*n)
//     //   ,In
//     // );
//
//     double cc = 1./In.T.mesure()*6.;
//     surfactant.addBilinear(
//       innerProduct(cc*s, r)
//       , interface
//       , 0
//       , In
//     );
//     surfactant.addLinear(
//       innerProduct(u0.expression(), cc*r)
//       , interface
//       , 0
//       , In
//     );
//
//     Fun_h funrhs(cutVh, In, fun_rhs);
//     Fun_h funrhsp(cutVh, In);
//     projection(funrhs, funrhsp, In, interface ,0);
//     surfactant.addLinear(
//       innerProduct(funrhs.expression(), r)
//       , interface
//       , In
//     );
//     intF = integralSurf(funrhsp, In, qTime) ;
//     intFnp = integralSurf(funrhs, In, qTime);
//
//
//     // for(int i=0;i<nbTime;++i) {
//     //
//     //   const R1 tq = In.map(qTime(i).x);
//     //   Fun_h funrhs0(cutVh, fun_rhs, tq.x);
//     //   Fun_h funrhsp(cutVh, fun_rhs, tq.x);
//     //
//     //   projection(funrhs0, funrhsp, *interface[i] );
//     //
//     //   // Fun_h funuh(cutVh, uh);
//     //   // std::cout << integralSurf(funrhs0,0,i) << "\t"
//     //   //           << integralSurf(funrhsp,0,i) << std::endl;
//     //
//     //   surfactant.addLinearFormInterface (
//     //     innerProduct(funrhsp.expression(), r)
//     //     , i, In
//     //   );
//     // }
//
// #else
//
//     surfactant.addBilinear(
//         -innerProduct(s, dt(r))
//         + innerProduct(epsilon_S*gradS(s), gradS(r))
//         , interface
//         , In
//         // , mapping
//     );
//
//     for(int i=0;i<nbTime;++i) {
//       ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
//       ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
//       surfactant.addBilinear(
//         - innerProduct(s, dx(r)*vx)
//         - innerProduct(s, dy(r)*vy)
//         , interface
//         , i
//         , In
//         // , mapping
//       );
//     }
//     // surfactant.addEdgeIntegral(
//     //   innerProduct(1e-2*jump(grad(s).t()*n), jump(grad(r).t()*n)),
//     //   In
//     // );
//
//
//     Fun_h funrhs (cutVh, In, fun_rhs);
//     Fun_h funrhsp(cutVh, In);
//     projection(funrhs, funrhsp, In, interface ,0);
//     surfactant.addLinear (
//       innerProduct(funrhs.expression(), r)
//       , interface
//       , In
//       // , mapping
//     );
//     intF   = integralSurf(funrhsp, In, qTime) ;
//     intFnp = integralSurf(funrhs , In, qTime);
//
//
//     double ccend = 1./In.T.mesure()*1./qTime[lastQuadTime].a;;
//
//     surfactant.addBilinear(
//       innerProduct(ccend*s, r)
//       , interface
//       , lastQuadTime
//       , In
//       // , mapping
//     );
//     double cc0 = 1./In.T.mesure()*1./qTime[0].a;;
//     surfactant.addLinear(
//       innerProduct(u0.expression(), cc0*r)
//       , interface
//       , 0
//       , In
//       // , mapping
//     );
//
//     FunTest D2un = grad(grad(r)*n)*n;
//
//     surfactant.addEdgeIntegral(
//       innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
//       // + innerProduct(1e-2*h_E*jump(D2un), h_E*jump(D2un)),
//       , In
//       // , mapping
//     );
//     // if stab interface, can have h in both stabilization
//     // surfactant.addBilinear(
//     //   innerProduct(1e-2*h_E*grad(s)*n, grad(r)*n)
//     //   , interface
//     //   , In
//     // );
//
//
//
// #endif
//     surfactant.solve();
//
//     KN_<double> dw(surfactant.rhs(SubArray(surfactant.nDoF, 0)));
//     uh = dw;
//     surfactant.saveSolution(uh);
//
//
//
//     // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
//     {
//       cout << "\n Features of the drop " << std::endl;
//       Fun_h funuh(cutVh, uh);
//
//       Rn sol2(cutVh.NbDoF(), 0.);
//       Fun_h funsol(cutVh, sol2);
//       sol2 += uh(SubArray(cutVh.NbDoF(), 0));
//       double q_0 = integralSurf(funsol, 0);
//       sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
//       double q_1 = integralSurf(funsol, 0, lastQuadTime);
//       if(iter==0){ q0_0 = q_0;q0_1 = q_1;
//                    qp_1 = q_1;
//                    q0_1 = integralSurf(u0, 0);
//                  }
//
//
//       // cout << " q0 - qt_1 \t" << q_1-qp_1 << std::endl;
//       // cout << " int F \t" << intF << std::endl;
//
//       outputData << setprecision(10);
//       outputData << GTime::current_time() << "\t"
//                  << (q_1-qp_1) << "\t"
//                  << fabs(q_1-qp_1) << "\t"
//                  << intF << "\t"
//                  << intFnp << "\t"
//                  << ((q_1 -qp_1) - intFnp) << "\t"
//                  << q0_1 - q_1 << "\t"
//                  << std::endl;
//       qp_1 = q_1;
//     }
//
//
//     {
//       Fun_h funsol(cutVh, fun_sol_surfactant, GTime::current_time());
//       funsol.v -= uh(SubArray(cutVh.NbDoF(), 0));
//
//       Rn sol2(cutVh.NbDoF(), 0.);
//       sol2 += uh(SubArray(cutVh.NbDoF(), 0));
//       // sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
//
//       ExpressionFunFEM<Mesh> sol(funsol, 0, op_id);
//       // cout << " error    ->    " << sqrt(integralSurf(sol*sol, cutVh, 0)) << std::endl;
//       // Fun_h funuh(cutVh, uh);
//       Fun_h funuh(cutVh, sol2);
//       // cout << " error    ->    " << L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1) << std::endl;
//       errL2 =  L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1);
//       std::cout << " || u-uex||_2 = " << errL2 << std::endl;
//     }
//
//
//   //  -----------------------------------------------------
//   //                     PLOTTING
//   //  -----------------------------------------------------
//   if(MPIcf::IamMaster() && iter%frequencyPlottingTheSolution == 0 && writeVTKFiles){
//     std::cout << " Plotting " << std::endl;
//
//     if(saveSurfactantVTK) {
//       KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
//       Fun_h sol(cutVh, sols);
//
//       // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
//       Paraview2 writer(cutVh, "surfactant_"+to_string(iterfig)+".vtk");
//
//       writer.add(sol , "surfactant", 0, 1);
//       writer.add(ls[0], "levelSet",  0, 1);
//       writer.add(ls[2], "levelSet2", 0, 1);
//     }
//     // if(saveLevelSetVTK) {
//     //   Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
//     //   writerLS.add(ls[0], "levelSet", 0, 1);
//     // }
//
//     iterfig++;
//   }
//   for(int i=0;i<nbTime;++i) delete mapping[i];
//   iter++;
//
//   }
//
//   myCout << 4./nx << "\t" << errL2 << std::endl;
//   nx *= 2;
//   ny *= 2;
// }
// return 0;
// }
