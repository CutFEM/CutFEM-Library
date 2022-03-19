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
#include "../num/gnuplot.hpp"
#include "projection.hpp"
#include "../util/redirectOutput.hpp"

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
  R fun_sol_surfactant(const R2 P,  const int i, const R t) {
    return P.y/1+2;;
  }
  R fun_rhs(const R2 P, const int cc, const R t) {
    return 0.;
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
  R fun_sol_surfactant(const R2 P,  const int i, const R t) {return P.x*P.y*exp(-4*t);}
  R fun_velocity(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    R b = M_PI/4*(cos(2*M_PI*t))/a;
    return (i==0) ? b*P.x : 0;
  }
  R fun_levelSet(const R2 P, const int i, const R t) {
    R a = (1+0.25*sin(2*M_PI*t));
    // a = 1;
    return sqrt(P.x*P.x/a + P.y*P.y) - 1;
  }
  R fun_rhs(const R2 P, const int cc, const R t) {
    R x = P.x,  y = P.y;
    R r = (4*x*y*exp(-4*t)*(sin(2*M_PI*t) + 4)*
           (2048*x*x*y*y + 336*y*y*y*y*pow(sin(2*M_PI*t),2) +
            52*y*y*y*y*pow(sin(2*M_PI*t),3) + 3*y*y*y*y*pow(sin(2*M_PI*t),4)
            + 64*x*x*x*x*sin(2*M_PI*t) + 960*y*y*y*y*sin(2*M_PI*t) + 1024*x*x*x*x
            + 1024*y*y*y*y + 144*x*x*y*y*pow(sin(2*M_PI*t),2)
            + 4*x*x*y*y*pow(sin(2*M_PI*t),3) + 1024*x*x*y*y*sin(2*M_PI*t)))
      /pow((y*y*pow(sin(2*t*M_PI),2) + 8*y*y*sin(2*t*M_PI) + 16*x*x + 16*y*y),3)
      - 4*x*y*exp(-4*t) + (x*y*M_PI*exp(-4*t)*cos(2*M_PI*t))/(sin(2*M_PI*t) + 4)
      + (x*y*y*y*M_PI*exp(-4*t)*cos(2*M_PI*t)*(sin(2*M_PI*t) + 4))/
      (y*y*pow(sin(2*M_PI*t),2) + 8*y*y*sin(2*M_PI*t) + 16*x*x + 16*y*y);

    return r;
  }
}

#define FORMULATION1;
// #define FORMULATION2;

using namespace NumericSurfactantEllipse2D ;
// using namespace NumericSurfactant2D0 ;
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
  int ny = 40;

  std::string pathOutpuFolder = "../../outputFiles/surfactant2/Ellipse/F1h40/";
  std::string pathOutpuFigure = "../../outputFiles/surfactant2/Ellipse/F1h40/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");

  for( int iii=0;iii<1;++iii){
    double errL2 = 0.;

  Mesh2 Th(nx, ny, -2, -2, 4., 4.);
  double meshSize = (4./nx);
  int divisionMeshsize = 4;
  double dT = meshSize/divisionMeshsize;
  double tfinal = dT;//0.25;
  GTime::total_number_iteration = int(tfinal/dT);
  dT = tfinal / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());


  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  // Create files to save results
  // std::string meshstr = to_string(nx)+to_string(ny);
  // std::string pathOutpuFolder = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  // std::string pathOutpuFigure = "../../outputFiles/surfactant2/Ellipse/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  // CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  cout << " path to the output files : "+pathOutpuFolder << std::endl;
  cout << " path to the vtk files : "+pathOutpuFigure << std::endl;

  std::ofstream outputData(pathOutpuFolder+"dataDrop_F13P.dat", std::ofstream::out);
  cout << " Creating the file \"dataDrop.dat\" to save data" << std::endl;
  cout << " Creating the file \"output.txt\" to save cout" << std::endl;




  cout << " nx        \t" << nx << std::endl;
  cout << " Mesh size \t" << meshSize << std::endl;
  cout << " Time Step \t" << dT << std::endl;
  cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  cout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  cout << " number of quadrature points in time : \t" << nbTime << std::endl;
  cout << " number of dof in time per time slab : \t" << ndfTime << std::endl;


  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveSurfactantVTK = true;
  const bool saveVelocityVTK = false;
  const bool saveLevelSetVTK = false;
  const int frequencyPlottingTheSolution = 10;

  Lagrange2 FEvelocity(2);
  FESpace2 VelVh  (Th, FEvelocity);
  vector<Fun_h> vel(nbTime);

  // levelSet stuff
  // ----------------------------------------------
  cout << "\n --------------------------------------- \n " << std::endl;
  cout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
  FESpace2 Lh  (Th, DataFE<Mesh2>::P1);
  FESpace2 Lh_k(Th, DataFE<Mesh2>::P2);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  for(int i=0;i<nbTime;++i) ls_k[i].init(Lh_k, fun_levelSet, 0.);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);
  projection(ls_k[0], ls[nbTime-1]);
  LevelSet2 levelSet(Lh_k);


  CReinitialization<Mesh2> reinitialization;
  reinitialization.number_iteration = 10;
  reinitialization.epsilon_diffusion = 1e-3;
  reinitialization.dt = dT/8;
  reinitialization.ON_OFF = "OFF";
  reinitialization.ODE_method = "Euler";
  reinitialization.mass_correction = "OFF";
  reinitialization.precision_correction = 1e-8;
  reinitialization.max_iteration_correction = 10;
  reinitialization.info();
  const int frequencyReinitialization = 1000;

  TimeInterface2 interface(nbTime);
  KN<Mapping2*> mapping(nbTime);

  CutFEM<Mesh2> surfactant(qTime);
  const R epsilon_S = 1.;
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  double q0_0, q0_1, qp_0, qp_1;
  double intF = 0, intFnp = 0;
  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;


    ls_k.begin()->swap(ls_k[nbTime-1]);
    // computation of the nbTime interfaces
    for(int i=0;i<nbTime;++i) {

      R tt = In.Pt(R1(qTime(i).x));
      ls_k[i].init(Lh_k, fun_levelSet, tt);
      ls[i].init(Lh, fun_levelSet, tt);

      mapping[i] = new Mapping2(VelVh, ls_k[i]);

      // projection(vel[i], vel[i]);
      vel[i].init(VelVh, fun_velocity, tt);

      // projection(ls_k[i], ls[i]);
      interface.init(i,Th,ls[i].v);

      // gnuplot::save(*interface[i], "interface"+to_string(i)+".dat");
      // if(i<nbTime-1) {
      //   LevelSet levelSet(ls_k[i], vel[i], vel[i+1], dt_levelSet);
      //   ls_k[i+1].init(levelSet.rhs);
      //   if(iter%frequencyReinitialization == 0 && i == 0)
      //    reinitialization.perform(ls_k[i], ls_k[i], *interface[i]);
      // }
    }

    // vel[i].init(VelVh, fun_velocity, In,tt);


    // Create the Active Cut Mesh for insoluble surfactant
    Mesh2 cutTh(interface);
	  FESpace2 cutVh(cutTh, interface, DataFE<Mesh2>::P1);
    cutVh.backSpace = &Lh_k;  // save backSpace to save solution

    // FESpace2 cutVh2(cutTh, interface, FEvelocity);
    // cutVh2.backSpace = &VelVh;
    // FESpace2 VelVh  (Th, FEvelocity);
    // Fun_h vel2(cutVh2, In, fun_velocity);


    surfactant.initSpace(cutVh, In);

    Rn datau0;
    surfactant.initialSolution(datau0);

    KN_<double> datas0(datau0(SubArray(cutVh.NbDoF(),0)));
    if(iter == 0)interpolate(cutVh, datas0, fun_init_surfactant);
    Rn uh(datau0);
    Fun_h u0(cutVh, datau0);
    Normal n;
    FunTest s(cutVh,1), r(cutVh,1);
    cout << " Number of dof of the problem \t" << surfactant.nDoF << std::endl;



    std::cout << " Begin assembly " << std::endl;
    double tt0 = CPUtime();



#ifdef FORMULATION1

    surfactant.addBilinear(
          // innerProduct(dt(s), r)
        + innerProduct(epsilon_S*gradS(s), gradS(r))
        , interface
        , In
    );

    // for(int i=0;i<nbTime;++i) {      // computation of the curvature
    //   ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
    //   ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
    //   // ExpressionFunFEM<Mesh2> vx(vel,0,op_id);
    //   // ExpressionFunFEM<Mesh2> vy(vel,1,op_id);
    //   surfactant.addBilinear(
    //         innerProduct(dx(s)*vx + dy(s)*vy, r)
    //       + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
    //       , interface
    //       ,i , In
    //   );
    // }


    surfactant.addEdgeIntegral(
      innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n)),
      In
    );

    // surfactant.addBilinearFormInterface(
    //   innerProduct(1e-2*grad(s)*n, grad(r)*n)
    //   ,In
    // );

    double cc = 1./In.T.mesure()*6.;
    surfactant.addBilinear(
      innerProduct(cc*s, r)
      , interface
      , 0
      , In
    );
    surfactant.addLinear(
      innerProduct(u0.expression(), cc*r)
      , interface
      , 0
      , In
    );

    Fun_h funrhs(cutVh, In, fun_rhs);
    Fun_h funrhsp(cutVh, In);
    projection(funrhs, funrhsp, In, interface ,0);
    surfactant.addLinear(
      innerProduct(funrhs.expression(), r)
      , interface
      , In
    );
    intF = integralSurf(funrhsp, In, qTime) ;
    intFnp = integralSurf(funrhs, In, qTime);


    // for(int i=0;i<nbTime;++i) {
    //
    //   const R1 tq = In.map(qTime(i).x);
    //   Fun_h funrhs0(cutVh, fun_rhs, tq.x);
    //   Fun_h funrhsp(cutVh, fun_rhs, tq.x);
    //
    //   projection(funrhs0, funrhsp, *interface[i] );
    //
    //   // Fun_h funuh(cutVh, uh);
    //   // std::cout << integralSurf(funrhs0,0,i) << "\t"
    //   //           << integralSurf(funrhsp,0,i) << std::endl;
    //
    //   surfactant.addLinearFormInterface (
    //     innerProduct(funrhsp.expression(), r)
    //     , i, In
    //   );
    // }

#else

    surfactant.addBilinear(
        -innerProduct(s, dt(r))
        + innerProduct(epsilon_S*gradS(s), gradS(r))
        , interface
        , In
        // , mapping
    );

    for(int i=0;i<nbTime;++i) {
      ExpressionFunFEM<Mesh2> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh2> vy(vel[i],1,op_id);
      surfactant.addBilinear(
        - innerProduct(s, dx(r)*vx)
        - innerProduct(s, dy(r)*vy)
        , interface
        , i
        , In
        // , mapping
      );
    }
    // surfactant.addEdgeIntegral(
    //   innerProduct(1e-2*jump(grad(s).t()*n), jump(grad(r).t()*n)),
    //   In
    // );


    Fun_h funrhs (cutVh, In, fun_rhs);
    Fun_h funrhsp(cutVh, In);
    projection(funrhs, funrhsp, In, interface ,0);
    surfactant.addLinear (
      innerProduct(funrhs.expression(), r)
      , interface
      , In
      // , mapping
    );
    intF   = integralSurf(funrhsp, In, qTime) ;
    intFnp = integralSurf(funrhs , In, qTime);


    double ccend = 1./In.T.mesure()*1./qTime[lastQuadTime].a;;

    surfactant.addBilinear(
      innerProduct(ccend*s, r)
      , interface
      , lastQuadTime
      , In
      // , mapping
    );
    double cc0 = 1./In.T.mesure()*1./qTime[0].a;;
    surfactant.addLinear(
      innerProduct(u0.expression(), cc0*r)
      , interface
      , 0
      , In
      // , mapping
    );

    FunTest D2un = grad(grad(r)*n)*n;

    surfactant.addEdgeIntegral(
      innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n))
      // + innerProduct(1e-2*h_E*jump(D2un), h_E*jump(D2un)),
      , In
      // , mapping
    );
    // if stab interface, can have h in both stabilization
    // surfactant.addBilinear(
    //   innerProduct(1e-2*h_E*grad(s)*n, grad(r)*n)
    //   , interface
    //   , In
    // );



#endif
    surfactant.solve();

    KN_<double> dw(surfactant.rhs(SubArray(surfactant.nDoF, 0)));
    uh = dw;
    surfactant.saveSolution(uh);



    // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
    {
      cout << "\n Features of the drop " << std::endl;
      Fun_h funuh(cutVh, uh);

      Rn sol2(cutVh.NbDoF(), 0.);
      Fun_h funsol(cutVh, sol2);
      sol2 += uh(SubArray(cutVh.NbDoF(), 0));
      double q_0 = integralSurf(funsol, 0);
      sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
      double q_1 = integralSurf(funsol, 0, lastQuadTime);
      if(iter==0){ q0_0 = q_0;q0_1 = q_1;
                   qp_1 = q_1;
                   q0_1 = integralSurf(u0, 0);
                 }


      // cout << " q0 - qt_1 \t" << q_1-qp_1 << std::endl;
      // cout << " int F \t" << intF << std::endl;

      outputData << setprecision(10);
      outputData << GTime::current_time() << "\t"
                 << (q_1-qp_1) << "\t"
                 << fabs(q_1-qp_1) << "\t"
                 << intF << "\t"
                 << intFnp << "\t"
                 << ((q_1 -qp_1) - intFnp) << "\t"
                 << q0_1 - q_1 << "\t"
                 << std::endl;
      qp_1 = q_1;
    }


    {
      Fun_h funsol(cutVh, fun_sol_surfactant, GTime::current_time());
      funsol.v -= uh(SubArray(cutVh.NbDoF(), 0));

      Rn sol2(cutVh.NbDoF(), 0.);
      sol2 += uh(SubArray(cutVh.NbDoF(), 0));
      // sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));

      ExpressionFunFEM<Mesh> sol(funsol, 0, op_id);
      // cout << " error    ->    " << sqrt(integralSurf(sol*sol, cutVh, 0)) << std::endl;
      // Fun_h funuh(cutVh, uh);
      Fun_h funuh(cutVh, sol2);
      // cout << " error    ->    " << L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1) << std::endl;
      errL2 =  L2normSurf(funuh, fun_sol_surfactant,GTime::current_time(),0,1);
      std::cout << " || u-uex||_2 = " << errL2 << std::endl;
    }


  //  -----------------------------------------------------
  //                     PLOTTING
  //  -----------------------------------------------------
  if(MPIcf::IamMaster() && iter%frequencyPlottingTheSolution == 0 && writeVTKFiles){
    std::cout << " Plotting " << std::endl;

    if(saveSurfactantVTK) {
      KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
      Fun_h sol(cutVh, sols);

      // Paraview2 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
      Paraview2 writer(cutVh, "surfactant_"+to_string(iterfig)+".vtk");

      writer.add(sol , "surfactant", 0, 1);
      writer.add(ls[0], "levelSet",  0, 1);
      writer.add(ls[2], "levelSet2", 0, 1);
    }
    // if(saveLevelSetVTK) {
    //   Paraview2 writerLS(Lh, pathOutpuFigure+"levelSet_"+to_string(iterfig)+".vtk");
    //   writerLS.add(ls[0], "levelSet", 0, 1);
    // }

    iterfig++;
  }
  for(int i=0;i<nbTime;++i) delete mapping[i];
  iter++;

  }

  myCout << 4./nx << "\t" << errL2 << std::endl;
  nx *= 2;
  ny *= 2;
}
return 0;
}
