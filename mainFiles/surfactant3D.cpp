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
namespace NumericSurfactant3D0 {
  R fun_levelSet(const R3 P, const int i, const double t) {
    R x = P.x,  y = P.y, z = P.z;
    return (x-z*z)*(x-z*z) + y*y + z*z - 1;
  }
  R fun_velocity(const R3 P, int i, double t) {
    R x = P.x,  y = P.y, z = P.z;
    R3 Q(0.1*x*cos(t),0.2*y*sin(t),0.2*z*cos(t));
    return Q[i];
  }
  R fun_rhs(const R3 P, const int i, const R t) {
    return 0;
  }
  R fun_init_surfactant(const R3 P, const int i) {
    R x = P.x,  y = P.y, z = P.z;
    return 1 + x*y*z;
  }

}

//#define FORMULATION1;
#define FORMULATION2;

using namespace NumericSurfactant3D0 ;

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

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;
  int nz = 20;

  Mesh3 Th(nx, ny, nz, -2, -2, -2, 4., 4., 4.);
  double meshSize = (4./nx);
  int divisionMeshsize = 2;
  double dT = meshSize/divisionMeshsize;
  double tfinal = 4;
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
  std::string meshstr = to_string(nx)+to_string(ny)+to_string(nz);
  std::string pathOutpuFolder = "../../outputFiles/surfactant3/ConservationExampleTest/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/";
  std::string pathOutpuFigure = "../../outputFiles/surfactant3/ConservationExampleTest/mesh"+meshstr+"/h"+to_string(divisionMeshsize)+"/paraview/";
  if(MPIcf::IamMaster()) {
    std::experimental::filesystem::create_directories(pathOutpuFigure);
  }
  CoutFileAndScreen myCout(pathOutpuFolder+"output.txt");
  myCout << " path to the output files : "+pathOutpuFolder << std::endl;
  myCout << " path to the vtk files : "+pathOutpuFigure << std::endl;

  std::ofstream outputData(pathOutpuFolder+"dataDrop1_3D_20_dt2.dat", std::ofstream::out);
  myCout << " Creating the file \"dataDrop1.dat\" to save data" << std::endl;
  myCout << " Creating the file \"output.txt\" to save cout" << std::endl;

  myCout << " nx        \t" << nx << std::endl;
  myCout << " Mesh size \t" << meshSize << std::endl;
  myCout << " Time Step \t" << dT << std::endl;
  myCout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  myCout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  myCout << " number of quadrature points in time : \t" << nbTime << std::endl;
  myCout << " number of dof in time per time slab : \t" << ndfTime << std::endl;


  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveSurfactantVTK = true;
  const int frequencyPlottingTheSolution = 1;

  // Space for the velocity field
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " The velocity is interpolated to a P2 space " << std::endl;
  Lagrange3 FEvelocity(2);
  FESpace3 VelVh  (Th, FEvelocity);
  vector<Fun_h> vel(nbTime);


  // levelSet stuff
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
  FESpace Lh  (Th, DataFE<Mesh>::P1);
  // FESpace Lh_k(Th, DataFE<Mesh>::P1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  // for(int i=0;i<nbTime;++i) ls_k[i].init(Lh_k, fun_levelSet, 0.);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet, 0.);

  TimeInterface3 interface(nbTime);


  CutFEM<Mesh> surfactant(qTime);
  const R epsilon_S = 1;
  const CutFEM_Parameter& h(Parameter::h);
  const CutFEM_Parameter& h_E(Parameter::meas);

  double q0_0, q0_1, qp_0, qp_1;
  double intF = 0, intFnp = 0;
  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {

    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    myCout << " ------------------------------------------------------------- "<< std::endl;
    myCout << " ITERATION \t : \t" << iter << std::endl;
    myCout << " TIME      \t : \t" << GTime::current_time() << std::endl;


    for(int i=0;i<nbTime;++i) {
      R tt = In.Pt(R1(qTime(i).x));
      // ls_k[i].init(Lh_k, fun_levelSet, tt);
      ls[i].init(Lh, fun_levelSet, tt);
      vel[i].init(VelVh, fun_velocity, tt);

      interface.init(i,Th,ls[i].v);
    }

    // Create the Active Cut Mesh for insoluble surfactant
    Mesh3 cutTh(interface);
	  FESpace3 cutVh(cutTh, interface, DataFE<Mesh3>::P1);
    cutVh.backSpace = &Lh;  // save backSpace to save solution

    surfactant.initSpace(cutVh, In);

    Rn datau0;
    surfactant.initialSolution(datau0);

    KN_<double> datas0(datau0(SubArray(cutVh.NbDoF(),0)));
    if(iter == 0)interpolate(cutVh, datas0, fun_init_surfactant);

    Rn uh(datau0);
    Fun_h u0(cutVh, datau0);
    Normal n;
    FunTest s(cutVh,1), r(cutVh,1);
    std::cout << " Number of dof of the problem \t" << surfactant.nDoF << std::endl;

    std::cout << " Begin assembly " << std::endl;
    double tt0 = CPUtime();

    #ifdef FORMULATION1

    surfactant.addBilinear(
      innerProduct(dt(s), r)
      + innerProduct(epsilon_S*gradS(s), gradS(r))
      , interface
      , In
    );

    for(int i=0;i<nbTime;++i) {
      ExpressionFunFEM<Mesh3> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh3> vy(vel[i],1,op_id);
      surfactant.addBilinear(
        innerProduct(dxS(s)*vx + dyS(s)*vy, r)
        + innerProduct(s*dxS(vel[i]) + s*dyS(vel[i]), r)
        , interface
        ,i , In
      );
    }

    surfactant.addEdgeIntegral(
      innerProduct(1e-2*h_E*jump(grad(s)*n), jump(grad(r)*n)),
      In
    );
    surfactant.addBilinear(
      innerProduct(1e-2*h_E*grad(s)*n, grad(r)*n)
      , interface
      , In
    );
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

    // Fun_h funrhs(cutVh, In, fun_rhs);
    // Fun_h funrhsp(cutVh, In);
    // projection(funrhs, funrhsp, In, interface ,0);
    // surfactant.addLinear(
    //   innerProduct(funrhs.expression(), r)
    //   , interface
    //   , In
    // );
    // intF = integralSurf(funrhsp, In, qTime) ;
    // intFnp = integralSurf(funrhs, In, qTime);


    #else
    surfactant.addBilinear(
      -innerProduct(s, dt(r))
      + innerProduct(epsilon_S*gradS(s), gradS(r))
      , interface
      , In
    );

    for(int i=0;i<nbTime;++i) {
      ExpressionFunFEM<Mesh3> vx(vel[i],0,op_id);
      ExpressionFunFEM<Mesh3> vy(vel[i],1,op_id);
      ExpressionFunFEM<Mesh3> vz(vel[i],2,op_id);
      surfactant.addBilinear(
        - innerProduct(s, dxS(r)*vx)
        - innerProduct(s, dyS(r)*vy)
        - innerProduct(s, dzS(r)*vz)
        , interface
        , i
        , In
      );
    }

    // Fun_h funrhs (cutVh, In, fun_rhs);
    // Fun_h funrhsp(cutVh, In);
    // projection(funrhs, funrhsp, In, interface ,0);
    // surfactant.addLinear (
    //   innerProduct(funrhs.expression(), r)
    //   , interface
    //   , In
    // );
    // intF = integralSurf(funrhsp, In, qTime) ;
    // intFnp = integralSurf(funrhs, In, qTime);

    double cc = 1./In.T.mesure()*6.;
    surfactant.addBilinear(
      innerProduct(cc*s, r)
      , interface
      , nbTime-1
      , In
    );
    surfactant.addLinear(
      innerProduct(u0.expression(), cc*r)
      , interface
      , 0
      , In
    );

    surfactant.addEdgeIntegral(
      innerProduct(1e-2*jump(grad(s)*n), jump(grad(r)*n)),
      In
    );
    // if stab interface, can have h in both stabilization
    // surfactant.addBilinear(
    //   innerProduct(1e-2*h_E*grad(s)*n, grad(r)*n)
    //   , interface
    //   , In
    // );

    #endif

    std::cout << " Time assembly  \t" << CPUtime() - tt0 << std::endl;

    surfactant.solve();

    KN_<double> dw(surfactant.rhs(SubArray(surfactant.nDoF, 0)));
    uh = dw;

    surfactant.saveSolution(uh);

    {
      myCout << "\n Features of the drop " << std::endl;
      Fun_h funuh(cutVh, uh);

      Rn sol2(cutVh.NbDoF(), 0.);
      Fun_h funsol(cutVh, sol2);
      sol2 += uh(SubArray(cutVh.NbDoF(), 0));
      double q_0 = integralSurf(funsol, 0);
      sol2 += uh(SubArray(cutVh.NbDoF(), cutVh.NbDoF()));
      double q_1 = integralSurf(funsol, 0, lastQuadTime);
      if(iter==0){ q0_0 = q_0;
                   q0_1 = q_1;
                   q0_1 = integralSurf(u0, 0);
                 }
      myCout << " q0 - qt_1 \t" << q_1-qp_1 << std::endl;
      myCout << " int F \t" << intF << std::endl;

      outputData << setprecision(10);
      outputData << GTime::current_time() << "\t"
                 << q_0 << "\t"
                 << q_1 << "\t"
                 << (q_1-qp_1) << "\t"
                 << fabs(q_1-qp_1) << "\t"
                 << q0_1 - q_1 << "\t"
                 << fabs(q0_1 - q_1) << "\t"
                 << std::endl;
                 qp_1 = q_1;
    }

    //  -----------------------------------------------------
    //                     PLOTTING
    //  -----------------------------------------------------
    if(MPIcf::IamMaster() && iter%frequencyPlottingTheSolution == 0 && writeVTKFiles){
      std::cout << " Plotting " << std::endl;

      if(saveSurfactantVTK) {
        KN_<double> sols(uh(SubArray(cutVh.NbDoF(), 0)));
        Fun_h sol(cutVh, sols);

        Paraview3 writer(cutVh, pathOutpuFigure+"surfactant_"+to_string(iterfig)+".vtk");
        writer.add(sol , "surfactant", 0, 1);
        writer.add(ls[0], "levelSet",  0, 1);
      }
      iterfig++;
    }

    iter++;
  }
}
