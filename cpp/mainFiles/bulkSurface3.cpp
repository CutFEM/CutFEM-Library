#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>

#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
#include "projection.hpp"
#include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"

double r0 = 0.5;
R fun_levelSet(const R3 P, const int i) { return sqrt((P.x-0.1)*(P.x-0.1) + (P.y-0.0)*(P.y-0.0) + (P.z-0.0)*(P.z-0.0)) - r0;}
R fun_velocity(const R3 P, const int i){
  if(i == 0)     return -0.5*(1+cos(M_PI*P.x))*sin(M_PI*P.y);
  else if(i==1)  return  0.5*(1+cos(M_PI*P.y))*sin(M_PI*P.x);
  else return 0;
}
R fun_init_surface(const R3 P, int i) { return 0.;}
R fun_w(double r) {
  return 0.5*(1-cos((r-r0)*M_PI/(0.5*0.3)));
}
R fun_init_bulk(const R3 P, int i) {
  double r = sqrt((P.x-0.1)*(P.x-0.1) + P.y*P.y + P.z*P.z);
  if(r > 1.5*r0) return 0.5*(1-P.x*P.x)*(1-P.x*P.x);
  else if(r0 < r && r <= 1.5*r0) return 0.5*(1-P.x*P.x)*(1-P.x*P.x)*fun_w(r);
  else return 0.;}



int main(int argc, char** argv ){

  const int d = 3;
  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> Space;
  typedef CutFESpace<Mesh> CutSpace;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;



  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 20;
  int ny = 20;
  int nz = 20;

  Mesh Kh(nx, ny, nz, -1., -1., -1., 2., 2., 2.);
  Space Vh(Kh, DataFE<Mesh>::P1);

  double hi = 2./(nx-1);
  int divisionMeshsize = 2;
  double Tend = 0.5;
  double dT = hi/divisionMeshsize;
  GTime::total_number_iteration = Tend/dT;
  dT = Tend / GTime::total_number_iteration;
  GTime::time_step = dT;

  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();
  const Uint lastQuadTime = nbTime-1;


  // Set parameters for paraview PLOTTING
  const bool writeVTKFiles = true;
  const bool saveSurfactantVTK = true;
  const int frequencyPlottingTheSolution = 1;


  // // Space for the velocity field
  // // ----------------------------------------------
  Lagrange3<Mesh> FEvelocity(1);
  Space VelVh(Kh, FEvelocity);
  Fun_h vel(VelVh, fun_velocity);

  // levelSet stuff
  Space Lh  (Kh, DataFE<Mesh>::P1);
  double dt_levelSet = dT/(nbTime-1);
  vector<Fun_h> ls_k(nbTime), ls(nbTime);
  for(int i=0;i<nbTime;++i) ls[i].init(Lh, fun_levelSet);



  // Declaration of the interface
  TimeInterface<Mesh> interface(qTime);

  // Stokes problem
  // ----------------------------------------------
  CutFEM<Mesh> bulkSurface(qTime);

  const R epsilon_surfactant = 1.;
  const R Bb = 1.;
  const R Bs = 1.;
  const R Bbs = 1.;
  const R Kb = 0.01;
  const R Ks = 1.;
  const R Kbs = 0.;


  int iter = 0, iterfig = 0;
  while( iter < GTime::total_number_iteration ) {
    double tbegin_iter = CPUtime();
    GTime::current_iteration = iter;
    const TimeSlab& In(Ih[iter]);

    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ------------------------------------------------------------- "<< std::endl;
    std::cout << " ITERATION \t : \t" << iter+1 << " / " << GTime::total_number_iteration << std::endl;
    std::cout << " TIME      \t : \t" << GTime::current_time() << std::endl;


    double t_ls0 = CPUtime();
    ls.begin()->swap(ls[lastQuadTime]);
    for(int i=0;i<nbTime;++i) {
      interface.init(i,Kh,ls[i]);

      if(i<lastQuadTime) {
        LevelSet::move(ls[i], vel, vel, dt_levelSet, ls[i+1]);
      }
    }


    ActiveMesh<Mesh> Kh1(Kh);
    Kh1.truncate(interface, -1);
    ActiveMesh<Mesh> Kh0(Kh);
    Kh0.createSurfaceMesh(interface);
    CutSpace Wh(Kh1, Vh);
    CutSpace Sh(Kh0, Vh);



    /*
    PROBLEM DEFINITION
    */
    Normal n;
    FunTest ub(Wh,1), vb(Wh,1);           // bulk function
    FunTest us(Sh,1), vs(Sh,1);     // surface function

    // Connect spaces to problem
    bulkSurface.initSpace(Wh, In);
    bulkSurface.add(Sh, In);

    // initialize the solution
    int idx_s0 = Wh.NbDoF()*In.NbDoF();
    Rn data_u0(bulkSurface.get_nb_dof(), 0.);
    Rn_ data_b0(data_u0(SubArray(Wh.NbDoF()   ,0)));
    Rn_ data_s0(data_u0(SubArray(Sh.NbDoF(),idx_s0)));

    if(iter == 0) {
      interpolate(Wh, data_b0, fun_init_bulk);
      interpolate(Sh, data_s0, fun_init_surface);
    } else {
      bulkSurface.initialSolution(data_u0);
    }
    Fun_h b0h(Wh, data_b0);
    Fun_h s0h(Sh, data_s0);

    // set uh for newton
    Rn data_uh(data_u0);
    Rn_ data_bh(data_uh(SubArray(Wh.NbDoF()*In.NbDoF(),0)));
    Rn_ data_sh(data_uh(SubArray(Sh.NbDoF()*In.NbDoF(),idx_s0)));
    Fun_h bh(Wh, In, data_bh);
    Fun_h sh(Sh, In, data_sh);

    double t0 = CPUtime();

    bulkSurface.addBilinear(
      innerProduct(Bb*dt(ub), vb)
      +innerProduct((vel.expression()*grad(ub)), vb)*Bb
      +innerProduct(grad(ub), grad(vb))*Bb*Kb
      , Kh1
      , In
    );

    bulkSurface.addBilinear(
      innerProduct(Bs*dt(us), vs)
      +innerProduct((vel.expression()*grad(us)), vs)*Bs
      +innerProduct(divS(vel)*us, vs)*Bs
      +innerProduct(gradS(us), gradS(vs))*Ks*Bs
      +innerProduct(Bb*ub-Bs*us, Bb*vb-Bs*vs)
      , interface
      , In
    );
    // R h3 = pow(meshSize,3);
    // const & h_E(Parameter::meas);
    bulkSurface.addFaceStabilization(
      innerProduct(1e-2*hi*jump(grad(ub)*n), jump(grad(vb)*n))
      ,Kh1
      ,In
    );
    bulkSurface.addBilinear(
      innerProduct(1e-2*hi*jump(grad(us)*n), jump(grad(vs)*n))
      , Kh0
      , innerFacet
      , In
    );

    // impose initial condition
    bulkSurface.addBilinear(
      innerProduct(ub,vb)*Bb
      , Kh1
    );
    bulkSurface.addBilinear(
      innerProduct(us, vs)*Bs
      , *interface(0)
      , In
      , 0
    );
    std::cout << " Time matrix assembly " << CPUtime()-t0 << std::endl;
    t0 = CPUtime();
    /*
    CONSTRUCTION LINEAR PART
    */
    // tt0 = CPUtime();
    // impose initial condition
    bulkSurface.addLinear(
      innerProduct(b0h.expression(), vb)*Bb
      , Kh1
    );
    bulkSurface.addLinear (
      innerProduct(s0h.expression(), vs)*Bs
      , *interface(0)
      , In , 0
    );
    std::cout << " Time RHS assembly " << CPUtime()-t0 << std::endl;

    bulkSurface.solve();

    Rn_ dw(bulkSurface.rhs_(SubArray(bulkSurface.get_nb_dof(), 0)));
    data_uh = dw;

    bulkSurface.saveSolution(data_uh);
    bulkSurface.cleanMatrix();
// //         break;
// //       }
//     // }
//     // // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
//     // {
//     //   Fun2_h funX(Wh, fun_x);
//     //   Fun2_h fun1(Wh, fun_1);
//     //   R tt0 = CPUtime();
//     //   double areaBubble   = integral(fun1, 0, 1) ;
//     //   double centerOfMass = integral(funX, 1, 1) / areaBubble ;
//     //
//     //   double Pb = integralSurf(fun1, 1);
//     //   double ra = sqrt(areaBubble / M_PI);
//     //   double Pa = 2*M_PI*ra;
//     //   double circularity = Pa / Pb;
//     //   double riseVelocity = integral(uh, 1, 1) / areaBubble;
//     //
//     //   double q = integralSurf(us, 0);
//     //
//     //   std::cout << "\n Features of the drop " << std::endl;
//     //   std::cout << " Time                   ->    " << GTime::current_time()+dT  << std::endl;
//     //   std::cout << " Center Of Mass         ->    " << centerOfMass << std::endl;
//     //   std::cout << " Circularity            ->    " << circularity << std::endl;
//     //   std::cout << " Rise velocity          ->    " << riseVelocity << std::endl;
//     //   std::cout << " Surfactant quantity    ->    " << q << std::endl;
//     //
//     //
//     //   outputData << GTime::current_time()+dT << "\t"
//     //              << centerOfMass << "\t"
//     //              << circularity << "\t"
//     //              << riseVelocity << "\t"
//     //              << areaBubble <<  "\t"
//     //              << q << std::endl;
//     //
//     // }
//     //
//     //
//     // // std::cout << " Set Velocity " << std::endl;
//     // for(int i=0;i<nbTime;++i) {
//     //   Fun2_h sol(Wh,In,data_uh);
//     //   set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
//     // }
//     //
//     // matlab::Export(bulkSurface.mat_, "matA.dat");
//     // return 0;
//  -----------------------------------------------------
//                     PLOTTING
//  -----------------------------------------------------
if(iter%frequencyPlottingTheSolution == 0 && writeVTKFiles && MPIcf::IamMaster()){
  std::cout << " Plotting " << std::endl;

  if(saveSurfactantVTK) {
    Paraview<Mesh> writer(Kh1, "bulk3_"+to_string(iterfig)+".vtk");
    writer.add(bh, "bulk", 0, 1);

    Paraview<Mesh> writerS(Kh0, "surface3_"+to_string(iterfig)+".vtk");
    writerS.add(sh, "surface", 0, 1);
    writerS.add(ls[0], "levelSet", 0, 1);

    Paraview<Mesh> writerLS(Kh, "levelSet_"+to_string(iterfig)+".vtk");
    writerLS.add(ls[0], "levelSet", 0, 1);
  }

  iterfig++;
}
std::cout << " Timer iteration computation \t" << CPUtime() - tbegin_iter << std::endl;
    iter += 1;
  }

}
