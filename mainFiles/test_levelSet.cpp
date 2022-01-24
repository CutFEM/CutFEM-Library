
#ifndef _TEST_LEVELSET_HPP
#define _TEST_LEVELSET_HPP

#include "../problem/levelSet.hpp"

#define EXAMPLE_1
#ifdef EXAMPLE_1
R fun_velocity(const R2 P,  const int i, const R t) {
  if(i==0) return -2.*sin(pi*P.x)*sin(pi*P.x)*cos(pi*P.y)*sin(pi*P.y);
  else return      2.*sin(pi*P.y)*sin(pi*P.y)*cos(pi*P.x)*sin(pi*P.x);
}
R fun_levelSet(const R2 P, const int i, const R t) {
  return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.75)*(P.y-0.75)) - 0.15;
}
#endif
#ifdef EXAMPLE_2
R fun_velocity(const R2 P,  const int i, const R t) {
  return (i==0);
}
R fun_levelSet(const R2 P, const int i, const R t) {
  return sqrt((P.x-t)*(P.x-t) + P.y*P.y) - 1;
}
#endif


int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 100;
  int ny = 100;
  double dt = 1./nx/4;
  double tid = 0.;
  int ifig = 0;
  Mesh2 Th(nx, ny, -3., -2., 6., 4.);
  Lagrange2 FEvelocity(2);
  FESpace2 Uh(Th, FEvelocity);
  Fun2_h velocity(Uh, fun_velocity, tid);
  Fun2_h velocity_p(Uh, fun_velocity, tid);

  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  LevelSet2 levelSet(Lh);
  Fun2_h ls(Lh, fun_levelSet, tid);

  if(MPIcf::IamMaster()){
    // Paraview2 writer(Lh, "testLevelSet"+to_string(ifig++)+".vtk");
    Paraview2 writer(Lh, "boxShear.vtk");

    writer.add(ls, "levelSet", 0, 1);
    writer.add(velocity, "velocity", 0, 2);
  }
return 0;
  for(int i=1;i<=500;++i){
    std::cout << "iteration => " << i << std::endl;

    levelSet.solve(ls, velocity_p, velocity, dt);
    ls.init(levelSet.rhs);

    tid += dt;
    velocity_p.v = velocity.v;
    velocity.init(Uh, fun_velocity, tid);


    if(MPIcf::IamMaster() && i%20==0){
      Paraview2 writer(Lh, "testLevelSet"+to_string(ifig++)+".vtk");
      writer.add(ls, "levelSet", 0, 1);
    }

  }

  return 0;

}
//       using namespace NumericLevelSet2D ;
//       const int numberOfIteration = 1;
//       vector<double> hv(numberOfIteration);
//       vector<double> error(numberOfIteration);
//       vector<double> convergence(numberOfIteration);
//
//       for(int i=0; i<numberOfIteration;++i) {
//
//         int nx = 51*pow(2,i);
//         Mesh2 Th(nx,nx,-1.70001,-1.70001,3.40002,3.40002);
//         FESpace2 Vh(Th, Pk);
//         FESpace2 VelVh(Th, DataFE<Mesh2>::P2);
//
//         const double h = Th[0].lenEdge(0);
//         const double dt = h/4;
//
//         GTime::time_step = dt;
//         GTime::total_number_iteration = 20;//int(0.25/dt);//int(h / dt);//20;
//         GTime::t0 = 0;
//         GTime::current_iteration = 0;
//
//
//         std::cout << "-------------------------------------------------------------" << st
// d::endl;
//         std::cout << " iteration : \t " << i << std::endl;
//         std::cout << " mesh size :\t" << h << std::endl;
//         std::cout << " nx :\t" << nx << std::endl;
//         std::cout << " time step :\t" << GTime::total_number_iteration << std::endl;
//
//
//         KN<double> velocityFieldp, velocityField;
//         FunVec2 vel, velp;
//         interpolate(VelVh,
//                     velocityFieldp,
//                     fun_velocity_field,
//                     GTime::current_time());
//
//         LevelSet2 levelSet(Vh, fun_levelSet);
//         levelSet.dt = dt;
//
//         // Rn ls;
//         // interpolate(Vh,ls,fun_levelSet);
//         // const FEMFunction<FESpace2> mySet(Vh, ls);
//         // mySet.toVTK("levelset0.vtk", "levelSet");
//
//
//         int iter = 0;
//         while( iter < GTime::total_number_iteration ) {
//           GTime::current_iteration = iter;
//
//           velp.init(VelVh, velocityFieldp);
//           interpolate(VelVh,
//                       velocityField,
//                       fun_velocity_field,
//                       GTime::current_time());
//
//           vel.init(VelVh, velocityField);
//           levelSet.solve(velp, vel);
//           velocityField.swap(velocityFieldp);
//
//
//           // // extract level set function
//           // const FEMFunction<FESpace2> mySet2(Vh, levelSet.rhs);
//           // mySet2.toVTK("levelset"+to_string(iter+1)+".vtk", "levelSet");
//
//           iter++;
//         }
//
//         Rn eh;
//         interpolate(Vh,
//                     eh,
//                     fun_solution,
//                     GTime::current_time());
//         eh -= levelSet.rhs;
//         Integral<Mesh2> integration(Vh);
//         std::cout << integration.L2norm(eh) << std::endl;
//
//       }
//
//       return 0;
//     }
//

#endif;
