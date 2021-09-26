#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "finiteElement.hpp"
#include "projection.hpp"
#include "../num/gnuplot.hpp"


R fun_velocity(const R2 P, int ci) {
  R2 R(P.y-0.5,0); return R[ci];
}
R2 fparam(double t){ return R2(2.+1./3*cos(t+1./3), 0.5+ 1./3*sin(t+1./3));}



int main(int argc, char** argv )
{
  const int d = 2;

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef FunFEM<Mesh> Fun_h;

  MPIcf cfMPI(argc,argv);

  Mesh Th(40, 10, 0., 0., 4., 1.);

  Marker marker(Th, fparam, 0, 2*M_PI, 50);

  gnuplot::save(Th);
  // gnuplot::saveNode(marker);

  Lagrange2 FEvelocity(2);
  FESpace Uh(Th, FEvelocity);
  Fun_h velocity(Uh, fun_velocity);


  marker.make_interface();

  gnuplot::save(marker, "marker.dat");


  CutFESpace2 Wh(Uh, marker, {1, -1});
  gnuplot::save(Wh, 1);

  // gnuplot::saveNormal(marker);
  // double dt = 0.2;
  // for(int i=0;i<10;++i){
  //   marker.move(velocity, dt);
  //   gnuplot::save(marker, "node_t"+to_string(i)+".dat");
  //
  // }




};
