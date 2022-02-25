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
#include "baseCutProblem.hpp"
#include "spline.hpp"

#include "../num/gnuplot.hpp"


R fun_velocity(const R2 P, int ci) {
  R2 R(P.y,0); return R[ci];
}
R2 fparam(double t){ return R2(cos(t), sin(t));}



int main(int argc, char** argv )
{
  const int d = 2;

  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface2 Interface;
  typedef FunFEM<Mesh> Fun_h;

  MPIcf cfMPI(argc,argv);

  Mesh Th(7, 7, -2., -2., 4., 4.);

  Marker marker(Th, fparam, 0, 2*M_PI, 10);
  marker.make_interface();
  gnuplot::save(marker, "interface.dat");

  gnuplot::saveMarker(marker, "marker.dat");
  gnuplot::save(Th);
  // gnuplot::saveNode(marker);

  // Spline spline(marker.T, marker.X);
  // Spline spline(marker.T, marker.Y);
  int n = marker.vertices_.size();
  vector<double> x(n);
  vector<double> y(n);
  vector<double> t(n);
  double h = 1/(n-1);
  for(int i=0;i<n;++i) {
    x[i] = marker.vertices_[i].x;
    y[i] = marker.vertices_[i].y;
    // std::cout << P[i] << std::endl;
    t[i] = i*h;
  }
  // Spline2D spline(t, x, y);
  // spline.gnuplot(20, "spline2.dat");

  Spline splineX(t, x);
  splineX.gnuplot(10, "splineX.dat");
  Spline splineY(t, y);
  splineY.gnuplot(10, "splineY.dat");


  Lagrange2 FEvelocity(2);
  FESpace Uh(Th, FEvelocity);
  Fun_h velocity(Uh, fun_velocity);


  gnuplot::saveMarker(marker, "marker_t0.dat");


  CutFESpace2 Wh(Uh, marker, {1, -1});
  gnuplot::save(Wh, 1);

  // gnuplot::saveNormal(marker);
  double dt = 0.1;
  for(int i=0;i<20;++i){
    marker.move(velocity, dt);
    gnuplot::saveMarker(marker, "marker_t"+to_string(i+1)+".dat");
  }




};
