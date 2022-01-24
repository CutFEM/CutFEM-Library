#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <omp.h>
#include "../util/cputime.h"
#include "cfmpi.hpp"
#include "baseCutProblem.hpp"
#include "levelSet.hpp"
#include "../num/gnuplot.hpp"
#include "../FESpace/limiter.hpp"

#define PROJECTION_TEST

#ifdef PROJECTION_TEST
R fun_solution(const R2 P, int elementComp) {
  return sin(pi*(P.x+P.y));
}


int main(int argc, char** argv ) {

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  // DEFINITION OF THE MESH
  // =====================================================
  Mesh2 Th(50, 50, 0., 0., 1, 1.);
  R2 P( 0.5001, 0.500005);
  gnuplot::save(Th);
  int k = find_triangle_contenant_p(Th, P);

  std::cout << k << std::endl;


  // double err0;
  // Mesh2 Th(2, 2, 0., 0., 1, 1.);
  // const typename Mesh2::Element& K(Th[0]);   // The reference element
  //
  // Filter_Box filter_box( R2(0.15001,0.2001), R2(0.5001, 0.7001));
  //
  // {
  //   std::ofstream plot;
  //   plot.open("element.dat", std::ofstream::out);
  //   for(int i=0;i<3;++i) {
  //     plot << K[i] << std::endl;
  //   }
  //   plot << K[0] << std::endl;
  //   plot << std::endl;
  //   plot << std::endl;
  //
  //   plot.close();
  //   plot.open("box.dat", std::ofstream::out);
  //   int ni = 10;
  //   double hi = 1./ni;
  //   for(int i = 0; i<=ni; ++i){
  //     for(int j = 0; j<=ni-i; ++j){
  //       R2 P(i*hi, j*hi);
  //       if(filter_box.contain(P)) {
  //         plot << P << std::endl;
  //       }
  //     }
  //   }
  //   plot.close();
  // }
  //
  // if(filter_box.intersect_with(K)) {
  //   std::cout << " Triangle intersect with box" << std::endl;
  // }
  // filter_box.cut_triangle(K);
  // filter_box.print();


  // int nx = 10;
  // int ny = 10;
  // Mesh2 Th(2, 2, 0., 0., 1, 1.);
  // for(int i=0;i<5;++i) {
  //
  //   Mesh2 Th(nx, ny, -1., -1., 2., 2.);
  //   FESpace2 Vh(Th, DataFE<Mesh2>::P0);
  //   Fun2_h Un(Vh, fun_solution);
  //
  //
  //   double errU = L2norm(Un,fun_solution,0,1);
  //   double cr =(i==0)? 0 :log(errU/err0)/log(1./2);
  //   err0 = errU;
  //   std::cout << 2./nx << "\t" << errU << "\t" << cr << std::endl;
  //
  //   nx*= 2;
  //   ny*= 2;
  // }



}
#endif
