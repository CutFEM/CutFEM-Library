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


typedef std::map<std::pair<int,int>,R> MatMap;


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef ExpressionFunFEM<Mesh2> Expr;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = MPIcf::Wtime();

  int nx = 2;
  int ny = 2;
  Mesh Th(nx, ny, 0., 0., 1, 1.);
  const Element& K(Th[0]);   // The reference element


  };
