#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#ifdef USE_OMP
#include "/usr/local/opt/libomp/include/omp.h"
#endif

#include "num/print_container.hpp"
#include "num/util.hpp"
#include "problem/baseProblem.hpp"
#include "FESpace/expression.hpp"
#include "FESpace/paraview.hpp"
#include "problem/generalNorm.hpp"
#include "problem/projection.hpp"
#include "solver/solver.hpp"