#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#ifdef USE_OMP
#include "/usr/local/opt/libomp/include/omp.h"
#endif

#include "../problem/baseProblem.hpp"
#include "../FESpace/paraview.hpp"
#include "../problem/generalNorm.hpp"