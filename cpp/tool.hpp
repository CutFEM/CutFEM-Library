#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#ifdef USE_OMP
#include "/usr/local/opt/libomp/include/omp.h"
#endif

#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <concepts>

#include "concept/function.hpp"
#include "common/global.hpp"
#include "common/Mesh3dn.hpp"
#include "num/print_container.hpp"
#include "num/util.hpp"
#include "problem/baseProblem.hpp"
#include "FESpace/expression.hpp"
#include "FESpace/paraview.hpp"
#include "problem/generalNorm.hpp"
#include "problem/projection.hpp"
#include "solver/solver.hpp"