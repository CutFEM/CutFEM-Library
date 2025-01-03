

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <mutex>
#include <thread>
#include <concepts>
#include <algorithm>
#include <numeric>
#include <numbers>

#include "concept/function.hpp"

#include "parallel/cfmpi.hpp"
#include "parallel/cfomp.hpp"

#include "common/global.hpp"
#include "common/Mesh3dn.hpp"
#include "common/time_interface.hpp"
#include "common/logger.hpp"

#include "num/DA.hpp"
#include "num/print_container.hpp"
#include "num/util.hpp"
#include "num/matlab.hpp"

#include "FESpace/expression.hpp"
#include "FESpace/integrationFunFEM.hpp"
#include "FESpace/paraview.hpp"
#include "FESpace/funfem_util.hpp"

#include "solver/solver.hpp"

#include "problem/baseProblem.hpp"
#include "problem/generalNorm.hpp"
#include "problem/projection.hpp"
#include "problem/solver_advection.hpp"
#include "problem/solver_curvature.hpp"
#include "problem/solver_stokes.hpp"
#include "problem/time_scheme.hpp"

