/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#ifndef _CPP_PARALLEL_CFOMP_HPP_
#define _CPP_PARALLEL_CFOMP_HPP_

#include <thread>

#include "../cutFEMConfig.h"

inline int cutfem_get_max_threads() { return std::thread::hardware_concurrency(); }

#ifdef USE_OMP
#include <omp.h>
#else

inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_max_threads() { return 1; }

#endif
#endif