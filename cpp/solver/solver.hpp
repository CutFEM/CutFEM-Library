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
#ifndef SOLVER_HPP_
#define SOLVER_HPP_
#include <cassert>
#include "cutFEMConfig.h"
#include "../num/util.hpp"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

struct ProblemOption {
   int order_space_element_quadrature_ = 5;
   int order_space_bord_quadrature_    = 5;
   int order_time_quadrature_          = 1;
   std::string solver_name_            = "mumps";
   bool clear_matrix_                  = true;
   int verbose_                        = 0;
};

namespace solver {

void umfpack(std::map<std::pair<int, int>, R> &, Rn &, bool);
void LAPACK(Rnm &a, Rn &b);
} // namespace solver
class Solver {

#ifdef USE_MPI
   double get_Time() const { return MPIcf::Wtime(); }
#else
   double get_Time() const { return CPUtime(); }
#endif

 public:
   int verbose_             = 0;
   bool clearMatrix_        = true;
   //   std::string reordering = "none";
   std::string solver_name_ = "default";

   Solver(const ProblemOption &option) {
      clearMatrix_ = option.clear_matrix_;
      solver_name_ = option.solver_name_;
      verbose_     = option.verbose_;
   }

   void solve(std::map<std::pair<int, int>, R> &A, Rn &b);
};

#endif

//   void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const int) ;
//   void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const int, const
//   MatrixOrdering*) ;
// void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const MatrixOrdering*) ;
// void MUMPS(std::map<std::pair<int,int>,R> &, Rn &) ;

//   void PETSC_Sequential(std::map<std::pair<int,int>,R> &, Rn &);
// //   void PETSC_Parallel(std::map<std::pair<int,int>,R> &, Rn &);
//   void PETSC(std::map<std::pair<int,int>,R> &, Rn &);

// namespace solver
// namespace toolPETSC {
//   inline int whatProc(const int i, const int * Iend, const int nproc);
//   inline void map2AIJ(std::map<std::pair<int,int>,R> & Amap,
// 		      KN<int>& Iv, KN<int>& Jv, KN<double>& a);
//   int getLocalSize(int sizeGlob, int&Istart, int& Iend);
//   void transform2PETSCformat(std::map<std::pair<int,int>,R> & Amap,
// 			     int Istart, int Iend);

// }