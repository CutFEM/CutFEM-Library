#ifndef SOLVER_HPP_
#define SOLVER_HPP_
#include <cassert>
#include "cutFEMConfig.h"
#include "../util/util.hpp"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#else
#include "../util/cputime.h"
#endif
// #include "../parallel/partitioner.hpp"

namespace solver {

  static void test_Prep() {
    #ifdef USE_MPI
    std::cout << " Seen USE_MPI" << std::endl;
    #endif
    #ifdef USE_LAPACK
    std::cout << " Seen USE_LAPACK" << std::endl;
    #endif
    #ifdef USE_UMFPACK
    std::cout << " Seen USE_UMFPACK" << std::endl;
    #endif

  }

    void umfpack(std::map<std::pair<int,int>,R> &, Rn &);
    void LAPACK(Rnm & a, Rn & b);
//   void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const int) ;
//   void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const int, const MatrixOrdering*) ;
  // void MUMPS(std::map<std::pair<int,int>,R> &, Rn &, const MatrixOrdering*) ;
  // void MUMPS(std::map<std::pair<int,int>,R> &, Rn &) ;

//   void PETSC_Sequential(std::map<std::pair<int,int>,R> &, Rn &);
// //   void PETSC_Parallel(std::map<std::pair<int,int>,R> &, Rn &);
//   void PETSC(std::map<std::pair<int,int>,R> &, Rn &);

}
// namespace toolPETSC {
//   inline int whatProc(const int i, const int * Iend, const int nproc);
//   inline void map2AIJ(std::map<std::pair<int,int>,R> & Amap,
// 		      KN<int>& Iv, KN<int>& Jv, KN<double>& a);
//   int getLocalSize(int sizeGlob, int&Istart, int& Iend);
//   void transform2PETSCformat(std::map<std::pair<int,int>,R> & Amap,
// 			     int Istart, int Iend);

// }


class Solver {


#ifdef USE_MPI
double get_Time() const {return MPIcf::Wtime();}
#else
double get_Time() const {return CPUtime();}
#endif

public :
int verbose = 0;
//   std::string reordering = "none";
  //std::string solver = "umfpack";

  void solve(std::map<std::pair<int,int>,R> & A, Rn & b) ;


};








#endif
