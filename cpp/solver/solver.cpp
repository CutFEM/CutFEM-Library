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
#include <fstream>
#include <iostream>
#include "../num/util.hpp"
#include "solver.hpp"
#include "../common/SparseMatMap.hpp"
#ifdef USE_UMFPACK
#include "umfpack.h"
#endif
#ifdef USE_LAPACK
#include "lapacke.h"
#endif

// #include "petsc.h"
#ifdef USE_MUMPS
#include "solverMumps.hpp"
#endif

void Solver::solve(std::map<std::pair<int, int>, R> &A, Rn &b) {

   double tsolver = this->get_Time();
   if (solver_name_ == "mumps") {
#ifdef USE_MUMPS
      MUMPS(*this, A, b);
#else
#ifdef USE_UMFPACK
      solver::umfpack(A, b, clearMatrix_);
#else
      // assert(0);
#endif
#endif
   } else if (solver_name_ == "umfpack") {
      if (verbose_ > 1)
         std::cout << " solve using umfpack" << std::endl;
#ifdef USE_UMFPACK
      solver::umfpack(A, b, clearMatrix_);
#else
         // assert(0);
#endif
   } else {
      // assert(0);
   }
   // solver::umfpack(A,b);

   tsolver = this->get_Time() - tsolver;

   if (this->verbose_ > 0)
      std::cout << " Real Time Solver \t \t " << tsolver << std::endl;
}

namespace solver {
#ifdef USE_UMFPACK
void umfpack(std::map<std::pair<int, int>, R> &Amap, Rn &b, bool clearMatrix) {
   const int n = b.size();
   KN<double> x(n);
   SparseMatrixRC<double> A(n, n, Amap);
   if (clearMatrix)
      Amap.clear();
   void *Symbolic, *Numeric;
   (void)umfpack_di_symbolic(n, n, A.p, A.j, A.a, &Symbolic, 0, 0);
   (void)umfpack_di_numeric(A.p, A.j, A.a, Symbolic, &Numeric, 0, 0);
   umfpack_di_free_symbolic(&Symbolic);
   (void)umfpack_di_solve(UMFPACK_At, A.p, A.j, A.a, x, b, Numeric, 0, 0);
   umfpack_di_free_numeric(&Numeric);
   b = x;
}
#endif

#ifdef USE_LAPACK
void LAPACK(Rnm &a, Rn &b) {
   lapack_int n    = a.N();
   lapack_int m    = a.M();
   lapack_int lda  = n;
   lapack_int ldb  = 1;
   lapack_int nrhs = 1;
   lapack_int info =
       LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, 1, a, lda, b, ldb);
}
#endif

} // namespace solver

//   void MUMPS(std::map<std::pair<int,int>,R> & Amap, Rn & rhs, const int nloc)
//   {

//     NoOrdering mapp;
//     MUMPS(Amap, rhs, nloc,&mapp);
//   }

//   void MUMPS(std::map<std::pair<int,int>,R> & Amap, Rn & rhs,
// 	     const int nloc, const MatrixOrdering * mapp) {

//     DMUMPS_STRUC_C mumps_par;
//     // Definition of the MPI communicator for MUMPS
//     //------------------------------------------------------
//     mumps_par.comm_fortran=USE_COMM_WORLD;
//     // Symmetry of the matrix
//     //------------------------------------------------------
//     mumps_par.sym = 0;
//     // Type of parallelism (par = 1 : host working - par = 0 : host not
//     working)
//     // -----------------------------------------------------
//     mumps_par.par = 1;

//     // Initialization of one instance of the package
//     //------------------------------------------------------
//     mumps_par.job = JOB_INIT;
//     dmumps_c(&mumps_par);

//     // Standard output (no output)
//     //------------------------------------------------------
//     mumps_par.ICNTL(1) = -1;
//     mumps_par.ICNTL(2) = -1;
//     mumps_par.ICNTL(3) = -1;
//     mumps_par.ICNTL(4) = -1;

//     // Definition of the matrix, distributed assembled format
//     //-------------------------------------------------------
//     mumps_par.ICNTL(5)  = 0;
//     mumps_par.ICNTL(18) = 3;

//     // Pre-treatment of the matrix (permutation, ordering)
//     // 7 : automatic choice
//     //------------------------------------------------------
//     mumps_par.ICNTL(28) = 1;
//     mumps_par.ICNTL(6) = 7;
//     mumps_par.ICNTL(7) = 3;
//     mumps_par.ICNTL(8) = 77;

//     // Increase MAXIS (cf. doc MUMPS) for extra fill-in
//     //-------------------------------------------------------
//     // mumps_par.ICNTL(14) = 30;

//     // Format of the right hand side
//     //-------------------------------------------------------
//     mumps_par.ICNTL(20) = 0;

//     // Evaluating the number of non-zero elements in the matrix
//     //--------------------------------------------------------
//     Uint nz_glob = 0, nz_loc = Amap.size();
//     MPIcf::AllReduce(nz_loc, nz_glob,MPI_SUM);

//     const int ndof_glob = rhs.size();
//     const int ndof_loc  = nloc;

//     if(MPIcf::IamMaster()) {
//       mumps_par.n  = ndof_glob;
//       mumps_par.nz = nz_glob;
//     }
//     mumps_par.nz_loc = nz_loc;

//     // Construction of the local matrices for MUMPS
//     //-------------------------------------------------------
//     KN<int> IRN_loc(mumps_par.nz_loc), JCN_loc(mumps_par.nz_loc);
//     Rn A_loc(mumps_par.nz_loc);
//     int k = 0;
//     for (std::map<std::pair<int,int>,R>::const_iterator q=Amap.begin();
// 	 q != Amap.end(); ++q, ++k)
//       {
// 	// !!!!!!     FORTRAN NUMBERING    !!!!!!!!!!
// 	IRN_loc(k) = mapp->iperm(q->first.first)  + 1;
// 	JCN_loc(k) = mapp->iperm(q->first.second) + 1;
// 	A_loc(k)   = q->second;
//       }
//     Amap.clear();

//     mumps_par.irn_loc = IRN_loc;
//     mumps_par.jcn_loc = JCN_loc;
//     mumps_par.a_loc = A_loc;

//     // Construction of the right hand side
//     // The full rhs is saved on the root.
//     // the rhs will also be used to save the solution
//     //-------------------------------------------------------
//     Rn rhs_sum(ndof_glob*(MPIcf::IamMaster()));                     // alloc
//     mem for global rhs MPIcf::Reduce(rhs, rhs_sum , MPI_SUM,
//     MPIcf::Master());

//     Rn rhsMapp;

//     if(MPIcf::IamMaster()) {
//       mapp->reorder(rhs_sum, rhsMapp);
//       mumps_par.rhs = rhsMapp;
//     }

//     R ierr =0;
//     // Analysis phase
//     //-------------------------------------------------------
//     //  macro s.t. indices match documentation
// #define INFO(I) info[(I)-1]
//     mumps_par.job = JOB_ANALYSIS;
//     R ta0 = MPI_Wtime();
//     dmumps_c(&mumps_par);
//     R ta1 = MPI_Wtime();
//     ierr = mumps_par.INFO(1);
//     if(ierr != 0) {
//       std::cout << " Error in analysis phase of MUMPS : ierr = " << ierr <<
//       std::endl; std::cout <<mumps_par.INFO(2)  << std::endl; std::cout
//       <<mumps_par.ICNTL(2) << std::endl; MPIcf::Barrier(); assert(0);
//     }

//     // Factorization phase
//     //-------------------------------------------------------
//     mumps_par.job = JOB_FACTORIZATION;
//     R tf0 = MPI_Wtime();
//     dmumps_c(&mumps_par);
//     R tf1 = MPI_Wtime();
//     ierr = mumps_par.INFO(1);

//     if(ierr != 0) {
//       std::cout << " Error in factorization phase of MUMPS : ierr = " << ierr
//       << std::endl; MPIcf::Barrier(); assert(0);
//     }

//     // Solving phase
//     //--------------------------------------------------------
//     R ts0 = MPI_Wtime();
//     mumps_par.job = JOB_SOLVE;
//     dmumps_c(&mumps_par);
//     R ts1 = MPI_Wtime();

//     ierr = mumps_par.INFO(1);
//     if(ierr != 0) {
//       std::cout << " Error in solving phase of MUMPS : ierr = " << ierr <<
//       std::endl; MPIcf::Barrier(); assert(0);
//     }

//     // Distribution of the DOFs of the solution on each processor
//     // (the right-hand side, stored on the host, stores the solution)
//     //--------------------------------------------------------
//     if(MPIcf::IamMaster())  {
//       mapp->inverseMapp(rhsMapp, rhs);
//     }
//     MPIcf::Bcast(rhs, MPIcf::Master());

//     // Terminate instance
//     //--------------------------------------------------------
//     mumps_par.job=JOB_END;
//     dmumps_c(&mumps_par);

//     // Get information about the resolution of the linear system
//     // 9 : Size of the space used to store the factor matrices.
//     // 16 : total size in million of bits of data allocated
//     //      during the factorization
//     //----------------------------------------------------------
//      R szlumn, szlumx, szlu = mumps_par.INFO(9);
//     MPIcf::AllReduce(szlu, szlumn, MPI_MIN);
//     MPIcf::AllReduce(szlu, szlumx, MPI_MAX);
//     R szwkmn, szwkmx, szwk = mumps_par.INFO(16);
//     MPIcf::AllReduce(szwk, szwkmn, MPI_MIN);
//     MPIcf::AllReduce(szwk, szwkmx, MPI_MAX);
//     if(MPIcf::IamMaster()) {
//     szlumn = int((szlumn*8.0)/(1024.0*1024.0));
//     szlumx = int((szlumx*8.0)/(1024.0*1024.0));

//     R ratio = ((R)(mumps_par.nz/mumps_par.n)/mumps_par.n)*100.0;
//     std::cout << " --------------------------------------------------------
//     \n"; std::cout << "                MUMPS DIRECT SOLVER               " <<
//     std::endl; std::cout << "
//     -------------------------------------------------------- \n"; std::cout
//     <<" STATISTICS OF THE GLOBAL MATRIX " << std::endl; std::cout << " Matrix
//     order                         " << mumps_par.n << std::endl; std::cout <<
//     " Number of non-zero entries           " << mumps_par.nz << std::endl;
//     std::cout << " Fill-in ratio percentage             " << ratio <<
//     std::endl; std::cout << "\n STATISTICS OF THE LU FACTORIZATION " <<
//     std::endl; std::cout << " Number of entries in the factors     " <<
//     mumps_par.infog[19] << std::endl; std::cout << "\n Storage of the factors
//     " << std::endl; std::cout << " Minimum memory                       " <<
//     szlumn << std::endl; std::cout << " Maximum memory " << szlumx <<
//     std::endl; std::cout << "\n Working memory for factorization   " <<
//     std::endl; std::cout << " Minimum memory                       " <<
//     szwkmn << std::endl; std::cout << " Maximum memory " << szwkmx <<
//     std::endl; std::cout << std::endl; std::cout << " Time of analysis phase
//     " << ta1 - ta0 << std::endl; std::cout << " Time of factorization phase
//     " << tf1 - tf0 << std::endl; std::cout << " Time for solving " << ts1 -
//     ts0 << std::endl << std::endl;
//   }

//   }

//   void MUMPS(std::map<std::pair<int,int>,R> & Amap, Rn & rhs) {

//     NoOrdering mapp;
//     MUMPS(Amap, rhs, &mapp);

//   }

//   void MUMPS(std::map<std::pair<int,int>,R> & Amap, Rn & rhs, const
//   MatrixOrdering * mapp) {

//     // MPIcf cfMPI;

//     // MPI_Init(nullptr, nullptr);                    // initialize MPI

//     DMUMPS_STRUC_C mumps_par;
//     // Definition of the MPI communicator for MUMPS
//     //------------------------------------------------------
//     mumps_par.comm_fortran=USE_COMM_WORLD;
//     // Symmetry of the matrix
//     //------------------------------------------------------
//     mumps_par.sym = 0;
//     // Type of parallelism (par = 1 : host working - par = 0 : host not
//     working)
//     // -----------------------------------------------------
//     mumps_par.par = 1;

//     // Initialization of one instance of the package
//     //------------------------------------------------------
//     mumps_par.job = JOB_INIT;
//     dmumps_c(&mumps_par);

//     // Standard output (no output)
//     //------------------------------------------------------
//     mumps_par.ICNTL(1) = -1;
//     mumps_par.ICNTL(2) = -1;
//     mumps_par.ICNTL(3) = -1;
//     mumps_par.ICNTL(4) = -1;

//     // Definition of the matrix, distributed assembled format
//     //-------------------------------------------------------
//     mumps_par.ICNTL(5)  = 0;
//     mumps_par.ICNTL(18) = 0;//3;

//     // Pre-treatment of the matrix (permutation, ordering)
//     // 7 : automatic choice
//     //------------------------------------------------------
//     mumps_par.ICNTL(28) = 1;
//     mumps_par.ICNTL(6) = 7;
//     mumps_par.ICNTL(7) = 3;
//     mumps_par.ICNTL(8) = 77;

//     // Increase MAXIS (cf. doc MUMPS) for extra fill-in
//     //-------------------------------------------------------
//     // mumps_par.ICNTL(14) = 30;

//     // Format of the right hand side
//     //-------------------------------------------------------
//     mumps_par.ICNTL(20) = 0;

//     // Evaluating the number of non-zero elements in the matrix
//     //--------------------------------------------------------
//     Uint nz_glob = Amap.size(), nz_loc = Amap.size();
//     const int ndof_glob = rhs.size();
//     const int ndof_loc  = rhs.size();

//     mumps_par.n  = ndof_glob;
//     mumps_par.nz = nz_glob;
//     // mumps_par.nz_loc = nz_loc;

//     // Construction of the local matrices for MUMPS
//     //-------------------------------------------------------
//     // KN<int> IRN_loc(mumps_par.nz_loc), JCN_loc(mumps_par.nz_loc);
//     // Rn A_loc(mumps_par.nz_loc);
//     KN<int> IRN_loc(mumps_par.nz), JCN_loc(mumps_par.nz);
//     Rn A_loc(mumps_par.nz);
//     int k = 0;
//     for (std::map<std::pair<int,int>,R>::const_iterator q=Amap.begin();
// 	 q != Amap.end(); ++q, ++k)
//       {
// 	// !!!!!!     FORTRAN NUMBERING    !!!!!!!!!!
// 	IRN_loc(k) = mapp->iperm(q->first.first)  + 1;
// 	JCN_loc(k) = mapp->iperm(q->first.second) + 1;
// 	A_loc(k)   = q->second;
//       }
//     Amap.clear();

//     // mumps_par.irn_loc = IRN_loc;
//     // mumps_par.jcn_loc = JCN_loc;
//     // mumps_par.a_loc = A_loc;

//     mumps_par.irn = IRN_loc;
//     mumps_par.jcn = JCN_loc;
//     mumps_par.a = A_loc;

//     // Construction of the right hand side
//     // The full rhs is saved on the root.
//     // the rhs will also be used to save the solution
//     //-------------------------------------------------------
//     // Rn rhsMapp;  mapp->reorder(rhs, rhsMapp);
//     //mumps_par.rhs = rhsMapp;
//     mumps_par.rhs = rhs;
//     R ierr =0;
//     // Analysis phase
//     //-------------------------------------------------------
//     //  macro s.t. indices match documentation
// #define INFO(I) info[(I)-1]
//     mumps_par.job = JOB_ANALYSIS;
//     R ta0 = MPI_Wtime();
//     dmumps_c(&mumps_par);
//     R ta1 = MPI_Wtime();
//     ierr = mumps_par.INFO(1);
//     if(ierr != 0) {
//       std::cout << " Error in analysis phase of MUMPS : ierr = " << ierr <<
//       std::endl; std::cout <<mumps_par.INFO(2)  << std::endl; std::cout
//       <<mumps_par.ICNTL(2) << std::endl;
//       // MPIcf::Barrier();
//       assert(0);
//     }

//     // Factorization phase
//     //-------------------------------------------------------
//     mumps_par.job = JOB_FACTORIZATION;
//     R tf0 = MPI_Wtime();
//     dmumps_c(&mumps_par);
//     R tf1 = MPI_Wtime();
//     ierr = mumps_par.INFO(1);

//     if(ierr != 0) {
//       std::cout << " Error in factorization phase of MUMPS : ierr = " << ierr
//       << std::endl;
//       // MPIcf::Barrier();
//       assert(0);
//     }

//     // Solving phase
//     //--------------------------------------------------------
//     R ts0 = MPI_Wtime();
//     mumps_par.job = JOB_SOLVE;
//     dmumps_c(&mumps_par);
//     R ts1 = MPI_Wtime();

//     ierr = mumps_par.INFO(1);
//     if(ierr != 0) {
//       std::cout << " Error in solving phase of MUMPS : ierr = " << ierr <<
//       std::endl;
//       // MPIcf::Barrier();
//       assert(0);
//     }

//     // Terminate instance
//     //--------------------------------------------------------
//     mumps_par.job=JOB_END;
//     dmumps_c(&mumps_par);

//   // Get information about the resolution of the linear system
//   // 9 : Size of the space used to store the factor matrices.
//   // 16 : total size in million of bits of data allocated
//   //      during the factorization
//   //----------------------------------------------------------

//     R szlumn, szlumx, szlu = mumps_par.INFO(9);
//     // MPIcf::AllReduce(szlu, szlumn, MPI_MIN);
//     // MPIcf::AllReduce(szlu, szlumx, MPI_MAX);
//     R szwkmn, szwkmx, szwk = mumps_par.INFO(16);
//     // MPIcf::AllReduce(szwk, szwkmn, MPI_MIN);
//     // MPIcf::AllReduce(szwk, szwkmx, MPI_MAX);
//     // if(MPIcf::IamMaster())
//       // {
//       // 	szlumn = int((szlumn*8.0)/(1024.0*1024.0));
//       // 	szlumx = int((szlumx*8.0)/(1024.0*1024.0));

//       // 	R ratio = ((R)(mumps_par.nz/mumps_par.n)/mumps_par.n)*100.0;
//       // 	std::cout << "
//       -------------------------------------------------------- \n";
//       // 	std::cout << "                MUMPS DIRECT SOLVER " <<
//       std::endl;
//       // 	std::cout << "
//       -------------------------------------------------------- \n";
//       // 	std::cout <<" STATISTICS OF THE GLOBAL MATRIX " << std::endl;
//       // 	std::cout << " Matrix order                         " <<
//       mumps_par.n << std::endl;
//       // 	std::cout << " Number of non-zero entries           " <<
//       mumps_par.nz << std::endl;
//       // 	std::cout << " Fill-in ratio percentage             " << ratio
//       << std::endl;
//       // 	std::cout << "\n STATISTICS OF THE LU FACTORIZATION " <<
//       std::endl;
//       // 	std::cout << " Number of entries in the factors     " <<
//       mumps_par.infog[19] << std::endl;
//       // 	std::cout << "\n Storage of the factors " << std::endl;
//       // 	std::cout << " Memory                       " << szlu <<
//       std::endl;
//       // 	std::cout << "\n Working memory for factorization   " <<
//       std::endl;
//       // 	std::cout << " Memory                       " << szwk <<
//       std::endl;
//       // 	std::cout << std::endl;
//       // 	std::cout << " Time of analysis phase               " << ta1 -
//       ta0 << std::endl;
//       // 	std::cout << " Time of factorization phase          " << tf1 -
//       tf0 << std::endl;
//       // 	std::cout << " Time for solving                     " << ts1 -
//       ts0 << std::endl << std::endl;
//       // }
//       // MPI_Finalize();

//   }

//  void PETSC_Sequential(std::map<std::pair<int,int>,R> & Amap, Rn & rhs) {

//    int argc = 0;
//    char ** argv = NULL;
//    PetscInitialize(&argc, &argv, NULL, NULL);             // initialize PETSC

//    Vec            x, b, D;         /* approx solution, RHS */
//    Mat            A, M;            /* linear system matrix */
//    KSP            solver;       /* linear solver context */
//    PC             pc;           /* preconditioner context */

//    const int size = rhs.size();

//  /*
//    Create vectors.  Note that we form 1 vector from scratch and
//    then duplicate as needed.
//  */

//    VecCreateSeq(PETSC_COMM_WORLD, size, &x);
//    VecDuplicate(x,&b);

//    /*
//      Assemble rhs
//    */
//    for(int i = 0; i < rhs.size(); ++i)
//      VecSetValues(b, 1, &i, &rhs(i), INSERT_VALUES);

//    VecAssemblyBegin(b);
//    VecAssemblyEnd(b);

//    MatCreate(PETSC_COMM_WORLD, &A);
//    MatSetSizes(A,size,size,size,size);
//    MatSetFromOptions(A);
//    MatSetUp(A);

//    /*
//      Assemble matrix from map
//    */
//    const long nbcoef = Amap.size();

//    int* Iv = new int[size+1];
//    int* Jv = new int[nbcoef];
//    R* a = new R[nbcoef];
//    int* nnz = new int[size];
//    R cmm=0;

//    int k=0;
//    for(int i=0;i<=size;++i) Iv[i]=0;        // pour les lignes vide
//    for (typename std::map<std::pair<int,int>,R>::const_iterator
//    q=Amap.begin();
// 	 q != Amap.end(); ++q, ++k)
//      {
// 	int i = q->first.first;
// 	Iv[i+1] = k+1;
// 	Jv[k]   = q->first.second;
// 	a[k]   = q->second;
// 	cmm = std::max(cmm , a[k]);
//      }
//    Amap.clear();
//    for(int i = 1; i <=size; ++i)  Iv[i] = std::max(Iv[i-1],Iv[i]);      //
//    pour les lignes vides for(int i = 0; i < size; ++i)  nnz[i] = Iv[i+1] -
//    Iv[i];             // nonzeros per row

//    assert(k==nbcoef);

//    MatSeqAIJSetPreallocationCSR(A,Iv,Jv,a);

//    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

//    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Create the linear solver and set various options
//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    KSPCreate(PETSC_COMM_WORLD,&solver);

//    /*
//      Set operators. Here the matrix that defines the linear system
//      also serves as the preconditioning matrix.
//    */
//    KSPSetOperators(solver,A,A);

//    //
//    KSPSetTolerances(solver,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

//    // KSPSetType(solver,KSPGMRES);
//    KSPSetType(solver,KSPBCGSL);
//    //KSPSetType(solver,KSPCG);
//    KSPSetInitialGuessNonzero(solver,PETSC_TRUE);
// /*
// 78:    * ILU preconditioner;
// 79:    * this will break down unless you add the Shift line,
// 80:    * or use the -pc_factor_shift_positive_definite option */
//    KSPGetPC(solver,&pc);
//    PCSetType(pc,PCILU);
//    PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE);

//    KSPSetFromOptions(solver);
//    KSPSetUp(solver);

//   /*
// 89:    * Now that the factorisation is done, show the pivots;
// 90:    * note that the last one is negative. This in itself is not an error,
// 91:    * but it will make the iterative method diverge.
// 92:    */
//    PCFactorGetMatrix(pc,&M);
//    VecDuplicate(b,&D);
//    MatGetDiagonal(M,D);

//    /*
//      Solve linear system
//    */
//    // KSPSetFromOptions(solver);

//    KSPSolve(solver,b,x);

//    PetscInt           its;
//    KSPConvergedReason reason;
//    KSPGetConvergedReason(solver,&reason);
//    if (reason==KSP_DIVERGED_INDEFINITE_PC) {
//      PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite
//      preconditioner;\n"); PetscPrintf(PETSC_COMM_WORLD,"Run the executable
//      again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option.\n");
//    } else if (reason<0) {
//      PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should
//      not happen.\n");
//    } else {
//      KSPGetIterationNumber(solver,&its);
//      PetscPrintf(PETSC_COMM_WORLD," Iterations %D\n",its);
//    }
//    // std::cout << " KSP solver()" << std::endl;
//    // KSPSolve(solver,b,x);
//    // PetscReal      normm;
//    // Vec u;
//    // VecDuplicate(x,&u);
//    // MatMult(A,x,u);
//    // VecAXPY(b,-1.0,u);
//    // VecNorm(u,NORM_2,&normm);

//    // // R t3 = MPIcf::Wtime();
//    // // std::cout << " Time PETSC solver \t" << t3 - t2 << std::endl;
//    // PetscInt its;
//    // KSPGetIterationNumber(solver,&its);
//    // // PetscPrintf(PETSC_COMM_WORLD,"iterations %D\n",its);
//    // PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations
//    %D\n",(double)normm,its);

//    // KSPView(solver,PETSC_VIEWER_STDOUT_WORLD);

//    // getchar();
//    for(int i = 0; i < rhs.size(); ++i)
//      VecGetValues(x, 1, &i, &rhs(i));

//   VecDestroy(&D);

//    VecDestroy(&x);
//    VecDestroy(&b);
//    MatDestroy(&A);
//    KSPDestroy(&solver);

//    delete[] Iv;
//    delete[] Jv;
//    delete[] a;
//    delete[] nnz;

//  }

//   void PETSC_Parallel(std::map<std::pair<int,int>,R> & Amap, Rn & rhs) {

//     Vec            x, b;         /* approx solution, RHS */
//     Mat            A;            /* linear system matrix */
//     KSP            solver;          /* linear solver context */
//     MatInfo        matinfo;
//     PC             pc;           /* preconditioner context */

//     /*
//       Get the rows for each processor
//     */
//     const int sizeGlob = rhs.size();
//     // int sizeLoc = rhs.size() / MPIcf::size();
//     // int leftDof = rhs.size() - MPIcf::size()*sizeLoc;
//     // if(MPIcf::IamMaster()) sizeLoc += leftDof;

//     // int Iend;
//     // MPIcf::Scan(sizeLoc, Iend, MPI_SUM);
//     // const int Istart   = Iend - sizeLoc ;
//     // Iend -= 1;

//     int Istart, Iend;
//     int sizeLoc = toolPETSC:: getLocalSize(sizeGlob, Istart, Iend);

//     /*
//       Create vectors.  Note that we form 1 vector from scratch and
//       then duplicate as needed.
//     */
//     VecCreateMPI(PETSC_COMM_WORLD, sizeLoc, sizeGlob, &x);
//     VecSetSizes(x,sizeLoc, sizeGlob);             // local size , global size
//     VecDuplicate(x,&b);

//     Rn rhs_sum(sizeGlob);
//     MPIcf::AllReduce(rhs, rhs_sum, MPI_SUM);
//     for(int i = Istart; i <= Iend; ++i)
//       VecSetValues(b, 1, &i, &rhs_sum(i), INSERT_VALUES);

//     VecAssemblyBegin(b);
//     VecAssemblyEnd(b);

//     /*
//       Create matrix.  When using MatCreate(), the matrix format can
//       be specified at runtime.

//       Performance tuning note:  For problems of substantial size,
//       preallocation of matrix memory is crucial for attaining good
//       performance. See the matrix chapter of the users manual for details.
//     */
//     MatCreate(PETSC_COMM_WORLD, &A);
//     // MatSetSizes(A,sizeLoc,sizeLoc,sizeGlob,sizeGlob);
//     MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,sizeGlob,sizeGlob);

//     MatSetFromOptions(A);
//     MatSetUp(A);

//     // R t0 = MPIcf::Wtime();

//     toolPETSC::transform2PETSCformat(Amap, Istart, Iend);

//     // R t1 = MPIcf::Wtime();
//     // std::cout << " Time transform matrix to PETSC format \t" << t1 - t0 <<
//     std::endl;

//     for(int i=Istart;i<=Iend;++i) Amap[std::make_pair(i,i)] += 0;

//   /*
//     Assemble matrix from map
//   */
//     const int n = sizeLoc;
//     const int m = n;
//     const long nbcoef = Amap.size();

//     int* d_nnz = new int[sizeLoc];
//     int* o_nnz = new int[sizeLoc];

//     KN<int> Iv(n+1); Iv = 0;
//     KN<int> Jv(nbcoef); Jv = 0;
//     R* a = new R[nbcoef];
//     R cmm=0;

//     int k=0;
//     for(int i=0;i<=n;++i) Iv[i]=0;          // pour les lignes vide
//     for(int i=0;i<n;++i) d_nnz[i]=0;        // pour les lignes vide
//     for(int i=0;i<n;++i) o_nnz[i]=0;        // pour les lignes vide

//     for (typename std::map<std::pair<int,int>,R>::const_iterator
//     q=Amap.begin();
// 	 q != Amap.end(); ++q, ++k)
//     {
//       const int i = q->first.first;
//       const int iloc = i - Istart;

//       Iv[iloc+1] = k+1;
//       Jv[k]      = q->first.second;
//       a[k]       = q->second;
//       cmm = std::max(cmm , a[k]);

//       if( Jv[k] < Istart || Jv[k] > Iend )
//       	o_nnz[iloc] += 1;
//       else
//       	d_nnz[iloc] += 1;

//     }
//     Amap.clear();

//     for(int i = 1; i<=n; ++i)  Iv[i] = std::max(Iv[i-1],Iv[i]);  // pour les
//     lignes vides

//     assert(k==nbcoef);
//     MatMPIAIJSetPreallocationCSR(A,Iv,Jv,a);
//     MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

//     // R t2 = MPIcf::Wtime();
//     // std::cout << " Time assembly Matrix in  PETSC solver \t" << t2 - t1 <<
//     std::endl;

//     /*
//       Create linear solver context
//     */
//     KSPCreate(PETSC_COMM_WORLD,&solver);

//     KSPSetOperators(solver,A,A);

//     // KSPGetPC(solver,&pc);

//     // std::cout << " preconditioner " << std::endl;
//     // PCSetType(pc,PCJACOBI);
//     KSPSetTolerances(solver,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

//     KSPSetFromOptions(solver);

//     // std::cout << " here" << std::endl;
//     // getchar();
//     // MPIcf::Barrier();
//     std::cout << " KSP solver()" << std::endl;

//     KSPSolve(solver,b,x);
//     PetscReal      normm;
//     Vec u;
//     VecDuplicate(x,&u);
//     MatMult(A,x,u);
//     VecAXPY(b,-1.0,u);
//     VecNorm(u,NORM_2,&normm);

//     // R t3 = MPIcf::Wtime();
//     // std::cout << " Time PETSC solver \t" << t3 - t2 << std::endl;
//     PetscInt its;
//     KSPGetIterationNumber(solver,&its);
//     // PetscPrintf(PETSC_COMM_WORLD,"iterations %D\n",its);
//     PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations
//     %D\n",(double)normm,its);

//     KSPView(solver,PETSC_VIEWER_STDOUT_WORLD);

//     rhs_sum = 0.0;
//     for(int i = Istart; i <= Iend; ++i)
//       VecGetValues(x, 1, &i, &rhs_sum(i));

//     MPIcf::AllReduce(rhs_sum, rhs , MPI_SUM);

//     VecDestroy(&x);
//     VecDestroy(&b);
//     MatDestroy(&A);
//     KSPDestroy(&solver);

//     delete[] a;
//     delete[] d_nnz;
//     delete[] o_nnz;

//   }

//   void PETSC(std::map<std::pair<int,int>,R> & Amap, Rn & b) {

//     // if(MPIcf::size() < 2)
//       PETSC_Sequential(Amap, b);
//     // else
//     //   PETSC_Parallel(Amap, b);

//   }

// }

// namespace toolPETSC {

//   inline int whatProc(const int i, const int * Iend, const int nproc) {
//     for(int ip = 0; ip < nproc; ip++)
//       if(i <= Iend[ip]) return ip;

//     assert(0 && "domain not found");
//     return 0;
//   }

//   inline void map2AIJ(std::map<std::pair<int,int>,R> & Amap,
// 		      KN<int>& Iv, KN<int>& Jv, KN<double>& a) {

//     Iv.resize(Amap.size());
//     Jv.resize(Amap.size());
//     a.resize(Amap.size());
//     int k = 0;
//     for (typename std::map<std::pair<int,int>,R>::const_iterator
//     q=Amap.begin();
// 	 q != Amap.end(); ++q, ++k)
//       {
// 	Iv(k) = q->first.first;
// 	Jv(k) = q->first.second;
// 	a(k)  = q->second;
//       }
//     Amap.clear();
//   }

//   int getLocalSize(int n, int& Istart, int&Iend) {

//     int sizeLoc = n / MPIcf::size();
//     int leftDof = n - MPIcf::size()*sizeLoc;
//     if(MPIcf::IamMaster()) sizeLoc += leftDof;

//     MPIcf::Scan(sizeLoc, Iend, MPI_SUM);
//     Istart   = Iend - sizeLoc ;
//     Iend -= 1;

//     return sizeLoc;
//   }

//   void transform2PETSCformat(std::map<std::pair<int,int>,R> & Amap,
// 			     int Il, int Ie) {

//     const int sizeLoc = Ie - Il + 1;
//     int Istart[MPIcf::size()], Iend[MPIcf::size()];
//     Istart[MPIcf::my_rank()] = Il;
//     Iend[MPIcf::my_rank()] = Ie;

//     for(int i=0; i<MPIcf::size(); ++i) {
//       MPIcf::Bcast(Istart[i], i, 1);
//       MPIcf::Bcast(Iend[i]  , i, 1);
//     }

//     std::map<std::pair<int,int>,R> m[MPIcf::size()];

//     for (typename std::map<std::pair<int,int>,R>::const_iterator
//     q=Amap.begin();
// 	 q != Amap.end(); ++q)
//     {
//       const int i = q->first.first;
//       const int ip = whatProc(i, Iend, MPIcf::size());
//       if( ip == MPIcf::my_rank()) continue;

//       m[ip][std::make_pair(i,q->first.second)] = q->second;
//       Amap.erase(q);
//     }

//     KN<int> nRecv(MPIcf::size()), nSend(MPIcf::size()); nRecv = 0;
//     for(int i=0;i<MPIcf::size();++i) nSend[i] = m[i].size();

//     for(int i=0;i<MPIcf::size();++i) {
//       MPI_Scatter(nSend, 1, MPI_INT, &nRecv(i),1,MPI_INT, i, MPI_COMM_WORLD);
//     }

//     KN<int> IvRecv[MPIcf::size()], Iv[MPIcf::size()];
//     KN<int> JvRecv[MPIcf::size()], Jv[MPIcf::size()];
//     KN<double> aRecv[MPIcf::size()], a[MPIcf::size()];
//     for(int i = 0; i < MPIcf::size(); ++i) {
//       if(i == MPIcf::my_rank()) continue;
//       IvRecv[i].resize(nRecv[i]);
//       JvRecv[i].resize(nRecv[i]);
//       aRecv[i].resize(nRecv[i]);

//       map2AIJ(m[i], Iv[i], Jv[i], a[i]);

//     }

//     MPI_Status status;
//     MPI_Request request_recv[MPIcf::size()];
//     MPI_Request request_send[MPIcf::size()];
//     MPI_Request request_recvJv[MPIcf::size()];
//     MPI_Request request_recvA[MPIcf::size()];
//     MPI_Request request_sendJv[MPIcf::size()];
//     MPI_Request request_sendA[MPIcf::size()];

//     for(int iproc=0; iproc<MPIcf::size(); ++iproc) {

//       if(iproc == MPIcf::my_rank()) continue;

//       MPI_Isend(Iv[iproc], nSend[iproc], MPI_INT, iproc, MPIcf::my_rank(),
//       MPI_COMM_WORLD, 		&request_send[iproc]);
//       MPI_Irecv(IvRecv[iproc], nRecv[iproc], MPI_INT, iproc, iproc,
//       MPI_COMM_WORLD, 		&request_recv[iproc]);

//       MPI_Isend(Jv[iproc], nSend[iproc], MPI_INT, iproc, MPIcf::my_rank() +
//       10, MPI_COMM_WORLD, 		&request_sendJv[iproc]);
//       MPI_Irecv(JvRecv[iproc], nRecv[iproc], MPI_INT, iproc, iproc + 10,
//       MPI_COMM_WORLD, 		&request_recvJv[iproc]);

//       MPI_Isend(a[iproc], nSend[iproc], MPI_DOUBLE, iproc, MPIcf::my_rank() +
//       20, MPI_COMM_WORLD, 		&request_sendA[iproc]);
//       MPI_Irecv(aRecv[iproc], nRecv[iproc], MPI_DOUBLE, iproc, iproc + 20,
//       MPI_COMM_WORLD, 		&request_recvA[iproc]);

//     }

//     for(int iproc=0; iproc<MPIcf::size(); ++iproc) {
//       if(iproc == MPIcf::my_rank()) continue;
//       MPI_Wait(&request_recv[iproc], &status);
//       MPI_Wait(&request_send[iproc], &status);
//       MPI_Wait(&request_recvJv[iproc], &status);
//       MPI_Wait(&request_sendJv[iproc], &status);
//       MPI_Wait(&request_recvA[iproc], &status);
//       MPI_Wait(&request_sendA[iproc], &status);

//     }

//     for(int iproc=0; iproc<MPIcf::size(); ++iproc) {

//       if(iproc == MPIcf::my_rank()) continue;

//       for(int k=0;k<nRecv[iproc]; ++k) {

// 	const int i = IvRecv[iproc](k);
// 	const int j = JvRecv[iproc](k);
// 	const R val = aRecv[iproc](k);

// 	Amap[std::make_pair(i, j)] += val;

//       }

//     }

//   }

// }
