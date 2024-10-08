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
#include "solverMumps.hpp"
#include "../num/matlab.hpp"
#include <thread>
#include "../num/print_container.hpp"
#include "../common/logger.hpp"

#define ICNTL(I) icntl[(I) - 1]
#define INFO(I) info[(I) - 1]

MUMPS::MUMPS(const Solver &s, matmap &AA, std::span<double> bb)
    :

      // verbose(s.verbose),
      // reordering(s.reordering),
      mat(AA), rhs(bb), cleanMatrix(s.clearMatrix_) {

    //    if(MPIcf::size() > 1) assert(0);
    LOG_INFO << "MUMPS solver is used" << logger::endl;
    initializeSetting();
    setDoF();
    saveMatrixToCSR();
    analyzeMatrix();
    factorizationMatrix();
    solvingLinearSystem();

    // if(verbose > 1)
    // info();
}

void MUMPS::setFormatMatrix() {
    // if(matrixFormat == centralized) {
    if (MPIcf::size() == 1) {
        mumps_par.ICNTL(18) = 0;
    } else {
        mumps_par.ICNTL(18) = 3; // distributed
    }
}

void MUMPS::initializeSetting() {
    mumps_par.comm_fortran = USE_COMM_WORLD_;
    // Symmetry of the matrix
    mumps_par.sym          = 0;
    // Type of parallelism (par = 1 : host working - par = 0 : host not working)
    mumps_par.par          = 1;

    // Initialization of one instance of the package
    mumps_par.job = JOB_INIT_;
    dmumps_c(&mumps_par);

    //------------------------------------------------------
    mumps_par.ICNTL(1) = -1;
    mumps_par.ICNTL(2) = -1;
    mumps_par.ICNTL(3) = -1;
    mumps_par.ICNTL(4) = -1;

    // Definition of the matrix, distributed assembled format
    //-------------------------------------------------------
    mumps_par.ICNTL(5)  = 0;
    //    setFormatMatrix();
    mumps_par.ICNTL(18) = 3; // 3;

    // Pre-treatment of the matrix (permutation, ordering)
    // 7 : automatic choice
    //------------------------------------------------------
    mumps_par.ICNTL(28) = 1;
    mumps_par.ICNTL(6)  = 7;
    mumps_par.ICNTL(7)  = 3;
    mumps_par.ICNTL(8)  = 77;

    // Increase MAXIS (cf. doc MUMPS) for extra fill-in
    //-------------------------------------------------------
    mumps_par.ICNTL(14) = 100;

    // Format of the right hand side
    //-------------------------------------------------------
    mumps_par.ICNTL(20) = 0;
}

void MUMPS::setDoF() {

    Uint nz_glob = 0, nz_loc = mat.size();
    MPIcf::AllReduce(nz_loc, nz_glob, MPI_SUM);

    if (MPIcf::IamMaster()) {
        mumps_par.n  = rhs.size();
        mumps_par.nz = nz_glob;
    }
    mumps_par.nz_loc = nz_loc;
}

void MUMPS::saveMatrixToCSR() {

    // Construction of the local matrices for MUMPS
    //-------------------------------------------------------
    IRN_loc.init(mumps_par.nz_loc);
    JCN_loc.init(mumps_par.nz_loc);
    A_loc.resize(mumps_par.nz_loc);

    int k = 0;
    for (std::map<std::pair<int, int>, R>::const_iterator q = mat.begin(); q != mat.end(); ++q, ++k) {
        // !!!!!!     FORTRAN NUMBERING    !!!!!!!!!!
        // IRN_loc(k) = mapp->iperm(q->first.first)  + 1;
        // JCN_loc(k) = mapp->iperm(q->first.second) + 1;
        IRN_loc(k) = q->first.first + 1;
        JCN_loc(k) = q->first.second + 1;
        A_loc[k]   = q->second;
    }
    if (cleanMatrix)
        mat.clear();

    mumps_par.irn_loc = IRN_loc;
    mumps_par.jcn_loc = JCN_loc;
    mumps_par.a_loc   = A_loc.data();

    // mumps_par.irn = IRN_loc;
    // mumps_par.jcn = JCN_loc;
    // mumps_par.a = A_loc;

    // Construction of the right hand side
    // The full rhs is saved on the root.
    //-------------------------------------------------------
    // rhsG.init(rhs.size()*(MPIcf::IamMaster()));             // alloc mem for
    // global rhs
    rhsG.resize(rhs.size()); // alloc mem for global rhs
    std::copy(rhs.begin(), rhs.end(), rhsG.begin());
    // rhsG = rhs;
    MPIcf::Reduce<double>(rhsG, rhs, MPI_SUM, MPIcf::Master());
    // MPI_Reduce(rhsG.data(), rhs.data(), rhs.size(), MPI_TYPE<double>::TYPE(), MPI_SUM, MPIcf::Master(),
    // MPI_COMM_WORLD);

    //   //    Rn rhsMapp;
    if (MPIcf::IamMaster()) {
        //      mapp->reorder(rhs_sum, rhsMapp);
        mumps_par.rhs = rhs.data(); // rhsMapp;
    }

    // mumps_par.rhs = rhs;
}

void MUMPS::analyzeMatrix() {
    timeAnalysis_ = MPIcf::Wtime();

    mumps_par.job = JOB_ANALYSIS_;
    dmumps_c(&mumps_par);

    R ierr        = mumps_info(1);
    timeAnalysis_ = MPIcf::Wtime() - timeAnalysis_;

    if (ierr != 0) {
        std::cout << " Error in analysis phase of MUMPS : ierr = " << ierr << std::endl;
        std::cout << mumps_par.INFO(2) << std::endl;
        std::cout << mumps_par.ICNTL(2) << std::endl;
        MPIcf::Barrier();
    }
}

void MUMPS::factorizationMatrix() {
    timeFactorization_ = MPIcf::Wtime();

    mumps_par.job = JOB_FACTORIZATION_;
    dmumps_c(&mumps_par);

    R ierr = mumps_info(1);

    timeFactorization_ = MPIcf::Wtime() - timeFactorization_;

    if (ierr != 0) {
        std::cout << " Error in factorization phase of MUMPS : ierr = " << ierr << std::endl;
        std::cout << " info(2) \t" << mumps_info(2) << std::endl;
        MPIcf::Barrier();
    }
}

void MUMPS::solvingLinearSystem() {

    timeSolving_ = MPIcf::Wtime();

    mumps_par.job = JOB_SOLVE_;
    dmumps_c(&mumps_par);

    R ierr = mumps_info(1);

    timeSolving_ = MPIcf::Wtime() - timeSolving_;

    if (ierr != 0) {
        std::cout << " Error in solving phase of MUMPS : ierr = " << ierr << std::endl;
        MPIcf::Barrier();
    }

    // Distribution of the DOFs of the solution on each processor
    // (the right-hand side, stored on the host, stores the solution)
    //--------------------------------------------------------
    // if(MPIcf::IamMaster())  {
    //   mapp->inverseMapp(rhsMapp, rhs);
    // }
    // if(!MPIcf::IamMaster()) rhs = 0.;
    // std::cout << rhs.size() << std::endl;
    // std::cout << rhs << std::endl;
    // getchar();
    // MPIcf::Barrier();

    MPIcf::Bcast(rhs, MPIcf::Master());
    // matlab::Export(rhs,
    // "sol"+to_string(MPIcf::my_rank())+"_"+to_string(MPIcf::size())+".dat");
    // MPIcf::Barrier();
}

void MUMPS::info() {
    // Get information about the resolution of the linear system
    // 9 : Size of the space used to store the factor matrices.
    // 16 : total size in million of bits of data allocated
    //      during the factorization
    //----------------------------------------------------------
    R szlumn, szlumx, szlu = mumps_par.INFO(9);
    MPIcf::AllReduce(szlu, szlumn, MPI_MIN);
    MPIcf::AllReduce(szlu, szlumx, MPI_MAX);
    R szwkmn, szwkmx, szwk = mumps_par.INFO(16);
    MPIcf::AllReduce(szwk, szwkmn, MPI_MIN);
    MPIcf::AllReduce(szwk, szwkmx, MPI_MAX);
    if (MPIcf::IamMaster()) {
        szlumn = int((szlumn * 8.0) / (1024.0 * 1024.0));
        szlumx = int((szlumx * 8.0) / (1024.0 * 1024.0));

        R ratio = ((R)(mumps_par.nz / mumps_par.n) / mumps_par.n) * 100.0;
        std::cout << " -------------------------------------------------------- \n";
        std::cout << "                MUMPS DIRECT SOLVER               " << std::endl;
        std::cout << " -------------------------------------------------------- \n";
        std::cout << " STATISTICS OF THE GLOBAL MATRIX " << std::endl;
        std::cout << " Matrix order                         " << mumps_par.n << std::endl;
        std::cout << " Number of non-zero entries           " << mumps_par.nz << std::endl;
        std::cout << " Fill-in ratio percentage             " << ratio << std::endl;
        std::cout << "\n STATISTICS OF THE LU FACTORIZATION " << std::endl;
        std::cout << " Number of entries in the factors     " << mumps_par.infog[19] << std::endl;
        std::cout << "\n Storage of the factors " << std::endl;
        std::cout << " Memory                       " << szlu << std::endl;
        std::cout << "\n Working memory for factorization   " << std::endl;
        std::cout << " Memory                       " << szwk << std::endl;
        std::cout << std::endl;
        std::cout << " Time of analysis phase               " << timeAnalysis_ << std::endl;
        std::cout << " Time of factorization phase          " << timeFactorization_ << std::endl;
        std::cout << " Time for solving                     " << timeSolving_ << std::endl << std::endl;
    }
}

#undef ICNTL
#undef INFO
