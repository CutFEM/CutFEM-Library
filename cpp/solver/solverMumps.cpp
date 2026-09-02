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
#include "../common/logger.hpp"
#include "solverMumps.hpp"

#include <sstream>
#include <stdexcept>

#define ICNTL(I) icntl[(I) - 1]
#define INFO(I) info[(I) - 1]
#define INFOG(I) infog[(I) - 1]
MUMPS::MUMPS(const Solver &s, matmap &AA, std::span<double> bb)
    : // verbose(s.verbose),
      //   reordering(s.reordering),
      mat(AA), rhs(bb), cleanMatrix(s.clearMatrix_) {

    //    if(MPIcf::size() > 1) assert(0);
    LOG_INFO << " Using sequential MUMPS  " << logger::endl;

    auto t0 = now();

    initializeSetting();
    std::size_t N = rhs.size();
    setDoF(N, 1);
    saveMatrixToCSR();
    auto t1            = now();
    timeSettingSolver_ = seconds(t0, t1);

    analyzeMatrix();
    factorizationMatrix();
    solvingLinearSystem();

    if (s.verbose_ > 0)
        info();
}

MUMPS::MUMPS(matmap &A, std::span<double> b, std::size_t nrhs, bool clean) : mat(A), rhs(b), cleanMatrix(clean) {
    auto t0 = now();
    initializeSetting();
    std::size_t N = rhs.size() / nrhs;
    setDoF(N, nrhs);
    saveMatrixToCSR();
    auto t1            = now();
    timeSettingSolver_ = seconds(t0, t1);
    analyzeMatrix();
    factorizationMatrix();
    solvingLinearSystem();
}

void MUMPS::setFormatMatrix() { mumps_par.ICNTL(18) = 0; }

void MUMPS::initializeSetting() {
    mumps_par.comm_fortran = USE_COMM_WORLD_;
    // Symmetry of the matrix
    mumps_par.sym          = 0;
    // Type of parallelism (par = 1 : host working - par = 0 : host not working)
    mumps_par.par          = 1;

    // Initialization of one instance of the package
    mumps_par.job = JOB_INIT_;
    dmumps_c(&mumps_par);
    mumps_initialized_ = true;

    //------------------------------------------------------
    mumps_par.ICNTL(1) = -1;
    mumps_par.ICNTL(2) = -1;
    mumps_par.ICNTL(3) = -1;
    mumps_par.ICNTL(4) = -1;

    // Definition of the matrix, distributed assembled format
    //-------------------------------------------------------
    mumps_par.ICNTL(5)  = 0;
    //    setFormatMatrix();
    mumps_par.ICNTL(18) = 0; // 3;

    // Pre-treatment of the matrix (permutation, ordering)
    // 7 : automatic choice
    //------------------------------------------------------
    mumps_par.ICNTL(28) = 1;
    mumps_par.ICNTL(6)  = 7;
    mumps_par.ICNTL(7)  = 3;
    mumps_par.ICNTL(8)  = 77;

    // Increase MAXIS (cf. doc MUMPS) for extra fill-in
    //-------------------------------------------------------
    mumps_par.ICNTL(14) = 70;

    // Format of the right hand side
    //-------------------------------------------------------
    mumps_par.ICNTL(20) = 0;
}

void MUMPS::setDoF(std::size_t n_matrix, std::size_t nrhs) {
    Uint nz_glob = mat.size(), nz_loc = mat.size();

    mumps_par.n      = n_matrix;
    mumps_par.nz     = nz_glob;
    mumps_par.nrhs   = nrhs;
    mumps_par.lrhs   = n_matrix;
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

    // mumps_par.irn_loc = IRN_loc;
    // mumps_par.jcn_loc = JCN_loc;
    // mumps_par.a_loc = A_loc;

    mumps_par.irn = IRN_loc;
    mumps_par.jcn = JCN_loc;
    mumps_par.a   = A_loc.data();

    // Construction of the right hand side
    // The full rhs is saved on the root.
    //-------------------------------------------------------
    mumps_par.rhs = rhs.data();
}

void MUMPS::analyzeMatrix() {
    // timeAnalysis_ = CPUtime();
    auto t0 = now();

    mumps_par.job = JOB_ANALYSIS_;
    dmumps_c(&mumps_par);

    auto t1       = now();
    timeAnalysis_ = seconds(t0, t1);

    LOG_DEBUG << " ordering actually used " << mumps_par.INFOG(7) << logger::endl;
    checkPhase("analysis");
}

void MUMPS::factorizationMatrix() {
    auto t0 = now();

    mumps_par.job = JOB_FACTORIZATION_;
    dmumps_c(&mumps_par);

    auto t1 = now();

    timeFactorization_ = seconds(t0, t1);
    checkPhase("factorization");
}

void MUMPS::solvingLinearSystem() {

    auto t0 = now();

    mumps_par.job = JOB_SOLVE_;
    dmumps_c(&mumps_par);

    auto t1      = now();
    timeSolving_ = seconds(t0, t1);
    checkPhase("solve");

    // Distribution of the DOFs of the solution on each processor
    // (the right-hand side, stored on the host, stores the solution)
    //--------------------------------------------------------
    // if(MPIcf::IamMaster())  {
    //   mapp->inverseMapp(rhsMapp, rhs);
    // }
}

void MUMPS::checkPhase(const char *phase) {
    const int ierr  = mumps_info(1);
    const int info2 = mumps_info(2);
    if (ierr == 0)
        return;

    std::ostringstream message;
    message << "MUMPS " << phase << " failed: INFO(1)=" << ierr << ", INFO(2)=" << info2;
    finalize();
    throw std::runtime_error(message.str());
}

void MUMPS::finalize() noexcept {
    if (!mumps_initialized_)
        return;
    mumps_par.job = JOB_END_;
    dmumps_c(&mumps_par);
    mumps_initialized_ = false;
}

void MUMPS::info() {
    //  Get information about the resolution of the linear system 9 : Size of the space used to store the factor
    //  matrices.16
    //      : total size in million of bits of data allocated during the factorization
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    R szlu = mumps_par.INFO(9);

    //  MPIcf::AllReduce(szlu, szlumn, MPI_MIN);
    //  MPIcf::AllReduce(szlu, szlumx, MPI_MAX);
    R szwk = mumps_par.INFO(16);
    //  MPIcf::AllReduce(szwk, szwkmn, MPI_MIN);
    //  MPIcf::AllReduce(szwk, szwkmx, MPI_MAX);
    if (MPIcf::IamMaster()) {
        szlu = int((szlu * 8.0) / (1024.0 * 1024.0));

        R ratio = ((R)(mumps_par.nz / mumps_par.n) / mumps_par.n) * 100.0;
        LOG_INFO << "\n"
                 << " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- \n "
                 << "                MUMPS DIRECT SOLVER               \n"
                 << " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- \n "
                 << " STATISTICS OF THE GLOBAL MATRIX \n"
                 << " Matrix order " << mumps_par.n << "\n"
                 << " Fill-in ratio percentage " << ratio << "\n"
                 << " Number of non-zero entries           " << mumps_par.nz
                 << "\n STATISTICS OF THE LU FACTORIZATION \n"
                 << " Number of entries in the factors " << mumps_par.infog[19] << "\n"
                 << "\n Storage of the factors \n"
                 << " Memory " << szlu << "\n"
                 << "\n Working memory for factorization   \n"
                 << " Memory                       " << szwk << "\n"
                 << "\n"
                 << " Time to set the solver " << timeSettingSolver_ << "\n"
                 << " Time of analysis phase " << timeAnalysis_ << "\n"
                 << " Time of factorization phase " << timeFactorization_ << "\n"
                 << " Time for solving                     " << timeSolving_ << "\n"
                 << logger::endl;
    }
}

#undef ICNTL
#undef INFO
