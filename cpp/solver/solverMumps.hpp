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
#ifndef SOLVER_MUMPS_HPP_
#define SOLVER_MUMPS_HPP_

#include "solver.hpp"
#include "dmumps_c.h"

enum MatrixFormat { centralized, distributed };
class MUMPS {

    typedef std::map<std::pair<int, int>, R> matmap;

    static const int JOB_INIT_          = -1;
    static const int JOB_ANALYSIS_      = 1;
    static const int JOB_FACTORIZATION_ = 2;
    static const int JOB_SOLVE_         = 3;
    static const int JOB_ALL_           = 6;
    static const int JOB_END_           = -2;
    static const int USE_COMM_WORLD_    = -987654;
    bool cleanMatrix                    = true;
    bool mumps_initialized_             = false;

  public:
    //  int verbose = 0;
    //  std::string reordering = "none";
    MatrixFormat matrixFormat = centralized;

    R timeSettingSolver_, timeAnalysis_, timeFactorization_, timeSolving_;

  private:
    matmap &mat;
    //  Rn &rhs;
    std::span<double> rhs;

    DMUMPS_STRUC_C mumps_par;
    KN<int> IRN_loc, JCN_loc;
    std::vector<double> A_loc, rhsG;

  public:
    // MUMPS(const Solver &, matmap &, Rn &);
    MUMPS(const Solver &, matmap &, std::span<double>);

    MUMPS(matmap &A, std::span<double> b, std::size_t nrhs, bool clean = true);

  private:
    void initializeSetting();
    void setFormatMatrix();
    void setDoF(std::size_t, std::size_t);
    void saveMatrixToCSR();
    void analyzeMatrix();
    void factorizationMatrix();
    void solvingLinearSystem();
    void checkPhase(const char *phase);
    void finalize() noexcept;
    void info();

    int mumps_info(int i) { return mumps_par.info[i - 1]; }
    int mumps_icntl(int i) { return mumps_par.icntl[i - 1]; }

  public:
    ~MUMPS() { finalize(); }
};

#endif
