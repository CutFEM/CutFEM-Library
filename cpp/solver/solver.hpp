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
// #include "../cutFEMConfig.h"
#include "../num/util.hpp"
#include "../parallel/cfmpi.hpp"

struct ProblemOption {
    int order_space_element_quadrature_ = 5;
    int algoim_bernstein_deg_           = 2;
    int algoim_surface_quad_deg_        = 5;
    int algoim_vol_quad_deg_            = 5;
    // Use Saye's multivariate-polynomial quadrature on MeshQuad2 cells.  The
    // default keeps the general smooth-level-set quadGen path for existing
    // drivers; polynomial moving-interface drivers can opt in explicitly.
    bool algoim_use_multipoly_           = false;
    // 1D quadrature strategy for the multipoly (triangle/tet) algoim path:
    // 0 = AlwaysGL, 1 = AlwaysTS, 2 = AutoMixed (algoim::QuadStrategy values).
    // AutoMixed falls back to tanh-sinh on base integrals with detected
    // vertical tangents, whose accuracy at moderate q is ~1e-5..1e-8 -- far
    // below the Gauss-Legendre branch.  Keep AutoMixed for robustness on wild
    // cuts; use AlwaysGL when the quadrature consistency error must be small
    // (e.g. pressure-robust Stokes at small nu).
    int algoim_quad_strategy_           = 2;
    // Post-process the Mesh2 (triangle multipoly) cut rules so that the
    // discrete divergence theorem  int_{K cap Omega} div F = int_{K cap Gamma}
    // F.n + int_{faces cap Omega} F.n  holds EXACTLY (to solver precision) for
    // all polynomial fields F up to degree algoim_ibp_degree_, elementwise.
    // Two-stage constrained least-squares correction: (1) the surface rule's
    // vector weights w*n are corrected against exact 1D face integrals using
    // divergence-free polynomial test fields; (2) the volume weights are
    // moment-fitted to the divergence-theorem moments defined by the corrected
    // boundary rule.  This removes the quadrature inconsistency that is
    // amplified by |p|/nu in pressure-robust Stokes (see
    // workfiles/src/stokes/algoim_quadrature_probe.cpp).
    bool algoim_ibp_consistent_         = false;
    int algoim_ibp_degree_              = 6;
    int order_space_bord_quadrature_    = 5;
    int order_time_quadrature_          = 3;
    std::string solver_name_            = "mumps";
    bool clear_matrix_                  = true;
    int verbose_                        = 0;

    // --- Iterative solver options (used e.g. by Eigen CG) ---
    int    it_maxit_    = 2000;
    double it_tol_      = 1e-10;
    bool   it_use_ic_   = true;   // incomplete Cholesky (SPD)
};

namespace solver {

// void umfpack(std::map<std::pair<int, int>, R> &, Rn &, bool);
void umfpack(std::map<std::pair<int, int>, R> &, std::span<double>, bool);

void LAPACK(Rnm &a, Rn &b);

void mumps(std::map<std::pair<int, int>, R> &A, std::span<double> b, std::size_t nrhs);

} // namespace solver
class Solver {
    double get_Time() const { return MPIcf::Wtime(); }

  public:
    int verbose_             = 0;
    bool clearMatrix_        = true;
    std::string solver_name_ = "default";

    // iterative options
    int    it_maxit_   = 2000;
    double it_tol_     = 1e-10;
    bool   it_use_ic_  = true;

    // populated after each CG solve; empty for direct solvers
    std::vector<double> residual_history_;

    Solver(const ProblemOption &option) {
        clearMatrix_ = option.clear_matrix_;
        solver_name_ = option.solver_name_;
        verbose_     = option.verbose_;

        it_maxit_  = option.it_maxit_;
        it_tol_    = option.it_tol_;
        it_use_ic_ = option.it_use_ic_;
    }

    void solve(std::map<std::pair<int, int>, R> &A, std::span<double> b);
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
