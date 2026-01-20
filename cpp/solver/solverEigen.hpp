#ifndef SOLVER_EIGEN_HPP_
#define SOLVER_EIGEN_HPP_

#include "solver.hpp"

#ifdef CUTFEM_USE_EIGEN
  #include <Eigen/Sparse>
  #include <Eigen/IterativeLinearSolvers>
  #include <vector>
  #include <iostream>
#endif

namespace solver {

#ifdef CUTFEM_USE_EIGEN

namespace detail {
using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Trip  = Eigen::Triplet<double>;

inline SpMat map_to_sparse(int n, const std::map<std::pair<int,int>, R>& Amap) {
    std::vector<Trip> trip;
    trip.reserve(Amap.size());

    for (const auto& it : Amap) {
        const int i = it.first.first;
        const int j = it.first.second;
        const double a = static_cast<double>(it.second);
        if (a != 0.0) trip.emplace_back(i, j, a);
    }

    SpMat A(n, n);
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    return A;
}

inline Eigen::VectorXd span_to_vec(std::span<double> b) {
    Eigen::VectorXd v((int)b.size());
    for (int i = 0; i < v.size(); ++i) v[i] = b[(size_t)i];
    return v;
}

inline void vec_to_span(const Eigen::VectorXd& x, std::span<double> b) {
    for (int i = 0; i < x.size(); ++i) b[(size_t)i] = x[i];
}
} // namespace detail

// Overwrites b with x, like UMFPACK/MUMPS wrappers do.
// Assumes SPD.
inline void eigen_cg(std::map<std::pair<int,int>, R>& Amap,
                     std::span<double> b,
                     bool clearMatrix,
                     double tol,
                     int maxit,
                     bool use_incomplete_cholesky,
                     int verbose)
{
    const int n = (int)b.size();

    auto A = detail::map_to_sparse(n, Amap);
    if (clearMatrix) Amap.clear();

    const auto rhs = detail::span_to_vec(b);

    if (use_incomplete_cholesky) {
        Eigen::ConjugateGradient<detail::SpMat,
                                 Eigen::Lower|Eigen::Upper,
                                 Eigen::IncompleteCholesky<double>> cg;
        cg.setTolerance(tol);
        cg.setMaxIterations(maxit);
        cg.compute(A);

        const Eigen::VectorXd x = cg.solve(rhs);
        detail::vec_to_span(x, b);

        if (verbose) {
            std::cout << "[Eigen CG+IC] iters=" << cg.iterations()
                      << " relres=" << cg.error() << "\n";
        }
    } else {
        Eigen::ConjugateGradient<detail::SpMat,
                                 Eigen::Lower|Eigen::Upper,
                                 Eigen::DiagonalPreconditioner<double>> cg;
        cg.setTolerance(tol);
        cg.setMaxIterations(maxit);
        cg.compute(A);

        const Eigen::VectorXd x = cg.solve(rhs);
        detail::vec_to_span(x, b);

        if (verbose) {
            std::cout << "[Eigen CG+Diag] iters=" << cg.iterations()
                      << " relres=" << cg.error() << "\n";
        }
    }
}

#endif // CUTFEM_USE_EIGEN

} // namespace solver

#endif
