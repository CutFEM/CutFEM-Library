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
#ifndef SOLVER_CG_HPP_
#define SOLVER_CG_HPP_

#include "solver.hpp"
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>

namespace solver {

// Solve A x = b for SPD A using preconditioned CG.
// On input,  b holds the right-hand side.
// On output, b holds the solution x (same convention as UMFPACK/MUMPS wrappers).
// Preconditioner: Jacobi (diagonal scaling) when use_jacobi=true, otherwise identity.
inline void cg(std::map<std::pair<int,int>, R>& Amap,
               std::span<double> b,
               bool clearMatrix,
               double tol,
               int maxit,
               bool use_jacobi,
               int verbose,
               std::vector<double>* history = nullptr)
{
    const int n = (int)b.size();

    // --- matrix-vector product y = A * x ---
    auto matvec = [&](const std::vector<double>& x, std::vector<double>& y) {
        std::fill(y.begin(), y.end(), 0.0);
        for (const auto& [ij, val] : Amap)
            y[ij.first] += val * x[ij.second];
    };

    // --- build Jacobi preconditioner (diagonal of A) ---
    std::vector<double> diag(n, 1.0);
    if (use_jacobi) {
        for (const auto& [ij, val] : Amap)
            if (ij.first == ij.second && val != 0.0)
                diag[ij.first] = val;
    }
    auto precond = [&](const std::vector<double>& r, std::vector<double>& z) {
        for (int i = 0; i < n; ++i)
            z[i] = r[i] / diag[i];
    };

    // --- helper: dot product ---
    auto dot = [&](const std::vector<double>& u, const std::vector<double>& v) {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += u[i] * v[i];
        return s;
    };

    // --- initialise ---
    // x_0 = 0, r_0 = b
    std::vector<double> x(n, 0.0);
    std::vector<double> r(b.begin(), b.end());
    std::vector<double> z(n), p(n), Ap(n);

    precond(r, z);
    p = z;

    double rz     = dot(r, z);
    double b_norm = std::sqrt(dot(r, r)); // ||b||  (r == b at step 0)
    if (b_norm == 0.0) {
        std::copy(x.begin(), x.end(), b.begin());
        return;
    }

    if (history) history->clear();

    int iter = 0;
    for (; iter < maxit; ++iter) {
        matvec(p, Ap);
        double pAp = dot(p, Ap);

        if (std::abs(pAp) < 1e-300) {
            if (verbose)
                std::cout << "[CG] breakdown: p'Ap ~ 0 at iter " << iter << "\n";
            break;
        }

        double alpha = rz / pAp;
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double res_norm = std::sqrt(dot(r, r));
        double relres   = res_norm / b_norm;
        if (history) history->push_back(relres);
        if (verbose > 1)
            std::cout << "[CG] iter=" << iter << " relres=" << relres << "\n";

        if (relres < tol)
            break;

        precond(r, z);
        double rz_new = dot(r, z);
        double beta   = rz_new / rz;
        rz            = rz_new;

        for (int i = 0; i < n; ++i)
            p[i] = z[i] + beta * p[i];
    }

    if (verbose) {
        double res_norm = std::sqrt(dot(r, r));
        std::cout << "[CG" << (use_jacobi ? "+Jacobi" : "") << "]"
                  << " iters=" << iter
                  << " relres=" << res_norm / b_norm << "\n";
    }

    if (clearMatrix) Amap.clear();

    std::copy(x.begin(), x.end(), b.begin());
}

} // namespace solver

#endif
