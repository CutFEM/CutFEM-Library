/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
*/

#ifndef SOLVER_CG_HPP_
#define SOLVER_CG_HPP_

#include "solver.hpp"
#include <cmath>
#include <vector>
#include <map>
#include <span>
#include <iostream>
#include <algorithm>

namespace solver {

/**
 * @brief Compressed Sparse Row (CSR) storage for efficient Mat-Vec products.
 */
struct CSRMatrix {
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_ptr;
    int n = 0;

    void clear() {
        values.clear();
        col_indices.clear();
        row_ptr.clear();
        n = 0;
    }
};

/**
 * @brief Utility to convert a coordinate map to CSR format.
 * While the conversion has overhead, the resulting CSR structure 
 * speeds up CG iterations by orders of magnitude.
 */
inline CSRMatrix mapToCSR(std::map<std::pair<int, int>, double>& Amap, int n) {
    CSRMatrix csr;
    csr.n = n;
    csr.row_ptr.reserve(n + 1);
    csr.values.reserve(Amap.size());
    csr.col_indices.reserve(Amap.size());

    int current_row = 0;
    csr.row_ptr.push_back(0);

    // We iterate through rows 0 to n-1
    // std::map with std::pair<int,int> is sorted by row then column
    auto it = Amap.begin();
    for (int i = 0; i < n; ++i) {
        while (it != Amap.end() && it->first.first == i) {
            csr.values.push_back(it->second);
            csr.col_indices.push_back(it->first.second);
            ++it;
        }
        csr.row_ptr.push_back(static_cast<int>(csr.values.size()));
    }
    return csr;
}

/**
 * @brief Solve Ax = b for SPD A using Preconditioned Conjugate Gradient.
 * * @param Amap         Input matrix in map format (will be converted to CSR internally).
 * @param b            On input: RHS. On output: Solution x.
 * @param clearMatrix  If true, clears the input Amap to save memory.
 * @param tol          Relative tolerance.
 * @param maxit        Maximum iterations.
 * @param use_jacobi   Enable diagonal scaling preconditioner.
 * @param verbose      Verbosity level (0: silent, 1: summary, 2: per-iter).
 * @param history      Optional pointer to store convergence history.
 */
inline void cg(std::map<std::pair<int, int>, double>& Amap,
               std::span<double> b,
               bool clearMatrix,
               double tol,
               int maxit,
               bool use_jacobi,
               int verbose,
               std::vector<double>* history = nullptr)
{
    const int n = static_cast<int>(b.size());
    
    // 1. Convert to efficient storage
    CSRMatrix A = mapToCSR(Amap, n);
    if (clearMatrix) Amap.clear();

    // 2. Build Jacobi Preconditioner (Inverse Diagonal)
    std::vector<double> inv_M(n, 1.0);
    if (use_jacobi) {
        for (int i = 0; i < n; ++i) {
            for (int jj = A.row_ptr[i]; jj < A.row_ptr[i + 1]; ++jj) {
                if (A.col_indices[jj] == i) {
                    double diag_val = A.values[jj];
                    inv_M[i] = (std::abs(diag_val) > 1e-18) ? 1.0 / diag_val : 1.0;
                    break;
                }
            }
        }
    }

    // 3. Helper Lambdas
    auto matvec = [&](const std::vector<double>& x_in, std::vector<double>& y_out) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int jj = A.row_ptr[i]; jj < A.row_ptr[i + 1]; ++jj) {
                sum += A.values[jj] * x_in[A.col_indices[jj]];
            }
            y_out[i] = sum;
        }
    };

    auto dot = [](const std::vector<double>& u, const std::vector<double>& v) {
        double s = 0.0;
        for (size_t i = 0; i < u.size(); ++i) s += u[i] * v[i];
        return s;
    };

    // 4. Initialisation
    std::vector<double> x(n, 0.0);
    std::vector<double> r(b.begin(), b.end()); // r = b - Ax_0 (with x_0=0)
    std::vector<double> z(n), p(n), Ap(n);

    // Apply preconditioner z = M^-1 r
    for (int i = 0; i < n; ++i) z[i] = r[i] * inv_M[i];
    p = z;

    double rz = dot(r, z);
    double b_norm_sq = dot(r, r);
    double b_norm = std::sqrt(b_norm_sq);

    if (b_norm < 1e-18) {
        std::fill(b.begin(), b.end(), 0.0);
        return;
    }

    if (history) history->clear();

    // 5. Main CG Loop
    int iter = 0;
    for (; iter < maxit; ++iter) {
        matvec(p, Ap);
        double pAp = dot(p, Ap);

        if (std::abs(pAp) < 1e-20) {
            if (verbose) std::cout << "[CG] Breakdown: pAp too small.\n";
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

        if (verbose > 1) {
            std::cout << "[CG] iter=" << iter << " relres=" << relres << "\n";
        }

        if (relres < tol) break;

        // Precondition: z = M^-1 r
        for (int i = 0; i < n; ++i) z[i] = r[i] * inv_M[i];

        double rz_new = dot(r, z);
        double beta = rz_new / rz;
        rz = rz_new;

        for (int i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
    }

    if (verbose) {
        std::cout << "[CG" << (use_jacobi ? "+Jacobi" : "") << "] "
                  << "iters=" << iter << " final_relres=" << std::sqrt(dot(r, r)) / b_norm << "\n";
    }

    // Copy result back to b
    std::copy(x.begin(), x.end(), b.begin());
}

} // namespace solver

#endif