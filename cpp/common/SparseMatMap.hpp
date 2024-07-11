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
/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef COMMON_SPARSE_MATRIX_HPP
#define COMMON_SPARSE_MATRIX_HPP

#include <map>

#include "RNM.hpp"
#include "../num/util.hpp"
#include <map>
#include <set>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cassert>

typedef std::map<std::pair<int, int>, R> Matrix;

void buil_CSR_array(int n, const std::map<std::pair<int, int>, R> &M, int32_t *p, int32_t *j, double *a);
void multiply(int N, const std::map<std::pair<int, int>, double> &A, const std::map<std::pair<int, int>, double> &B,
              std::map<std::pair<int, int>, double> &C);
void eraseAndSetRow(int N, std::map<std::pair<int, int>, double> &A, std::span<R> b, std::map<int, double> &dof2rm);
void eraseAndSetRow(int N, std::map<std::pair<int, int>, double> &A, std::span<R> b, int, int, double);
void eraseRow(int N, std::map<std::pair<int, int>, double> &A, std::span<R> b, std::set<int> &dof2rm);

template <typename Vec>
void multiply(int N, int M, const std::map<std::pair<int, int>, double> &A, const Vec &rhs, std::span<R> b) {
    std::fill(b.begin(), b.end(), 0);
    assert(b.size() == N && rhs.size() == M);
    auto itA = A.begin();
    while (itA != A.end()) {
        b[itA->first.first] += itA->second * rhs[itA->first.second];
        itA++;
    }
}

template <typename Vec> void multiply(int N, int M, const std::vector<Matrix> &A, const Vec &rhs, std::span<R> b) {
    std::fill(b.begin(), b.end(), 0);
    assert(b.size() == N && rhs.size() == M);

    auto itA = A.begin()->begin();
    while (itA != A.begin()->end()) {
        b[itA->first.first] += itA->second * rhs[itA->first.second];
        itA++;
    }
}

template <class R> struct VirtualMatrice {
  public:
    int N, M;
    VirtualMatrice(int nn, int mm) : N(nn), M(mm) {}
    VirtualMatrice(int nn) : N(nn), M(nn) {}
    virtual void addMatMul(std::span<R> x, std::span<R> y) const = 0;
    virtual ~VirtualMatrice() {}
};

template <class R> class SparseMatrixRC : public VirtualMatrice<R> {
  public:
    int n, m, nbcoef;
    // int32_t *p;
    // int32_t *j;
    int64_t *p;
    int64_t *j;
    double *a;
    SparseMatrixRC(int nn, int mm, const std::map<std::pair<int, int>, R> &M);

    R *pcoef(int ii, int jj);
    const R &operator()(int i, int j) const {
        R *p = pcoef(i, j);
        assert(p);
        return *p;
    }
    R &operator()(int i, int j) {
        R *p = pcoef(i, j);
        assert(p);
        return *p;
    }
    R &operator[](std::pair<int, int> ij) {
        R *p = pcoef(ij.first, ij.second);
        assert(p);
        return *p;
    }
    const R &operator[](std::pair<int, int> ij) const {
        R *p = pcoef(ij.first, ij.second);
        assert(p);
        return *p;
    }

    std::vector<double> getRow(int i) const {
        int n = p[i + 1] - p[i];
        std::vector<double> v(n);
        v.resize(n);
        for (int k = p[i], j = 0; k < p[i + 1]; ++k, ++j) {
            v[j] = a[k];
        }
        return v;
    }
    void addMatMul(std::span<R> x, std::span<R> Ax) const {
        for (int i = 0; i < n; ++i)
            for (int k = p[i]; k < p[i + 1]; ++k)
                Ax[i] += a[k] * x[j[k]];
    }
    ~SparseMatrixRC() {
        delete[] p;
        delete[] j;
        delete[] a;
    }

  private:
    SparseMatrixRC(const SparseMatrixRC &);
    void operator=(const SparseMatrixRC &);
};

void multiply(const SparseMatrixRC<double> &A, const SparseMatrixRC<double> &B,
              std::map<std::pair<int, int>, double> &C);

template <class R> class MatriceMap : VirtualMatrice<R> {
  public:
    typedef std::map<std::pair<int, int>, R> Map;

    const Map &m;
    MatriceMap(int in, int im, const Map &mm) : VirtualMatrice<R>(in, im), m(mm) {}
    void addMatMul(std::span<R> x, std::span<R> Ax) const;
};

template <class R> std::vector<R> operator*(const std::map<std::pair<int, int>, R> &A, std::span<R> x) {
    std::vector<R> Ax(x.size());
    auto last_index = A.end();
    last_index--;
    int row_end = last_index->first.first;
    assert(row_end < x.size());
    for (const auto &aij : A) {
        int i = aij.first.first;
        int j = aij.first.second;
        Ax[i] += aij.second * x[j];
    }
    return Ax;
}

template <class R> void MatriceMap<R>::addMatMul(std::span<R> x, std::span<R> Ax) const {
    for (typename Map::const_iterator k = this->m.begin(); k != this->m.end(); ++k) {
        int i = k->first.first;
        int j = k->first.second;
        R aij = k->second;
        Ax[i] += aij * x[j];
    }
}

template <class R>
SparseMatrixRC<R>::SparseMatrixRC(int nn, int mm, const std::map<std::pair<int, int>, R> &M)
    : VirtualMatrice<R>(nn, mm), n(nn), m(mm), nbcoef(M.size()), p(new long long[nn + 1]), j(new long long[nbcoef]),
      a(new R[nbcoef]) {
    R cmm = 0;

    int k = 0;
    for (int i = 0; i <= n; ++i)
        p[i] = 0; // pour les lignes vide
    for (typename std::map<std::pair<int, int>, R>::const_iterator q = M.begin(); q != M.end(); ++q, ++k) {
        int i      = q->first.first;
        p[i + 1]   = k + 1;
        this->j[k] = q->first.second;
        this->a[k] = q->second;
        cmm        = std::max(cmm, this->a[k]);
    }
    for (int i = 1; i <= n; ++i)
        p[i] = std::max(p[i - 1], p[i]); // pour les lignes vides
    assert(k == this->nbcoef);
}

template <class R>
SparseMatrixRC<R>::SparseMatrixRC(const SparseMatrixRC<R> &A)
    : VirtualMatrice<R>(A.n, A.m), n(A.n), m(A.m), nbcoef(A.nbcoef), p(new long long[n + 1]), j(new long long[nbcoef]),
      a(new R[nbcoef]) {
    R cmm = 0;

    for (int i = 0; i <= n; ++i)
        p[i] = 0; // pour les lignes vide
    for (int k = 0; k < nbcoef; ++k) {
        int i      = A.i[k];
        p[i + 1]   = k + 1;
        this->j[k] = A.j[k];
        this->a[k] = A.a[k];
        cmm        = std::max(cmm, this->a[k]);
    }
    for (int i = 1; i <= n; ++i)
        p[i] = std::max(p[i - 1], p[i]); // pour les lignes vides
    std::cout << " nb coef = " << nbcoef << " c max = " << cmm << std::endl;
}

template <class R> R *SparseMatrixRC<R>::pcoef(int ii, int jj) {
    if (ii < 0 && jj >= n)
        return 0;
    int k0 = p[ii], k1 = p[ii + 1] - 1;
    while (k0 <= k1) {
        int km = (k0 + k1) / 2;
        int aa = 0;
        if (j[km] > jj)
            aa = -1;
        else if (j[km] < jj)
            aa = 1;
        if (aa < 0)
            k1 = km - 1;
        else if (aa > 0)
            k0 = km + 1;
        else {
            return a + km;
        }
    }
    return 0;
}

#endif