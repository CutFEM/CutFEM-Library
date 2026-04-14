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
#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include <map>
#include <list>

#include "../parallel/cfmpi.hpp"
#include "../parallel/cfomp.hpp"
#include "../common/geometry.hpp"
#include "../common/interface_levelSet.hpp"
#include "../common/interface_levelSetP2.hpp"
#include "../common/SparseMatMap.hpp"
#include "../num/DA.hpp"
#include "../FESpace/FESpace.hpp"
#include "../FESpace/expression.hpp"
#include "../FESpace/interpolation.hpp"
#include "../FESpace/restriction.hpp"
#include "../FESpace/integrationFunFEM.hpp"
#include "../FESpace/finiteElement.hpp"
#include "../FESpace/macroElement.hpp"
#include "../FESpace/limiter.hpp"
#include "../num/util.hpp"
#include "itemVF.hpp"
#include "mapping.hpp"
#include "../solver/solver.hpp"

// Base class for problem.
// contain info about the linear system
template <typename Mesh> class ShapeOfProblem {
    typedef std::map<std::pair<int, int>, R> Matrix;
    typedef typename GFESpace<Mesh>::FElement FElement;

  public:
    // The right hand side vector
    // Can be reassign user should not use it!
    std::vector<double> rhs_;

    // // matrix is on a std::map form
    // Matrix mat_;

    // For openmp;
    int thread_count_max_;
    Matrix mat_;
    std::vector<int> index_i0_{0};
    std::vector<int> index_j0_{0};
    std::vector<Matrix> local_contribution_matrix_;

    // pointer on a std::map
    // the user can give is own std::map
    // can be use for newton, matrix that wanna be saved by the user etc
    Matrix *pmat_;

    // std::map used to save CutFEM solution on the background mesh
    // for time dependent problem
    // std::pair(domain, dof_on backSpace) => value
    Matrix mapU0_;

    // local map matrix
    // to reduce the numer of access to the global std::map
    // when the integral is computed on an elements_to_integrate
    // the local contribution is added to the global matrix
    // pointed by pmat
    // Matrix local_contribution_matrix_;
    Rnm loc_mat;
    // number of degree of freedom of the problem
    // This is never modified after initialization
    // => not when adding lagrange multiplier
    Ulint nb_dof_;

    // number of degrees of freedom in a time slab
    // only for time problems
    int nb_dof_time_;

    // index where degree of freedom of consider space starts
    // this make possible to use as many space as we need
    // it is saved in the map
    // given the space it return the where the index start
    std::map<const GFESpace<Mesh> *, int> mapIdx0_;
    // int index_i0_ = 0, index_j0_ = 0;

  public:
    ShapeOfProblem() : nb_dof_(0), nb_dof_time_(1), thread_count_max_(omp_get_max_threads()) {
        set_multithreading_tool();
    };

  private:
    void set_multithreading_tool() {
        // mat_.resize(thread_count_max_);
        local_contribution_matrix_.resize(thread_count_max_);
        index_i0_.resize(thread_count_max_);
        index_j0_.resize(thread_count_max_);
        set_map();
    }

  public:
    long get_size() const { return nb_dof_; }
    long get_nb_dof() const { return nb_dof_; }
    long get_nb_dof_time() const { return nb_dof_time_; }
    int get_num_threads() const { return thread_count_max_; }

    void set_num_thread(int nn) {
        thread_count_max_ = nn;
        set_multithreading_tool();
    }

    void set_map(Matrix &A) {
        assert(thread_count_max_ == 1 && " There are multiple threads. You cannot set only one matrix");
        pmat_ = &A;
    }

    void set_map() { pmat_ = &mat_; }
    void cleanBuildInMatrix() { mat_.clear(); }

    void init(int n) {
        nb_dof_ = n;
        std::fill(rhs_.begin(), rhs_.end(), 0.);
        rhs_.resize(nb_dof_, 0.);
    }
    void init(int n, int nt) {
        nb_dof_      = n;
        nb_dof_time_ = nt;
        std::fill(rhs_.begin(), rhs_.end(), 0.);
        rhs_.resize(nb_dof_);
    }

    void initIndex(const FElement &FKu, const FElement &FKv) {
        int thread_id              = omp_get_thread_num();
        this->index_i0_[thread_id] = mapIdx0_[&FKv.Vh];
        this->index_j0_[thread_id] = mapIdx0_[&FKu.Vh];
    }

    /* This function returns the reference to the element at the specified location in the matrix.
     * Parameters: i, j
     * Return value: reference to the element at the specified location
     */
    R &operator()(int i, int j) { return (*pmat_)[std::make_pair(i + index_i0_[0], j + index_j0_[0])]; }

    // This function returns the i-th element of the right-hand side vector.
    R &operator()(int i) { return rhs_[i + index_i0_[0]]; }

  protected:
    // This function returns a reference to the element (i,j) of the
    // local_contribution_matrix_ member variable. The element (i,j) of
    // the contribution matrix is the element of the matrix that
    // corresponds to the local contribution of the i-th cell on the
    // processor to the j-th cell in the global matrix. The element
    // index_i0_[0] is the index of the first cell on the processor. The
    // element index_j0_[0] is the index of the first cell in the global
    // matrix. This function is used to add a contribution to the
    // element (i,j) of the local_contribution_matrix_.
    double &addToLocalContribution(int i, int j) {
        return local_contribution_matrix_[0][std::make_pair(i + index_i0_[0], j + index_j0_[0])];
    }

    // double &addToLocalContribution_Opt(int i, int j) { return loc_mat(i, j); }

    // This function returns a reference to a double value in the local
    // contribution matrix. The purpose of this function is to allow
    // the caller to modify the value. The key to the local
    // contribution matrix is the pair (i+index_i0_[thread_id],
    // j+index_j0_[thread_id]). The value of the key is the double
    // that is returned by the function. The local contribution matrix
    // is a map of pairs to doubles. The pairs are keys, and the
    // doubles are values. The local contribution matrix is stored
    // in a vector of maps of pairs to doubles. The local contribution
    // matrix is chosen based on the value of thread_id. The index_i0_
    // and index_j0_ are vectors of integers that are used to
    // translate i and j to the correct row and column in the
    // local contribution matrix. This function is used to
    // modify the value of a double in the local contribution
    // matrix.
    double &addToLocalContribution(int i, int j, int thread_id) {
        return local_contribution_matrix_[thread_id]
                                         [std::make_pair(i + index_i0_[thread_id], j + index_j0_[thread_id])];
    }

    void addLocalContribution() {
        int thread_id              = omp_get_thread_num();
        this->index_i0_[thread_id] = 0;
        this->index_j0_[thread_id] = 0;
        auto &A(*pmat_);

#pragma omp critical
        for (const auto &aij : local_contribution_matrix_[thread_id]) {
            A[std::make_pair(aij.first.first, aij.first.second)] += aij.second;
        }
        local_contribution_matrix_[thread_id].clear();
    }

    void addLocalContributionLagrange(int nend) {
        this->index_j0_[0] = 0;
        this->index_i0_[0] = 0;
#pragma omp critical
        for (auto q = local_contribution_matrix_[0].begin(); q != this->local_contribution_matrix_[0].end(); ++q) {
            (*this)(q->first.first, nend) += q->second;
            (*this)(nend, q->first.first) += q->second;
        }
        this->local_contribution_matrix_[0].clear();
    }

  public:
    void applyPreconditioning(std::map<std::pair<int, int>, R> &P) {
        int N = nb_dof_;
        SparseMatrixRC<double> Pl(N, N, P);
        { // P*A -> DF
            SparseMatrixRC<double> A(N, N, mat_);
            multiply(Pl, A, mat_);
        }
        { // (A*P)* -> DF
            SparseMatrixRC<double> A(N, N, mat_);
            multiply(A, Pl, mat_);
        }

        std::vector<R> x(N, 0.);

        multiply(N, N, P, rhs_, x);

        rhs_.resize(N);

        rhs_ = x;
    }
    void recoverSolution(std::map<std::pair<int, int>, R> &P) {
        int N = nb_dof_;
        std::vector<R> x(N, 0.);
        multiply(N, N, P, rhs_, x);
        rhs_ = x;
    }
    void leftPreconditioning(std::map<std::pair<int, int>, R> &P) {

        int N = nb_dof_;
        SparseMatrixRC<double> Pl(N, N, P);
        SparseMatrixRC<double> A(N, N, mat_);
        multiply(Pl, A, mat_);

        std::vector<R> x(N, 0.);

        multiply(N, N, P, rhs_, x);

        rhs_.resize(N);

        rhs_ = x;
    }
    void addMatMul(std::span<double> uuh) {
        assert(uuh.size() == nb_dof_);
        MatriceMap<double> A(nb_dof_, nb_dof_, mat_);
        A.addMatMul(uuh, rhs_);
    }
    // void gather_map() {
    //     for (int i = 1; i < thread_count_max_; ++i) {
    //         auto &A(mat_[i]);
    //         for (const auto &aij : A) {
    //             mat_[0][aij.first] += aij.second;
    //         }
    //         A.clear();
    //     }
    // }
};

static ProblemOption defaultProblemOption;

template <typename Mesh> class QuadratureOfProblem {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename FElement::QFB QFB;
    typedef QuadratureFormular1d QFT;
    typedef typename QFT::QP QPT;

    // CANNOT CHANGE HERE, NEED TO CHANGE IN PROBLEMOPTION
    int order_space_element_quadrature_ = 5;
    int order_space_bord_quadrature_    = 5;

    const QFT *time_quadrature_formular = nullptr;

  public:
    QuadratureOfProblem() {}
    QuadratureOfProblem(const ProblemOption &option) {
        order_space_element_quadrature_ = option.order_space_element_quadrature_;
        order_space_bord_quadrature_    = option.order_space_bord_quadrature_;
    }
    QuadratureOfProblem(const QFT &qt, const ProblemOption &option) : time_quadrature_formular(&qt) {
        order_space_element_quadrature_ = option.order_space_element_quadrature_;
        order_space_bord_quadrature_    = option.order_space_bord_quadrature_;
    }

    const QF &get_quadrature_formular_K() const;
    const QFB &get_quadrature_formular_dK() const;

    const QF &get_quadrature_formular_cutK() const;
    const QFB &get_quadrature_formular_cutFace() const;

    int get_nb_quad_point_time() const { return (time_quadrature_formular) ? time_quadrature_formular->n : 1; }
    QPT get_quadrature_time(int itq) const {
        assert(itq < get_nb_quad_point_time());
        return (time_quadrature_formular) ? time_quadrature_formular->at(itq) : QPT(1., 0.);
    }

    int get_quad_order_space() const {return order_space_element_quadrature_; }
    int get_quad_order_border() const {return order_space_bord_quadrature_; }
};

#endif
