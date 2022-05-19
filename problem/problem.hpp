#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include <map>
#include <list>
#include "../common/geometry.hpp"
// #include "../common/marker.hpp"
#include "../common/interface_levelSet.hpp"
#include "../common/GTime.hpp"
#include "../common/SparseMatMap.hpp"
#include "../common/DA.hpp"
#include "../FESpace/FESpace.hpp"
#include "../FESpace/expression.hpp"
#include "../FESpace/interpolation.hpp"
#include "../FESpace/restriction.hpp"
#include "../FESpace/integrationFunFEM.hpp"
#include "../FESpace/finiteElement.hpp"
#include "../FESpace/macroElement.hpp"
#include "../FESpace/limiter.hpp"
#include "../util/util.hpp"
#include "../util/cputime.h"
#include "itemVF.hpp"
#include "../solver/solver.hpp"



struct CBorder {
  CBorder() {}
};
const  CBorder boundary;
struct CFacet {
  CFacet() {}
};
const CFacet innerFacet;
struct CRidge {
  CRidge() {}
};
const CRidge innerRidge;


// Base class for problem.
// contain info about the linear system
template<typename Mesh>
class ShapeOfProblem {
  typedef std::map<std::pair<int,int>,R> Matrix;

  typedef typename GFESpace<Mesh>::FElement FElement;

public:
  // The right hand side vector
  // Can be reassign user should not use it!
  KN<double> rhs_;

  // matrix is on a std::map form
  Matrix mat_;


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
  Matrix local_contribution_matrix_;

  // number of degree of freedom of the problem
  // This is never modify after initialization
  // => not when adding lagrange multiplier
  Ulint nb_dof_;

  // number of degrees of freedom in a time slab
  // only for time problems
  int nb_dof_time_;

  // index where degree of freedom of consider space starts
  // this make possible to use as many space as we need
  // it is saved in the map
  // given the space it return the where the index start
  std::map<const GFESpace<Mesh>*, int> mapIdx0_;
  int index_i0_ = 0, index_j0_ = 0;


public :
  ShapeOfProblem() : nb_dof_(0), nb_dof_time_(1) { pmat_ = &mat_;};
  ShapeOfProblem(int n) : nb_dof_(n), rhs_(n), nb_dof_time_(1) {pmat_ = &mat_; rhs_=0.0;}

  // return the number of degrees of freedom
  // of the problem
  long get_size() const {return nb_dof_;}
  long get_nb_dof() const {return nb_dof_;}
  long get_nb_dof_time() const {return nb_dof_time_;}

  void set_map(std::map<std::pair<int,int>,R>& A) {
    pmat_ = &A;
  }
  void set_map() {
    pmat_ = &mat_;
  }
  // void set_to_built_in_map() {
  //   pmat_ = &mat_;
  // }
  void cleanMatrix() {
    mat_.clear();
  }
  void init(int n) { nb_dof_=n; rhs_.init(nb_dof_);}
  void init(int n, int nt) { nb_dof_=n; nb_dof_time_=nt; rhs_.init(nb_dof_);}


  void initIndex(const FElement& FKu, const FElement& FKv) {
    this->index_i0_ = mapIdx0_[&FKv.Vh];
    this->index_j0_ = mapIdx0_[&FKu.Vh];
  }
  Matrix& get_matrix() {return *pmat_;}
  Rn& get_rhs() {return rhs_;}
  Rn_ get_solution() {
    return Rn_(rhs_(SubArray(nb_dof_,0)));
  }


  R & operator()(int i, int j) { return (*pmat_)[std::make_pair(i+index_i0_,j+index_j0_)]; }
  R & operator()(int i) { return rhs_[i+index_i0_];}

 protected :
  double & addToLocalContribution(int i, int j) { return local_contribution_matrix_[std::make_pair(i+index_i0_,j+index_j0_)]; }
  void addLocalContribution() {
    this->index_i0_ = 0;
    this->index_j0_ = 0;
    for (auto q=local_contribution_matrix_.begin(); q != this->local_contribution_matrix_.end(); ++q) {
      (*this)(q->first.first,q->first.second) += q->second;
    }
    this->local_contribution_matrix_.clear();
  }
  void addLocalContributionLagrange(int nend) {
    this->index_j0_ = 0;
    this->index_i0_ = 0;
    for (auto q=local_contribution_matrix_.begin(); q != this->local_contribution_matrix_.end(); ++q) {
      (*this)(q->first.first, nend) += q->second;
      (*this)(nend, q->first.first) += q->second;
    }
    this->local_contribution_matrix_.clear();
  }



public:

  void applyPreconditioning(std::map<std::pair<int,int>,R>& P) {
    int N = nb_dof_;
    SparseMatrixRC<double> Pl(N,N,P);
    {// P*A -> DF
      SparseMatrixRC<double> A(N,N,mat_);
      multiply(Pl, A, mat_);
    }
    {// (A*P)* -> DF
      SparseMatrixRC<double> A(N,N,mat_);
      multiply(A, Pl, mat_);
    }

    Rn x(N, 0.);

    multiply(N, N, P, rhs_, x);

    rhs_.resize(N);

    rhs_ = x;

  }
  void recoverSolution(std::map<std::pair<int,int>,R>& P) {
    int N = nb_dof_;
    Rn x(N, 0.);
    multiply(N, N, P, rhs_, x);
    rhs_=x;
  }
  void leftPreconditioning(std::map<std::pair<int,int>,R>& P) {

    int N = nb_dof_;
    SparseMatrixRC<double> Pl(N,N,P);
    SparseMatrixRC<double> A(N,N,mat_);
    multiply(Pl, A, mat_);

    Rn x(N, 0.);

    multiply(N, N, P, rhs_, x);

    rhs_.resize(N);

    rhs_ = x;

  }


};



static ProblemOption defaultProblemOption;

template<typename Mesh>
class QuadratureOfProblem{
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;
  typedef QuadratureFormular1d QFT;
  typedef typename QFT::QP QPT;


  int order_space_element_quadrature_ = 5;
  int order_space_bord_quadrature_    = 5;

  const QFT * time_quadrature_formular = nullptr;

public:
  QuadratureOfProblem() {}
  QuadratureOfProblem(const ProblemOption& option) {
    order_space_element_quadrature_ = option.order_space_element_quadrature_;
    order_space_bord_quadrature_    = option.order_space_bord_quadrature_;
  }
  QuadratureOfProblem(const QFT& qt, const ProblemOption& option) : time_quadrature_formular(&qt){
    order_space_element_quadrature_ = option.order_space_element_quadrature_;
    order_space_bord_quadrature_    = option.order_space_bord_quadrature_;
  }

  const QF&  get_quadrature_formular_K()  const;
  const QFB& get_quadrature_formular_dK() const;

  const QF&  get_quadrature_formular_cutK()   const;
  const QFB& get_quadrature_formular_cutFace() const;

  int get_nb_quad_point_time() const {
    return (time_quadrature_formular)? time_quadrature_formular->n : 1;
  }
  QPT get_quadrature_time(int itq) const {
    assert(itq < get_nb_quad_point_time());
    return (time_quadrature_formular)? time_quadrature_formular->at(itq): QPT(1.,0.);
  }


};



#endif
