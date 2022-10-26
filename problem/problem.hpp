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
#include "mapping.hpp"
#include "../solver/solver.hpp"
#include "/usr/local/opt/libomp/include/omp.h"



struct CBorder {
  CBorder() {}
};
const  CBorder INTEGRAL_BOUNDARY;
struct CFacet {
  CFacet() {}
};
const CFacet INTEGRAL_INNER_FACET;
const CFacet INTEGRAL_INNER_EDGE_2D;
const CFacet INTEGRAL_INNER_FACE_3D;

struct CRidge {
  CRidge() {}
};
const CRidge INTEGRAL_INNER_RIDGE;
const CRidge INTEGRAL_INNER_NODE_2D;
const CRidge INTEGRAL_INNER_EDGE_3D;

struct CExtension {
  CExtension() {}
};
const CExtension INTEGRAL_EXTENSION;


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

  // // matrix is on a std::map form
  // Matrix mat_;

  // For openmp;
  const int thread_count_max_;
  std::vector<Matrix> mat_;
  std::vector<int> index_i0_ = {0};
  std::vector<int> index_j0_ = {0};
  std::vector<Matrix> local_contribution_matrix_;

  // pointer on a std::map
  // the user can give is own std::map
  // can be use for newton, matrix that wanna be saved by the user etc
  std::vector<Matrix*> pmat_;

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
  // int index_i0_ = 0, index_j0_ = 0;


public :
  ShapeOfProblem() : nb_dof_(0), nb_dof_time_(1)
  ,thread_count_max_ (1)
  {
     set_multithreading_tool();
  };
  ShapeOfProblem(int np) : nb_dof_(0), nb_dof_time_(1)
  ,thread_count_max_ (np)
  {
    set_multithreading_tool();
  }

private:
  void set_multithreading_tool() {
    // thread_count_max_ = np;
    mat_.resize(thread_count_max_);
    local_contribution_matrix_.resize(thread_count_max_);
    index_i0_.resize(thread_count_max_);
    index_j0_.resize(thread_count_max_);
    set_map();
  }


public:
  long get_size() const {return nb_dof_;}
  long get_nb_dof() const {return nb_dof_;}
  long get_nb_dof_time() const {return nb_dof_time_;}
  int get_nb_thread() const {return thread_count_max_;}

  void set_map(Matrix& A) {
    assert(thread_count_max_ == 1 && " There are multiple threads. You cannot set only one matrix");
    pmat_.resize(0);
    pmat_.push_back(&A);
  }
  void set_map(std::vector<Matrix>& M) {
    assert(thread_count_max_ <= M.size() && " There are multiple threads. Not good size of vector");
    pmat_.resize(0);
    for(auto& A : M) {
      pmat_.push_back(&A);
    }
  }
  void set_map() {
    pmat_.resize(0);
    for(auto& A : mat_) {
      pmat_.push_back(&A);
    }
  }
  void cleanBuildInMatrix() {
    for(auto& A : mat_) A.clear();
  }

  void init(int n) { nb_dof_=n; rhs_.init(nb_dof_);}
  void init(int n, int nt) { nb_dof_=n; nb_dof_time_=nt; rhs_.init(nb_dof_);}


  void initIndex(const FElement& FKu, const FElement& FKv) {
    int thread_id = omp_get_thread_num();
    this->index_i0_[thread_id] = mapIdx0_[&FKv.Vh];
    this->index_j0_[thread_id] = mapIdx0_[&FKu.Vh];
  }

  // Matrix& get_matrix() {return *(pmat_[0]);}
  // Rn& get_rhs() {return rhs_;}
  // Rn_ get_solution() {
  //   return Rn_(rhs_(SubArray(nb_dof_,0)));
  // }


  R & operator()(int i, int j) { return (*pmat_[0])[std::make_pair(i+index_i0_[0],j+index_j0_[0])]; }
  // R & operator()(int i, int j, int thread_id) { return (*pmat_[thread_id])[std::make_pair(i+index_i0_[thread_id],j+index_j0_[thread_id])]; }
  R & operator()(int i) { return rhs_[i+index_i0_[0]];}

 protected :
  double & addToLocalContribution(int i, int j) {
    return local_contribution_matrix_[0][std::make_pair(i+index_i0_[0],j+index_j0_[0])];
  }
  double & addToLocalContribution_Opt(int i, int j) {
    return  loc_mat(i,j);
  }
  double & addToLocalContribution(int i, int j, int thread_id) {
    return local_contribution_matrix_[thread_id][std::make_pair(i+index_i0_[thread_id],j+index_j0_[thread_id])];
    // return (*pmat_[thread_id])[std::make_pair(i+index_i0_[thread_id],j+index_j0_[thread_id])];

  }
  void addLocalContribution() {
    int thread_id = omp_get_thread_num();
    this->index_i0_[thread_id] = 0;
    this->index_j0_[thread_id] = 0;
    auto& A(*pmat_[thread_id]);
    // for (auto q=local_contribution_matrix_[thread_id].begin(); q != this->local_contribution_matrix_[thread_id].end(); ++q) {
      // (*this)(q->first.first,q->first.second, thread_id) += q->second;
      for (const auto& aij : local_contribution_matrix_[thread_id]) {

      A[std::make_pair(aij.first.first+index_i0_[thread_id],aij.first.second+index_j0_[thread_id])] += aij.second;
    }
    // std::cout << thread_id << "\t" << this->local_contribution_matrix_.size() << std::endl;
    // this->local_contribution_matrix_[thread_id].clear();
    local_contribution_matrix_[thread_id].clear();
  }
  // double & addLocalContribution(int i, int j, int id_thread=0) {
  //   return openmp_mat_[id_thread][std::make_pair(i+openmp_index_i0_[id_thread],j+openmp_index_j0_[id_thread])];
  // }
  void addLocalContribution_Opt(const FElement& FK) {
    for(int i = 0; i < FK.NbDoF(); ++i) {
      for(int j = 0; j < FK.NbDoF(); ++j) {
        (*this)(FK.loc2glb(i),FK.loc2glb(j)) += loc_mat(i,j);
      }
    }
    this->index_i0_[0] = 0;
    this->index_j0_[0] = 0;
  }
  void addLocalContributionLagrange(int nend) {
    this->index_j0_[0] = 0;
    this->index_i0_[0] = 0;
    for (auto q=local_contribution_matrix_[0].begin(); q != this->local_contribution_matrix_[0].end(); ++q) {
      (*this)(q->first.first, nend) += q->second;
      (*this)(nend, q->first.first) += q->second;
    }
    this->local_contribution_matrix_[0].clear();
  }



public:

  void applyPreconditioning(std::map<std::pair<int,int>,R>& P) {
    int N = nb_dof_;
    SparseMatrixRC<double> Pl(N,N,P);
    {// P*A -> DF
      SparseMatrixRC<double> A(N,N,mat_[0]);
      multiply(Pl, A, mat_[0]);
    }
    {// (A*P)* -> DF
      SparseMatrixRC<double> A(N,N,mat_[0]);
      multiply(A, Pl, mat_[0]);
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
    SparseMatrixRC<double> A(N,N,mat_[0]);
    multiply(Pl, A, mat_[0]);

    Rn x(N, 0.);

    multiply(N, N, P, rhs_, x);

    rhs_.resize(N);

    rhs_ = x;

  }
  void addMatMul(const KN_<R> & uuh) {
    assert(uuh.size()==nb_dof_);
    MatriceMap<double> A(nb_dof_, nb_dof_, mat_[0]);
    A.addMatMul(uuh, rhs_);
  }
  void gather_map() {
     for(int i=1;i<thread_count_max_;++i){
       auto& A(mat_[i]);
       for(const auto& aij : A) {
         mat_[0][aij.first] += aij.second;
       }
       A.clear();
     }
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
