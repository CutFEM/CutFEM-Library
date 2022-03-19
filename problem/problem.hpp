#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>
#include <map>
#include <list>
//#include "../parallel/cfmpi.hpp"
#include "../common/geometry.hpp"
#include "../common/marker.hpp"
#include "../common/timeInterface.hpp"
#include "../common/fracture.hpp"
#include "../common/GTime.hpp"
#include "../common/SparseMatMap.hpp"
#include "../common/DA.hpp"
#include "../FESpace/expression.hpp"
#include "../FESpace/interpolation.hpp"
// #include "../FESpace/interpolationCutFEM.hpp"
// #include "../FESpace/integration.hpp"
#include "../FESpace/restriction.hpp"
#include "../FESpace/CutFESpace.hpp"
#include "../FESpace/integrationFunFEM.hpp"
// #include "../toolFEM/TimeIntegration.hpp"
//#include "../num/OperationMV.hpp"
#include "../FESpace/paraview.hpp"
#include "../solver/solver.hpp"
#include "../num/matlab.hpp"
#include "../util/util.hpp"
#include "../util/cputime.h"
#include "finiteElement.hpp"
#include "macroElement.hpp"


struct CBorder {
  CBorder() {}
};
const  CBorder boundary;
struct CHyperFace {
  CHyperFace() {}
};
const CHyperFace innerEdge;
const CHyperFace innerFace;


class ShapeOfLinProblem {
public :

  Rn rhs;
  Ulint nDoF;
  std::map<std::pair<int,int>,R> mat;
public :
  ShapeOfLinProblem() : nDoF(0) {};
  ShapeOfLinProblem(int n) : nDoF(n), rhs(n) {rhs=0.0;}
  ShapeOfLinProblem(int n, int nn) : nDoF(n), rhs(nn) {rhs=0.0;}

  long size() {return nDoF;}

  R & operator()(int i, int j) { return mat[std::make_pair(i,j)]; }
  R & operator()(int i) { return rhs[i];}

};
class ShapeOfNonLinProblem {
  public :
  Ulint nDoF;
  Ulint nt;
  std::map<std::pair<int,int>,R> DF;     // the matrix
  std::map<std::pair<int,int>,R> NL;     // for optimizition when doing Newton
  std::map<std::pair<int,int>,R> *pmat;  // for optimization with Newton
  std::map<std::pair<int,int>,R> mapU0;  // (domain, indexInBackMesh)
  std::map<std::pair<int,int>,R> localContributionMatrix;
  Rn F;

  //Alias for connection
  std::map<std::pair<int,int>,R>& mat = DF;
  Rn& rhs = F;

  //index for sum of FESpace
  int index_i0 = 0, index_j0 = 0;
public :
  ShapeOfNonLinProblem() : nDoF(0), nt(0){ pmat = &DF;};
  ShapeOfNonLinProblem(long n) : nDoF(n),nt(0) {
    pmat = &DF;
    F.resize(n); F=0.;
  }

  long size() {return nDoF;}

  R & operator()(int i, int j) { return (*pmat)[std::make_pair(i+index_i0,j+index_j0)]; }
  R & addToLocalContribution(int i, int j) { return localContributionMatrix[std::make_pair(i+index_i0,j+index_j0)]; }

  R & operator()(int i) { return F[i+index_i0];}

  void addMatMul(const KN_<R> & uuh) {
    assert(uuh.size()==nDoF);
    MatriceMap<double> A(nDoF, nDoF, DF);
    A.addMatMul(uuh, F);
  }

  void setMatrixTo(std::map<std::pair<int,int>,R>& A) {
    pmat = &A;
  }
  void reset() {
    pmat = &DF;
  }
  void saveMatrix() {
    pmat = &NL;
  }
  void cleanMatrix() {
    DF.clear(); NL.clear();
    pmat = &DF;
  }
  void resetIndex() {
    this->index_i0 = 0;
    this->index_j0 = 0;
  }

  void addMapToMap(){
    for (typename map<pair<int,int>,R>::const_iterator q=DF.begin(); q != DF.end(); ++q) {
      NL[std::make_pair(q->first.first,q->first.second)] += q->second;
    }
  }
  void addLocalContribution() {
    resetIndex();
    for (auto q=localContributionMatrix.begin(); q != this->localContributionMatrix.end(); ++q) {
      (*this)(q->first.first,q->first.second) += q->second;
    }
    this->localContributionMatrix.clear();

  }


  void preconditionning(std::map<std::pair<int,int>,R>& P) {

    int N = nDoF;
    SparseMatrixRC<double> Pl(N,N,P);
    SparseMatrixRC<double> A(N,N,DF);
    multiply(Pl, A, DF);

    Rn x(N, 0.);

    multiply(N, N, P, rhs, x);

    rhs.resize(N);

    rhs = x;

  }


};


// Base class for problem.
// contain info about the linear system
class ShapeOfProblem {
  protected :

  // The right hand side vector
  KN<double> rhs_;

  // matrix is on a std::map form
  std::map<std::pair<int,int>,R>  mat_;

  // pointer on a std::map
  // the user can give is own std::map
  // can be use for newton, matrix that wanna be saved by the user etc
  std::map<std::pair<int,int>,R> *pmat_;

  // std::map used to save CutFEM solution on the background mesh
  // for time dependent problem
  // std::pair(domain, dof_on backSpace) => value
  std::map<std::pair<int,int>,R> mapU0_;

  // local map matrix
  // to reduce the numer of access to the global std::map
  // when the integral is computed on an elements_to_integrate
  // the local contribution is added to the global matrix
  // pointed by pmat
  std::map<std::pair<int,int>,R> local_contribution_matrix_;

  // number of degree of freedom of the problem
  // This is never modify after initialization
  // => not when adding lagrange multiplier
  Ulint nb_dof_;

  //index where degree of freedom of consider space starts
  // this make possible to use as many space as we need
  int index_i0_ = 0, index_j0_ = 0;

public :
  ShapeOfProblem() : nb_dof_(0) {};
  ShapeOfProblem(int n) : nb_dof_(n), rhs_(n) {rhs_=0.0;}


  // return the number of degrees of freedom
  // of the problem
  long get_size() const {return nb_dof_;}

  void set_map(std::map<std::pair<int,int>,R>& A) {
    pmat_ = &A;
  }
  void set_to_buildin_map() {
    pmat_ = &mat_;
  }
  void CleanMatrix() {
    mat_.clear();
  }


protected :
  double & operator()(int i, int j) { return (*pmat_)[std::make_pair(i,j)]; }
  double & addToLocalContribution(int i, int j) { return local_contribution_matrix_[std::make_pair(i+index_i0_,j+index_j0_)]; }
  double & operator()(int i) { return rhs_[i];}


};



#endif
