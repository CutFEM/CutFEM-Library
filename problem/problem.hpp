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
#include "../num/paraview.hpp"
#include "../solver/solver.hpp"
#include "../num/matlab.hpp"
#include "../util/util.hpp"
#include "../util/cputime.h"
#include "finiteElement.hpp"
#include "macroElement.hpp"


struct CBorder {
  CBorder() {}
};
const CBorder boundary;
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

};




#endif
