#include <map>
#include <cassert>
#include "RNM.hpp"
using namespace std;
#include "../util/util.hpp"
#include <map>
#include <set>
#include <iostream>
#include <iterator>
#include <algorithm>



void buil_CSR_array(int n,const  map<pair<int,int>,R> & M, int32_t* p, int32_t* j, double * a) ;
void multiply( int N, const std::map<std::pair<int,int>,double>& A, const std::map<std::pair<int,int>,double>& B, std::map<std::pair<int,int>,double>& C);
void eraseAndSetRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, std::map<int, double>& dof2rm);
void eraseAndSetRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, int , int, double);
void eraseRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, std::set<int>& dof2rm);

// void multiply( int N, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b);
void multiply( int N, int M, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b);


template<class R>
struct  VirtualMatrice { public:
    int N,M;
    VirtualMatrice(int nn,int mm): N(nn),M(mm) {}
    VirtualMatrice(int nn): N(nn),M(nn) {}
  virtual void addMatMul(const KN_<R> &  x, KN_<R> & y) const =0;
  virtual ~VirtualMatrice(){}
};




template<class R>
class SparseMatrixRC: public VirtualMatrice<R> {
public:
  int n,m,nbcoef;
  int32_t *p;
  int32_t *j;
  double  *a;
  SparseMatrixRC(int nn,int mm,const  map<pair<int,int>,R> & M);

  R *pcoef(int ii,int jj);
  const R & operator()(int i,int j) const {R *p=pcoef(i,j);assert(p);return *p;}
  R & operator()(int i,int j) {R *p=pcoef(i,j);assert(p);return *p;}
  R & operator[](pair<int,int> ij) {R *p=pcoef(ij.first,ij.second);assert(p);return *p;}
  const R & operator[](pair<int,int> ij) const  {R *p=pcoef(ij.first,ij.second);assert(p);return *p;}

  void getRow(int i, Rn& v) const {
    int n = p[i+1] - p[i];
    v.resize(n);
    for(int k=p[i], j=0; k<p[i+1];++k, ++j) {
      v[j] = a[k];
    }
  }
  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {
    for(int i=0;i<n;++i)
    for(int k=p[i];k<p[i+1];++k)
    Ax[i] += a[k]*x[j[k]];
  }
  ~SparseMatrixRC() { delete [] p; delete [] j; delete [] a;}
private:
  SparseMatrixRC(const SparseMatrixRC &);
  void operator=(const SparseMatrixRC &);
};

void multiply(const SparseMatrixRC<double>& A, const SparseMatrixRC<double>& B, std::map<std::pair<int,int>,double>& C);


template<class R>
class MatriceMap: VirtualMatrice<R> { public:
    typedef map<pair<int,int>,R> Map;

    const Map &m;
    MatriceMap(int in,int im,const Map & mm) :   VirtualMatrice<R>(in,im), m(mm) {}
    void addMatMul(const KN_<R> & x, KN_<R> & Ax) const ;
};

template<class R>
KN<R> operator*(const map<pair<int,int>,R>& A, const KN_<R> & x) {
  KN<R> Ax(x.size());
  auto last_index = A.end();
  last_index--;
  int row_end = last_index->first.first;
  assert(row_end < x.size());
  for (const auto& aij : A) {
      int i= aij.first.first;
      int j= aij.first.second;
      Ax[i] += aij.second*x[j];
  }
  return Ax;
}

template<class R>
void MatriceMap<R>::addMatMul(const KN_<R> & x, KN_<R> & Ax) const {
    for (typename Map::const_iterator k=this->m.begin(); k != this->m.end(); ++k)
    {
        int i= k->first.first;
        int j= k->first.second;
        R  aij= k->second;
        Ax[i] += aij*x[j];
    }
}


template<class R>
SparseMatrixRC<R>::SparseMatrixRC(int nn,int mm,const  map<pair<int,int>,R> & M)
:    VirtualMatrice<R>(nn,mm),
n(nn), m(mm),
nbcoef(M.size()),
p(new int[nn+1]),
j(new int[nbcoef]),
a(new R[nbcoef]) {
  R cmm=0;

  int k=0;
  for(int i=0;i<=n;++i)
  p[i]=0; // pour les lignes vide
  for (typename map<pair<int,int>,R>::const_iterator q=M.begin(); q != M.end(); ++q, ++k) {
    int i =q->first.first;
    p[i+1]=k+1;
    this->j[k]=q->first.second;
    this->a[k]=q->second;
    cmm=max(cmm,this->a[k]);
  }
  for(int i=1;i<=n;++i)
  p[i] = max(p[i-1],p[i]); // pour les lignes vides
  assert(k==this->nbcoef);
}



template<class R>
SparseMatrixRC<R>::SparseMatrixRC(const SparseMatrixRC<R> & A)
:   VirtualMatrice<R>(A.n,A.m),
n(A.n),
m(A.m),
nbcoef(A.nbcoef),
p(new int[n+1]),
j(new int[nbcoef]),
a(new R[nbcoef])
{
    R cmm=0;

    for(int i=0;i<=n;++i)
        p[i]=0; // pour les lignes vide
    for (int k=0;k<nbcoef;++k)
    {
        int i =A.i[k];
        p[i+1]=k+1;
        this->j[k]=A.j[k];
        this->a[k]=A.a[k];
        cmm=max(cmm,this->a[k]);
    }
    for(int i=1;i<=n;++i)
        p[i] = max(p[i-1],p[i]); // pour les lignes vides
    cout <<  " nb coef = " << nbcoef << " c max = "<< cmm << endl;
}

template<class R>
R * SparseMatrixRC<R>::pcoef(int ii,int jj) {
    if(ii<0 && jj>=n) return 0;
    int k0=p[ii],k1=p[ii+1]-1;
    while (k0<=k1)
    {
        int km=(k0+k1)/2;
        int aa=0;
        if ( j[km] > jj)
            aa=-1;
        else if (j[km]<jj)
            aa=1;
        if(aa<0) k1=km-1;
        else if (aa>0) k0=km+1;
        else { return a+km; }
    }
    return 0;
}
