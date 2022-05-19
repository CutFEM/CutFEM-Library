#include <map>
#include <cassert>
#include "RNM.hpp"
using namespace std;
#include "HeapSort.hpp"
#include "../util/util.hpp"
#include <map>
#include <set>
#include <iostream>
#include <iterator>
#include <algorithm>




void multiply( int N, const std::map<std::pair<int,int>,double>& A, const std::map<std::pair<int,int>,double>& B, std::map<std::pair<int,int>,double>& C);
void eraseRow( int N, std::map<std::pair<int,int>,double>& A, Rn& b, std::map<int, double>& dof2rm);
// void multiply( int N, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b);
void multiply( int N, int M, const std::map<std::pair<int,int>,double>& A, const Rn& rhs, Rn& b);



// pour afficher une pair
template<typename A,typename B>
ostream & operator<<(ostream & f, pair<A,B> ab)
{
    f << ab.first << " : " <<  ab.second ;
    return f;
}

// --------------------------------------//
// pour afficher un container de la stl  //
// --------------------------------------//
template<typename T>
void  show(char * s,const T & l,const char * separateur="\n")
{
    cout << s <<  separateur;
    for (typename T::const_iterator i=l.begin(); i != l.end(); ++i)
        cout << '\t' << * i <<  separateur;
    cout << endl;
}


#ifndef InternalError
typedef void (* TypeofInternalErrorRoutine)(const char *) ;
static TypeofInternalErrorRoutine &InternalErrorRoutinePtr()
{
  static TypeofInternalErrorRoutine routine=0;
  return routine;
}

static void InternalError(const char * str) {
  if (InternalErrorRoutinePtr() ) (*InternalErrorRoutinePtr())(str);
  cerr << str;
  exit(1);
}
inline void  SetInternalErrorRoutine(TypeofInternalErrorRoutine f)
{
  InternalErrorRoutinePtr()=f;
}
#endif


template<class R>
struct  VirtualMatrice { public:
    int N,M;
    VirtualMatrice(int nn,int mm): N(nn),M(mm) {}
    VirtualMatrice(int nn): N(nn),M(nn) {}
  //  y += A x
  virtual void addMatMul(const KN_<R> &  x, KN_<R> & y) const =0;
  virtual void addMatTransMul(const KN_<R> &  , KN_<R> & ) const
    { InternalError("VirtualMatrice::addMatTransMul not implemented "); }
  virtual bool WithSolver() const {return false;} // by default no solver
  virtual void Solve( KN_<R> &  ,const KN_<R> & ) const
    { InternalError("VirtualMatrice::solve not implemented "); }

#ifdef VersionFreeFempp
  virtual bool ChecknbLine  (int n) const= 0;
  virtual bool ChecknbColumn  (int m) const =0;
#else
  virtual bool ChecknbLine  (int n) const {return true;}
  virtual bool ChecknbColumn  (int m) const {return true;}
#endif
  struct  plusAx { const VirtualMatrice * A; const KN_<R>   x;
   plusAx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y)
      { ffassert(B->ChecknbColumn(y.N())); }
    };

   plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);}

  struct  plusAtx { const VirtualMatrice * A; const KN_<R>   x;
   plusAtx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y)
    {ffassert(B->ChecknbLine(y.N()));} };

  struct  solveAxeqb { const VirtualMatrice * A; const KN_<R>   b;
   solveAxeqb( const VirtualMatrice * B,const KN_<R> &  y) :A(B),b(y)
    {ffassert(B->ChecknbColumn(y.N()));} };

  virtual ~VirtualMatrice(){}
};




template<class R>
class SparseMatrix: public VirtualMatrice<R> {
public:
    //  data struct  k = 0, nbcoef -1 :   a_ij = a[k] i=i[k], j=j[k]
    int n,m,nbcoef;
    int *i;
    int *j;
    R *a;
    template<class Mesh>  SparseMatrix(Mesh & Th);
    SparseMatrix(int nn,int mm, map<pair<int,int>,R> M);

    SparseMatrix();

    R *pcoef(int ii,int jj);
    const R & operator()(int i,int j) const {R *p=pcoef(i,j);assert(p);return *p;}
    R & operator()(int i,int j) {R *p=pcoef(i,j);assert(p);return *p;}
    R & operator[](pair<int,int> ij) {R *p=pcoef(ij.first,ij.second);assert(p);return *p;}
    const R & operator[](pair<int,int> ij) const  {R *p=pcoef(ij.first,ij.second);assert(p);return *p;}

    typedef  typename VirtualMatrice<R>::plusAx plusAx;
    void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const
    { for(int k=0;k<nbcoef;++k)
        Ax[i[k]] += a[k]*x[j[k]];
    }
    //---mon petit ajout-----
    void addMatMulT(const  KN_<R>  & x, KN_<R> & Ax) const
    { for(int k=0;k<nbcoef;++k)
        Ax[j[k]] += a[k]*x[i[k]];
    }
    //--------------------------

    plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
    ~SparseMatrix() { delete [] i; delete [] j; delete [] a;}
private:
    SparseMatrix(const SparseMatrix &);
    void operator=(const SparseMatrix &);
};

// Constructeur par défaut


template<class R>
SparseMatrix<R>::SparseMatrix()
:   VirtualMatrice<R>(1,1),
n(1),
m(1),
nbcoef(0),
i(NULL),
j(NULL),
a(NULL)
{}

// UMFPACK data struct  transp ....
template<class R>
class SparseMatrixRC: public VirtualMatrice<R> {
public:
  int n,m,nbcoef;
  int *p;
  int *j;
  R *a;
  SparseMatrixRC(int nn,int mm,const  map<pair<int,int>,R> & M);
  SparseMatrixRC(const SparseMatrix<R> & M);
  SparseMatrixRC(int nn,int mm,const  Rnm & M);
  // SparseMatrixRC(int nn,int mm, KN_<int>& pp, KN_<int>& jj, KN_<double>& aa);


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
    typedef  typename VirtualMatrice<R>::plusAx plusAx;
    void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const
    {
        for(int i=0;i<n;++i)
            for(int k=p[i];k<p[i+1];++k)
                Ax[i] += a[k]*x[j[k]];
    }
    plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
    ~SparseMatrixRC() { delete [] p; delete [] j; delete [] a;}
private:
    SparseMatrixRC(const SparseMatrixRC &);
    void operator=(const SparseMatrixRC &);
};

void multiply(const SparseMatrixRC<double>& A, const SparseMatrixRC<double>& B, std::map<std::pair<int,int>,double>& C);


// template<class R>
// class MatrixDiag : public VirtualMatrice<R> { public:
//     KN_<R> d;
//     MatrixDiag(const KN_<R> dd) :   VirtualMatrice<R>(dd.N(),dd.N()), d(dd) {}
//     typedef typename VirtualMatrice<R>::plusAx plusAx;
//     void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const
//     { Ax += DotStar_KN_<R>(x,d);  }
//     plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
// };

template<class R>
class MatriceMap: VirtualMatrice<R> { public:
    typedef  typename  VirtualMatrice<R>::plusAx plusAx;
    typedef map<pair<int,int>,R> Map;

    const Map &m;
    MatriceMap(int in,int im,const Map & mm) :   VirtualMatrice<R>(in,im), m(mm) {}
    void addMatMul(const KN_<R> & x, KN_<R> & Ax) const ;
    plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};

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
template<class Mesh>
SparseMatrix<R>::SparseMatrix(Mesh & Th) :VirtualMatrice<R>(Th.nv,Th.nv)
{
    const int nve = Mesh::Element::nv;
    int end=-1;
    int nn= Th.nt*nve;
    int mm= Th.nv;
    KN<int> head(mm);
    KN<int> next(nn);
    KN<int> mark(mm);
    int color=0;
    mark=color;

    head = end;
    for (int p=0;p<nn;++p)
    {
        int s= Th(p/nve,p%nve);
        next[p]=head[s];
        head[s]=p;
    }
    nbcoef =0;
    n=mm;
    m=mm;
    int kk=0;
    for (int step=0;step<2;++step) //   2 etapes une pour le calcul du nombre de coef
        // l'autre pour la construction
    {
        for (int ii=0;ii<mm;++ii)
        {
            color++;
            int ki=nbcoef;
            for (int p=head[ii];p!=end;p=next[p])
            {
                int k=p/nve;
                for (int l=0;l<nve;++l)
                {
                    int jj= Th(k,l);
                    if( mark[jj] != color ) // c'est un nouveau sommet de l'ensemble
                    {
                        if (step==1)
                        {
                            i[nbcoef]=ii;
                            j[nbcoef]=jj;
                            a[nbcoef++]=R();
                        }
                        else
                            nbcoef++;
                    }
                    mark[jj]=color; // on colorie le sommet j;
                }
            }
            int ki1=nbcoef;
            if(step==1)
                HeapSort(j+ki,ki1-ki);
        }

        if(step==0)
        {
            cout << " Allocation des tableaux  " << nbcoef << endl;
            i= new int[nbcoef];
            j= new int[nbcoef];
            a= new R[nbcoef];
            kk=nbcoef;
            nbcoef=0;
        }


    }
    /*
     for (int k=0;k<nbcoef;k++)
     cout<< k << "  :: " << i[k] << " " << j[k] << " " << endl;
     */
}

template<class R>
SparseMatrix<R>::SparseMatrix(int nn,int mm, map<pair<int,int>,R> M)
:   VirtualMatrice<R>(nn,mm),
n(nn),
m(mm),
nbcoef(M.size()),
i(new int[nbcoef]),
j(new int[nbcoef]),
a(new R[nbcoef])
{
    R cmm=0;

    int k=0;
    for (typename map<pair<int,int>,R>::const_iterator p=M.begin(); p != M.end(); ++p, ++k)
	{
        this->i[k]=p->first.first;
        this->j[k]=p->first.second;
        this->a[k]=p->second;
        //assert(i[k]<n && j[k] <m && i[k]>=0 && j[k]>=0);
        cmm=max(cmm,this->a[k]);
	}
    //cout <<  " nb coef = " << nbcoef << " c max = "<< cmm << endl;
    assert(k==this->nbcoef);
}

template<class R>
SparseMatrixRC<R>::SparseMatrixRC(int nn,int mm,const  map<pair<int,int>,R> & M)
:    VirtualMatrice<R>(nn,mm),
n(nn),
m(mm),
nbcoef(M.size()),
p(new int[nn+1]),
j(new int[nbcoef]),
a(new R[nbcoef])
{
    R cmm=0;

    int k=0;
    for(int i=0;i<=n;++i)
        p[i]=0; // pour les lignes vide
    for (typename map<pair<int,int>,R>::const_iterator q=M.begin(); q != M.end(); ++q, ++k)
    {
        int i =q->first.first;
	p[i+1]=k+1;
        this->j[k]=q->first.second;
        this->a[k]=q->second;
        //assert(i[k]<n && j[k] <m && i[k]>=0 && j[k]>=0);
        cmm=max(cmm,this->a[k]);
    }
    for(int i=1;i<=n;++i)
        p[i] = max(p[i-1],p[i]); // pour les lignes vides
    // cout <<  " nb coef = " << nbcoef << " c max = "<< cmm << endl;
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
SparseMatrixRC<R>::SparseMatrixRC(int nn,int mm, const Rnm & M)
  :    VirtualMatrice<R>(nn,mm),
       n(nn),
       m(mm),
       nbcoef(M.size()),
       p(new int[nn+1]),
       j(new int[nbcoef]),
       a(new R[nbcoef])
{
  R cmm=0;

  int k=0;
  for(int i=0;i<=n;++i)
    p[i]=0; // pour les lignes vide
  for (int ii = 0; ii < n; ++ii) {
    for (int jj = 0; jj < n; ++jj, ++k)
    {
      p[ii+1]=k+1;
      this->j[k]=jj;
      this->a[k]=M(ii,jj);
      cmm=std::max(cmm,this->a[k]);
    }
  }
  for(int i=1;i<=n;++i)
    p[i] = std::max(p[i-1],p[i]); // pour les lignes vides
  assert(k==this->nbcoef);
}

// SparseMatrixRC(int nn,int mm, KN_<int>& pp, KN_<int>& jj, KN_<double>& aa)
// :   VirtualMatrice<R>(nn,mm),
//     n(nn),
//     m(mm),
//     nbcoef(aa.size()),
//     p(pp),
//     j(jj),
//     a(aa) {}
//
//


template<class R>
R * SparseMatrix<R>::pcoef(int ii,int jj)
{
    if(ii<0 && jj>=n) return 0;
    int k0=0,k1=nbcoef-1;
    while (k0<=k1)
    {
        int km=(k0+k1)/2;
        int aa=0;
        if ( i[km] > ii)
            aa=-1;
        else if (i[km]<ii)
            aa=1;
        else  if ( j[km] > jj)
            aa=-1;
        else if (j[km]<jj)
            aa=1;
        if(aa<0) k1=km-1;
        else if (aa>0) k0=km+1;
        else { return a+km; }
    }
    return 0;
}
template<class R>
R * SparseMatrixRC<R>::pcoef(int ii,int jj)
{
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


// void addMapToMap(const std::map<std::pair<int,int>,R> & A,
//                        std::map<std::pair<int,int>,R> & NL);
                       // {
  //
  // for (typename map<pair<int,int>,R>::const_iterator q=A.begin(); q != A.end(); ++q) {
  //
  //   NL[std::make_pair(q->first.first,q->first.second)] += q->second;
  //
  // }
// }
