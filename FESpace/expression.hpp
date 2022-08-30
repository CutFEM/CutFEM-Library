#ifndef _EXPRESSION_HPP
#define _EXPRESSION_HPP
#include "FESpace.hpp"
#include <list>

struct Normal_Component {
  virtual int component() const =0;
};
struct Normal_Component_X : public Normal_Component {
  virtual int component() const  {return 0;}
};
struct Normal_Component_Y : public Normal_Component {
  virtual int component() const {return 1;}
};
struct Normal_Component_Z : public Normal_Component {
  virtual int component() const {return 2;}
};

struct BaseVector{};
struct Normal   : public BaseVector {
  Normal_Component_X x;
  Normal_Component_Y y;
  Normal_Component_Z z;
  // int operator[](int i) const {return idx[i];}
};
struct Conormal : public BaseVector {
  Normal_Component_X x;
  Normal_Component_Y y;
  Normal_Component_Z z;
};
struct Tangent {
  // static const int idx[3];
  // int operator[](int i) const {return idx[i];}
};
struct Projection {
  Normal normal;
  KN<int> operator()(int i, int j) const {
    KN<int> ar(2); ar(0) = i; ar(1) = j;
    return ar;
  }
};

class ExpressionVirtual;
template<typename M> class ExpressionFunFEM;


template<typename M>
class FunFEM : public FunFEMVirtual {
public:
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Rd Rd;

// private:
// double * data = nullptr;
// KN_<double> v;
// bool alloc = false;
// double* databf = nullptr;
bool alloc = false;
double* databf = nullptr;
// using FunFEMVirtual::data;


public:
  FESpace const *  Vh = nullptr;
  TimeSlab const * In = nullptr;
  public :
  // FunFEM() : v(data, 0) {}
  FunFEM() : FunFEMVirtual () {}
  FunFEM(const FESpace& vh) : FunFEMVirtual (vh.NbDoF()),
  Vh(&vh), alloc(true) ,databf(new double[10*vh[0].NbDoF()*vh.N*4]){}
  // : Vh(&vh), data(new double[vh.NbDoF()]), v(data,vh.NbDoF()), alloc(true) ,databf(new double[10*vh[0].NbDoF()*vh.N*4]){v=0.;}
  FunFEM(const FESpace& vh, const TimeSlab& in) : FunFEMVirtual (vh.NbDoF()*in.NbDoF()),
  Vh(&vh), In(&in), alloc(true) ,databf(new double[10*vh[0].NbDoF()*vh.N*4]){}
  // : Vh(&vh), In(&in), data(new double[vh.NbDoF()*in.NbDoF()]), v(data,vh.NbDoF()*in.NbDoF()), alloc(true) ,databf(new double[10*vh[0].NbDoF()*vh.N*4]){v=0.;}
  FunFEM(const FESpace& vh, Rn_& u) :FunFEMVirtual (u),
  alloc(false), Vh(&vh), databf(new double[10*vh[0].NbDoF()*vh.N*4]) { }
  // : v(u), alloc(false), Vh(&vh), databf(new double[10*vh[0].NbDoF()*vh.N*4]) { }
  FunFEM(const FESpace& vh, const TimeSlab& in, Rn_& u) : FunFEMVirtual (u),
   alloc(false), Vh(&vh), In(&in), databf(new double[10*vh[0].NbDoF()*vh.N*4]) {}
  // : v(u), alloc(false), Vh(&vh), In(&in), databf(new double[10*vh[0].NbDoF()*vh.N*4]) {}

  FunFEM(const FESpace& vh, R (*f)(const Rd, int i) )
  : FunFEMVirtual (vh.NbDoF()),
  alloc(true), Vh(&vh),
    // data(new double[vh.NbDoF()]),
    // v(data, vh.NbDoF()),
    databf(new double[10*vh[0].NbDoF()*vh.N*4])
     {
       interpolate(*Vh, this->v, f);
     }
  FunFEM(const FESpace& vh, R (*f)(const Rd, int i, int dd) )
  : FunFEMVirtual (vh.NbDoF()),
  alloc(true), Vh(&vh),
    // data(new double[vh.NbDoF()]),
    // v(data, vh.NbDoF()),
    databf(new double[10*vh[0].NbDoF()*vh.N*4])
     {
       interpolate(*Vh, this->v, f);
     }
  FunFEM(const FESpace& vh, R (*f)(const Rd, int i, int dd, R tt), R tid )
  : FunFEMVirtual (vh.NbDoF()),
  alloc(true), Vh(&vh),
  // data(new double[vh.NbDoF()]),
  // v(data, vh.NbDoF()) ,
  databf(new double[10*vh[0].NbDoF()*vh.N*4])
  {
    interpolate(*Vh, this->v, f, tid);
  }

  FunFEM(const FESpace& vh, R (*f)(const Rd, int i, R tt), R tid )
  : FunFEMVirtual (vh.NbDoF()),
  alloc(true), Vh(&vh),
    // data(new double[vh.NbDoF()]),
    // v(data, vh.NbDoF()) ,
    databf(new double[10*vh[0].NbDoF()*vh.N*4])
    {
      interpolate(*Vh, this->v, f, tid);
    }

    FunFEM(const FESpace& vh, const TimeSlab& in, R (*f)(const Rd, int i, R tt) )
    : FunFEMVirtual (vh.NbDoF()*in.NbDoF()),
    alloc(true), Vh(&vh), In(&in),
      // data(new double[vh.NbDoF()*in.NbDoF()]),
      // v(data, vh.NbDoF()*in.NbDoF()) ,
      databf(new double[10*vh[0].NbDoF()*vh.N*4])
      {
        interpolate(*Vh,*In, this->v, f);
      }

      FunFEM(const FESpace& vh, const ExpressionVirtual& fh);
      FunFEM(const FESpace& vh, const ExpressionVirtual& fh1, const ExpressionVirtual& fh2);


    void init(const FESpace& vh) {
      assert(!data);
      Vh = &vh;
      data = new double[Vh->NbDoF()];
      v.set(data, Vh->NbDoF());
      alloc = true;
      v = 0.;
      if(!databf) databf = new double[10*vh[0].NbDoF()*vh.N*4];
    }
    void init(const KN_<R> & a) {
      assert(v.size() == a.size());
      for(int i=0;i<v.size();++i) data[i] = a(i);
    }

    void init(const FESpace& vh, R (*f)(const Rd, int i)) {
      assert(!data);
      Vh = &vh;
      data = new double[Vh->NbDoF()];
      v.set(data, Vh->NbDoF());
      alloc = true;
      interpolate(*Vh, v, f);

      if(!databf) databf = new double[10*vh[0].NbDoF()*vh.N*4];
    }

    void init(const FESpace& vh, R (*f)(const Rd, int i,int doma)) {
      assert(!data);
      Vh = &vh;
      data = new double[Vh->NbDoF()];
      v.set(data, Vh->NbDoF());
      alloc = true;
      interpolate(*Vh, v, f);
      if(!databf) databf = new double[10*vh[0].NbDoF()*vh.N*4];

    }

    void init(const FESpace& vh, R (*f)(const Rd, int, R), R tid) {
      // assert(!data);
      if(data) delete [] data;
      Vh = &vh;
      data = new double[Vh->NbDoF()];
      v.set(data, Vh->NbDoF());
      alloc = true;
      interpolate(*Vh, v, f, tid);
      if(!databf) databf = new double[10*vh[0].NbDoF()*vh.N*4];
    }

    void init(const FESpace& vh, const TimeSlab& in, R (*f)(const Rd, int, R)) {
      // assert(!data);
      if(data) delete [] data;
      Vh = &vh;
      In = &in;
      data = new double[Vh->NbDoF()*In->NbDoF()];
      v.set(data, Vh->NbDoF()*In->NbDoF());
      alloc = true;
      interpolate(*Vh,*In, v, f);
      if(!databf) databf = new double[10*vh[0].NbDoF()*vh.N*4];
    }

    void swap(FunFEM& f) {
      assert(v.size() == f.v.size());
      double * temp = data;
      data = f.data;
      f.data = temp;
      v.set(data, Vh->NbDoF());
      f.v.set(f.data, f.Vh->NbDoF());
    }

  double& operator()(int i) { return v(i); }
  double operator()(int i) const { return v(i); }
  operator Rn() const {return Rn(v);}


  const KN_<double>& getArray() const {return v;}

  double eval(const int k, const R* x,            int cu=0, int op=0) const ;
  double eval(const int k, const R* x, const R t, int cu=0, int op=0, int opt=0) const ;
  void eval ( R* u, const int k) const ;

  double evalOnBackMesh(const int kb, int dom, const R* x,            int cu, int op         ) const ;
  double evalOnBackMesh(const int kb, int dom, const R* x, const R t, int cu, int op, int opt ) const ;

  int size() const { return Vh->NbDoF()*((In)?In->NbDoF():1);}
  int size(int k) const { return (*Vh)[k].NbDoF();}
  void print() const;
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return Vh->idxElementFromBackMesh(kb, dd);
  }
  std::vector<int> getAllDomainId(int k) const {
    return Vh->getAllDomainId(k);
  };

  const FESpace& getSpace() const { return *Vh;}
  BasisFctType getBasisFctType() const {return Vh->basisFctType;}
  int getPolynomialOrder() const {return Vh->polynomialOrder;}
  std::list<ExpressionFunFEM<M>> expression(int n = -1)const;
  std::list<ExpressionFunFEM<M>> expression(int n, int i0)const;

  ~FunFEM(){
    if(databf) delete [] databf;
    if(alloc) delete [] data; }

private:
  FunFEM(const FunFEM& f);
  void operator=(const FunFEM& f) ;
};


class ParameterCutFEM;
template<typename M>
class ExpressionFunFEM;
class ExpressionProduct;
class ExpressionDif;

template<typename M>
ExpressionFunFEM<M> dx(const ExpressionFunFEM<M>& u);
template<typename M>
ExpressionFunFEM<M> dy(const ExpressionFunFEM<M>& u);
template<typename M>
ExpressionFunFEM<M> dz(const ExpressionFunFEM<M>& u);
template<typename M>
ExpressionFunFEM<M> dt(const ExpressionFunFEM<M>& u);

class ExpressionVirtual {
public:
  int cu, op, opt;
  KN<int> ar_normal;
  int domain = -1;

  ExpressionVirtual() : cu(0) , op(0), opt(0) {}
  ExpressionVirtual(int cc, int opp) : cu(cc) , op(opp), opt(0) {}
  ExpressionVirtual(int cc, int opp, int oppt) : cu(cc) , op(opp), opt(oppt) {}
  ExpressionVirtual(int cc, int opp, int oppt, int dom) : cu(cc) , op(opp), opt(oppt), domain(dom){}
  virtual R operator()(long i) const = 0;
  virtual R eval(const int k, const R* x, const R* normal = nullptr) const = 0;
  virtual R eval(const int k, const R* x, const R t, const R* normal = nullptr) const = 0;
public:
  virtual R evalOnBackMesh(const int k, int dom, const R* x,            const R* normal = nullptr) const = 0;
  virtual R evalOnBackMesh(const int k, int dom, const R* x, const R t, const R* normal = nullptr) const = 0;
  virtual int idxElementFromBackMesh(int kb, int dd=0) const = 0;
  virtual std::vector<int> getAllDomainId(int k) const {assert(0); };
public:
  R GevalOnBackMesh(const int k, int dom, const R* x, const R* normal) const {
    int theDomain = (domain == -1)? dom : domain;
    return evalOnBackMesh(k,theDomain,x,normal);
  }
  R GevalOnBackMesh(const int k, int dom, const R* x, const R t, const R* normal) const {
    int theDomain = (domain == -1)? dom : domain;
    return evalOnBackMesh(k,theDomain,x,t,normal);
  }


  virtual int size() const {assert(0);}

  ExpressionVirtual& operator= (const ExpressionVirtual& L) {
    cu = L.cu;
    op = L.op;
    opt = L.opt;
    ar_normal.init(L.ar_normal);
    return *this;
  }

  double computeNormal(const R* normal) const {
    if(normal == nullptr) return 1.;
    R val = 1;
    for(int i=0;i<ar_normal.size();++i) val *= normal[ar_normal(i)];
    return val;
  }
  void addNormal(int i) {
    int l = ar_normal.size();
    ar_normal.resize(l+1);
    ar_normal(l) = i;
  }



  virtual void whoAmI() const {std::cout << " I am virtual class Expression" << std::endl;}
  ~ExpressionVirtual(){}

};

template<typename M>
class ExpressionFunFEM : public ExpressionVirtual {

protected:
public:
  const FunFEM<M> & fun;
public:
  ExpressionFunFEM(const FunFEM<M> & fh, int cc, int opp)
  : ExpressionVirtual(cc, opp), fun(fh){}
  ExpressionFunFEM(const FunFEM<M> & fh, int cc, int opp, int oppt)
  : ExpressionVirtual(cc, opp, oppt), fun(fh) {}
  ExpressionFunFEM(const FunFEM<M> & fh, int cc, int opp, int oppt, int dom)
  : ExpressionVirtual(cc, opp, oppt,dom), fun(fh) {}


  R operator()(long i) const { return fun(i);}
  virtual int size() const {fun.size();}

  void whoAmI() const {std::cout << " I am class ExpressionFunFEM" << std::endl;}


  R eval(const int k, const R* x, const R* normal) const {
    return fun.eval(k,x,cu,op) * computeNormal(normal);
  }
  R eval(const int k, const R* x, const R t, const R* normal) const {
    return fun.eval(k,x,t,cu,op, opt) * computeNormal(normal);
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal) const {
    return fun.evalOnBackMesh(k, dom,x,cu,op) * computeNormal(normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal) const {
    return fun.evalOnBackMesh(k,dom,x,t,cu,op, opt) * computeNormal(normal);
  }
public:
  const GFESpace<M>& getSpace() const {return *fun.Vh;}

  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }


  ExpressionFunFEM(const ExpressionFunFEM& L) : ExpressionVirtual(L.cu, L.op, L.opt,L.domain), fun(L.fun)  {
    ar_normal.init(L.ar_normal);
  }

  ExpressionFunFEM& operator= (const ExpressionFunFEM& L) {
    cu = L.cu;
    op = L.op;
    opt = L.opt;
    ar_normal.init(L.ar_normal);
    fun(L.fun);
    domain=L.domain;
    return *this;
  }

  ExpressionFunFEM operator*(const Normal& n) {
    ExpressionFunFEM ff(*this);
    ff.addNormal(cu);
    return ff;
  }


  friend ExpressionFunFEM dx<M>(const ExpressionFunFEM<M>& u);
  friend ExpressionFunFEM dy<M>(const ExpressionFunFEM<M>& u);
  friend ExpressionFunFEM dz<M>(const ExpressionFunFEM<M>& u);
  friend ExpressionFunFEM dt<M>(const ExpressionFunFEM<M>& u);

};


class ExpressionMultConst : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
  const double c;
  const bool nx, ny, nz;
  const R2 p;
public:
  ExpressionMultConst(const ExpressionVirtual & fh1, const double& cc)
  : fun1(fh1) , c(cc), nx(false), ny(false), nz(false), p(R2(1,1)){
  }
  ExpressionMultConst(const ExpressionVirtual & fh1, const Normal_Component_X& nnx)
  : fun1(fh1) , c(1.), nx(true), ny(false), nz(false), p(R2(1,1)){
  }
  ExpressionMultConst(const ExpressionVirtual & fh1, const Normal_Component_Y& nny)
  : fun1(fh1) , c(1.), nx(false), ny(true), nz(false), p(R2(1,1)){
  }
  ExpressionMultConst(const ExpressionVirtual & fh1, const Normal_Component_Z& nnz)
  : fun1(fh1) , c(1.), nx(false), ny(true), nz(true), p(R2(1,1)){
  }
  // ExpressionMultConst(const ExpressionVirtual & fh1, const R2& v)
  // : fun1(fh1) , c(1.), nx(false), ny(false), nz(false), p(R2(v)){
  // }

  R operator()(long i) const {return c*fun1(i);}

  R eval(const int k, const R* x, const R* normal)const  {
    double compN = ((nx)?normal[0]:1)*((ny)?normal[1]:1)*((nz)?normal[2]:1);
    return fun1.eval(k,x, normal) * c * compN * p[0];
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    double compN = ((nx)?normal[0]:1)*((ny)?normal[1]:1)*((nz)?normal[2]:1);
    return fun1.eval(k,x,t, normal) * c * compN * p[0];
  }

  R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
    double compN = ((nx)?normal[0]:1)*((ny)?normal[1]:1)*((nz)?normal[2]:1);
    return fun1.evalOnBackMesh(k,dom,x,normal) * c * compN * p[dom==1];
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    double compN = ((nx)?normal[0]:1)*((ny)?normal[1]:1)*((nz)?normal[2]:1);
    return fun1.evalOnBackMesh(k,dom,x,t,normal) * c * compN * p[dom==1];
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionMultConst(){}
};
ExpressionMultConst operator*(const ExpressionVirtual& f1, double cc);
ExpressionMultConst operator*(double cc, const ExpressionVirtual& f1);
ExpressionMultConst operator*(const ExpressionVirtual& f1, const Normal_Component_X& cc);
ExpressionMultConst operator*(const ExpressionVirtual& f1, const Normal_Component_Y& cc);
ExpressionMultConst operator*(const ParameterCutFEM& v, const ExpressionVirtual& f1);


class ExpressionAbs : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
public:
  ExpressionAbs(const ExpressionVirtual & fh1)
  : fun1(fh1) {
  }

  R operator()(long i) const {return fabs(fun1(i));}

  R eval(const int k, const R* x, const R* normal)const  {
    return fabs(fun1.eval(k,x, normal));
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    return fabs(fun1.eval(k,x,t, normal));
  }

  R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
    return fabs(fun1.evalOnBackMesh(k,dom,x,normal));
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    return fabs(fun1.evalOnBackMesh(k,dom,x,t,normal));
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionAbs(){}
};
ExpressionAbs fabs(const ExpressionVirtual& f1);

 class ExpressionSqrt : public ExpressionVirtual {
   const ExpressionVirtual & fun1;
 public:
   ExpressionSqrt(const ExpressionVirtual & fh1)
   : fun1(fh1){
   }

   R operator()(long i) const {return sqrt(fun1(i));}

   R eval(const int k, const R* x, const R* normal)const  {
     return sqrt(fun1.eval(k,x, normal));
   }
   R eval(const int k, const R* x, const R t, const R* normal)const  {
     return sqrt(fun1.eval(k,x,t,normal));
   }

   R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
     return sqrt(fun1.evalOnBackMesh(k,dom,x,normal));
   }
   R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
     return sqrt(fun1.evalOnBackMesh(k,dom,x,t,normal));
   }
   int idxElementFromBackMesh(int kb, int dd=0) const {
     return fun1.idxElementFromBackMesh(kb, dd);
   }
   ~ExpressionSqrt(){}
 };
  ExpressionSqrt sqrt(const ExpressionVirtual& f1);

class ExpressionProduct : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
  const ExpressionVirtual & fun2;
public:
  ExpressionProduct(const ExpressionVirtual & fh1, const ExpressionVirtual & fh2)
  : fun1(fh1), fun2(fh2) {}

  R operator()(long i) const {return fun1(i)*fun2(i);}

  R eval(const int k, const R* x, const R* normal)const  {
    return fun1.eval(k,x,normal) * fun2.eval(k,x,normal);
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    return fun1.eval(k,x,t,normal) * fun2.eval(k,x,t,normal);
  }


  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,normal) * fun2.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,t,normal) * fun2.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionProduct(){}
};
ExpressionProduct operator*(const ExpressionVirtual& f1, const ExpressionVirtual& f2);

class ExpressionDivision : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
  const ExpressionVirtual & fun2;
public:
  ExpressionDivision(const ExpressionVirtual & fh1, const ExpressionVirtual & fh2)
  : fun1(fh1), fun2(fh2) {}

  R operator()(long i) const {return fun1(i)/fun2(i);}

  R eval(const int k, const R* x, const R* normal)const  {
    double v = fun2.eval(k,x,normal);
    assert(fabs(v) > 1e-15);
    return fun1.eval(k,x,normal) / v;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    double v = fun2.eval(k,x,t,normal);
    assert(fabs(v) > 1e-15);
    return fun1.eval(k,x,t,normal) / v;
  }


  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    double v = fun2.evalOnBackMesh(k,dom,x,normal);
    assert(fabs(v) > 1e-15);
    return fun1.evalOnBackMesh(k,dom,x,normal) / v;
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    double v = fun2.evalOnBackMesh(k,dom,x,t,normal);
    assert(fabs(v) > 1e-15);
    return fun1.evalOnBackMesh(k,dom,x,t,normal) / v;
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDivision(){}
};
ExpressionDivision operator/(const ExpressionVirtual& f1, const ExpressionVirtual& f2);


class ExpressionSum : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
  const ExpressionVirtual & fun2;
public:
  ExpressionSum(const ExpressionVirtual & fh1, const ExpressionVirtual & fh2)
  : fun1(fh1), fun2(fh2) {
  }
  R operator()(long i) const {return fun1(i) + fun2(i);}

  R eval(const int k, const R* x, const R* normal)const  {
    return fun1.eval(k,x,normal) + fun2.eval(k,x,normal);
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    return fun1.eval(k,x,t,normal) + fun2.eval(k,x,t,normal);
  }


  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,normal) + fun2.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,t,normal) + fun2.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionSum(){}
};
ExpressionSum operator+(const ExpressionVirtual& f1, const ExpressionVirtual& f2);

class ExpressionDif : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
  const ExpressionVirtual & fun2;
public:
  ExpressionDif(const ExpressionVirtual & fh1, const ExpressionVirtual & fh2)
  : fun1(fh1), fun2(fh2) {}

  R operator()(long i) const {return fun1(i)-fun2(i);}

  R eval(const int k, const R* x, const R* normal)const  {
    return fun1.eval(k,x,normal) - fun2.eval(k,x,normal);
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    return fun1.eval(k,x,t,normal) - fun2.eval(k,x,t,normal);
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,normal) - fun2.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    return fun1.evalOnBackMesh(k,dom,x,t,normal) - fun2.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDif(){}
};
ExpressionDif operator-(const ExpressionVirtual& f1, const ExpressionVirtual& f2);


class ExpressionNormal2 : public ExpressionVirtual {
  typedef Mesh2 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> uxnx, uyny;
  double c0 = 1;

public:
  ExpressionNormal2(const FunFEM<M> & fh1, const Normal n)
  : fun(fh1) , uxnx(fh1,0,op_id,0,0), uyny(fh1,1,op_id,0,0)
  {
    assert(fh1.Vh->N !=1);
    uxnx.addNormal(0); uyny.addNormal(1);
  }
  ExpressionNormal2(const FunFEM<M> & fh1, const Tangent t)
  : fun(fh1) , uxnx(fh1,0,op_id,0,0), uyny(fh1,1,op_id,0,0)
  {
    assert(fh1.Vh->N !=1);
    uxnx.addNormal(1); uyny.addNormal(0);
    c0 = -1;
  }
  ExpressionNormal2(const FunFEM<M> & fh1, const Conormal n)
  : fun(fh1) , uxnx(fh1,0,op_id,0,0), uyny(fh1,1,op_id,0,0)
  {
    assert(fh1.Vh->N !=1);
    uxnx.addNormal(0); uyny.addNormal(1);
  }

  R operator()(long i) const {assert(0);return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return c0*uxnx.evalOnBackMesh(k,dom,x,normal) + uyny.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    assert(normal);
    return c0*uxnx.evalOnBackMesh(k,dom,x,t,normal) + uyny.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionNormal2(){}
};
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Normal& n);
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Tangent& n);
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Conormal& n);


class ExpressionNormal3 : public ExpressionVirtual {
  typedef Mesh3 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> uxnx, uyny, uznz;

public:
  ExpressionNormal3(const FunFEM<M> & fh1)
  : fun(fh1) , uxnx(fh1,0,op_id,0,0), uyny(fh1,1,op_id,0,0), uznz(fh1,2,op_id,0,0)
  {
    assert(fh1.Vh->N !=1);
    uxnx.addNormal(0); uyny.addNormal(1); uznz.addNormal(2);
  }

  R operator()(long i) const {assert(0); return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return uxnx.evalOnBackMesh(k,dom,x,normal) + uyny.evalOnBackMesh(k,dom,x,normal)
    + uznz.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    assert(normal);
    return uxnx.evalOnBackMesh(k,dom,x,t,normal) + uyny.evalOnBackMesh(k,dom,x,t,normal)
    + uznz.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionNormal3(){}
};
ExpressionNormal3 operator*(const FunFEM<Mesh3>& f1, const Normal& n);



// divS for 2d
class ExpressionDSx2 : public ExpressionVirtual {
  typedef Mesh2 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> dxu1, dxu1nxnx, dyu1nxny;

public:
  ExpressionDSx2(const FunFEM<M> & fh1)
  : fun(fh1) , dxu1(fh1,0,op_dx,0,0), dxu1nxnx(fh1,0,op_dx,0,0), dyu1nxny(fh1,0,op_dy,0,0)
  {
    dxu1nxnx.addNormal(0); dxu1nxnx.addNormal(0);
    dyu1nxny.addNormal(0); dyu1nxny.addNormal(1);
  }

  R operator()(long i) const {assert(0); return  0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dxu1.evalOnBackMesh(k,dom,x,normal) - dxu1nxnx.evalOnBackMesh(k,dom,x,normal) - dyu1nxny.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    assert(normal);
    return dxu1.evalOnBackMesh(k,dom,x,t,normal) - dxu1nxnx.evalOnBackMesh(k,dom,x,t,normal) - dyu1nxny.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDSx2(){}
};
ExpressionDSx2 dxS(const FunFEM<Mesh2>& f1);

class ExpressionDSy2 : public ExpressionVirtual {
  typedef Mesh2 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> dxu2, dxu2nxny, dyu2nyny;

public:
  ExpressionDSy2(const FunFEM<M> & fh1)
  : fun(fh1) , dxu2(fh1,1,op_dy,0,0), dxu2nxny(fh1,1,op_dx,0,0), dyu2nyny(fh1,1,op_dy,0,0)
  {
    dxu2nxny.addNormal(0); dxu2nxny.addNormal(1);
    dyu2nyny.addNormal(1); dyu2nyny.addNormal(1);
  }

  R operator()(long i) const {assert(0); return  0; };

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dxu2.evalOnBackMesh(k,dom,x,normal) - dxu2nxny.evalOnBackMesh(k,dom,x,normal) - dyu2nyny.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x,const R t, const R* normal)const  {
    assert(normal);
    return dxu2.evalOnBackMesh(k,dom,x,t,normal) - dxu2nxny.evalOnBackMesh(k,dom,x,t,normal) - dyu2nyny.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDSy2(){}
};
ExpressionDSy2 dyS(const FunFEM<Mesh2>& f1);

class ExpressionDivS2 : public ExpressionVirtual {
  typedef Mesh2 M;
  const FunFEM<M>& fun;
  const ExpressionDSx2 dx;
  const ExpressionDSy2 dy;
  // const ExpressionDSz<M> dz;

public:
  ExpressionDivS2(const FunFEM<M> & fh1) : fun(fh1), dx(dxS(fh1)),dy(dyS(fh1)){}

    R operator()(long i) const {assert(0); return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dx.evalOnBackMesh(k,dom,x,normal) + dy.evalOnBackMesh(k,dom,x,normal) ;
  }
  R evalOnBackMesh(const int k, const int dom, const R* x,const R t, const R* normal)const  {
    assert(normal);
    return dx.evalOnBackMesh(k,dom,x,t,normal) + dy.evalOnBackMesh(k,dom,x,t,normal) ;
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDivS2(){}
};
ExpressionDivS2 divS(const FunFEM<Mesh2>& f1);


// divs for 3D
class ExpressionDSx3 : public ExpressionVirtual {
  typedef Mesh3 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> dxu1, dxu1nxnx, dyu1nxny , dzu1nxnz;

public:
  ExpressionDSx3(const FunFEM<M> & fh1)
  : fun(fh1) , dxu1(fh1,0,op_dx,0,0), dxu1nxnx(fh1,0,op_dx,0,0), dyu1nxny(fh1,0,op_dy,0,0) , dzu1nxnz(fh1,0,op_dz,0,0)
  {
    dxu1nxnx.addNormal(0); dxu1nxnx.addNormal(0);
    dyu1nxny.addNormal(0); dyu1nxny.addNormal(1);
    dzu1nxnz.addNormal(0); dzu1nxnz.addNormal(2);
  }

  R operator()(long i) const {assert(0); return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dxu1.evalOnBackMesh(k,dom,x,normal) - dxu1nxnx.evalOnBackMesh(k,dom,x,normal) - dyu1nxny.evalOnBackMesh(k,dom,x,normal) - dzu1nxnz.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    assert(normal);
    return dxu1.evalOnBackMesh(k,dom,x,t,normal) - dxu1nxnx.evalOnBackMesh(k,dom,x,t,normal) - dyu1nxny.evalOnBackMesh(k,dom,x,t,normal)- dzu1nxnz.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDSx3(){}
};
ExpressionDSx3 dxS(const FunFEM<Mesh3>& f1);

class ExpressionDSy3 : public ExpressionVirtual {
  typedef Mesh3 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> dxu2, dxu2nxny, dyu2nyny, dzu2nynz;

public:
  ExpressionDSy3(const FunFEM<M> & fh1)
  : fun(fh1) , dxu2(fh1,1,op_dy,0,0), dxu2nxny(fh1,1,op_dx,0,0), dyu2nyny(fh1,1,op_dy,0,0), dzu2nynz(fh1,1,op_dz,0,0)
  {
    dxu2nxny.addNormal(0); dxu2nxny.addNormal(1);
    dyu2nyny.addNormal(1); dyu2nyny.addNormal(1);
    dzu2nynz.addNormal(1); dzu2nynz.addNormal(2);
  }

  R operator()(long i) const {assert(0);return  0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dxu2.evalOnBackMesh(k,dom,x,normal) - dxu2nxny.evalOnBackMesh(k,dom,x,normal) - dyu2nyny.evalOnBackMesh(k,dom,x,normal)- dzu2nynz.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x,const R t, const R* normal)const  {
    assert(normal);
    return dxu2.evalOnBackMesh(k,dom,x,t,normal) - dxu2nxny.evalOnBackMesh(k,dom,x,t,normal) - dyu2nyny.evalOnBackMesh(k,dom,x,t,normal)- dzu2nynz.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDSy3(){}
};
ExpressionDSy3 dyS(const FunFEM<Mesh3>& f1);

class ExpressionDSz3 : public ExpressionVirtual {
  typedef Mesh3 M;
  const FunFEM<M>& fun;
  ExpressionFunFEM<M> dxu3, dxu3nxnz, dyu3nynz, dzu3nznz;

public:
  ExpressionDSz3(const FunFEM<M> & fh1)
  : fun(fh1) , dxu3(fh1,2,op_dz,0,0), dxu3nxnz(fh1,2,op_dx,0,0), dyu3nynz(fh1,2,op_dy,0,0), dzu3nznz(fh1,2,op_dz,0,0)
  {
    dxu3nxnz.addNormal(0); dxu3nxnz.addNormal(2);
    dyu3nynz.addNormal(1); dyu3nynz.addNormal(2);
    dzu3nznz.addNormal(2); dzu3nznz.addNormal(2);
  }

  R operator()(long i) const {assert(0); return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dxu3.evalOnBackMesh(k,dom,x,normal) - dxu3nxnz.evalOnBackMesh(k,dom,x,normal) - dyu3nynz.evalOnBackMesh(k,dom,x,normal)- dzu3nznz.evalOnBackMesh(k,dom,x,normal);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x,const R t, const R* normal)const  {
    assert(normal);
    return dxu3.evalOnBackMesh(k,dom,x,t,normal) - dxu3nxnz.evalOnBackMesh(k,dom,x,t,normal) - dyu3nynz.evalOnBackMesh(k,dom,x,t,normal)- dzu3nznz.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionDSz3(){}
};
ExpressionDSz3 dzS(const FunFEM<Mesh3>& f1);

class ExpressionDivS3 : public ExpressionVirtual {
  typedef Mesh3 M;
  const FunFEM<M>& fun;
  const ExpressionDSx3 dx;
  const ExpressionDSy3 dy;
  const ExpressionDSz3 dz;

public:
  ExpressionDivS3(const FunFEM<M> & fh1) : fun(fh1), dx(dxS(fh1)),dy(dyS(fh1)),dz(dzS(fh1)){}

  R operator()(long i) const {assert(0); return 0;};

  R eval(const int k, const R* x, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
    assert(0);
    return 0;
  }

  R evalOnBackMesh(const int k, const int dom, const R* x, const R* normal)const  {
    assert(normal);
    return dx.evalOnBackMesh(k,dom,x,normal) + dy.evalOnBackMesh(k,dom,x,normal) + dz.evalOnBackMesh(k,dom,x,normal) ;
  }
  R evalOnBackMesh(const int k, const int dom, const R* x,const R t, const R* normal)const  {
    assert(normal);
    return dx.evalOnBackMesh(k,dom,x,t,normal) + dy.evalOnBackMesh(k,dom,x,t,normal) + dz.evalOnBackMesh(k,dom,x,t,normal);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun.idxElementFromBackMesh(kb, dd);
  }

  ~ExpressionDivS3(){}
};
ExpressionDivS3 divS(const FunFEM<Mesh3>& f1);


class ExpressionAverage {
  public:
  const ExpressionVirtual & fun1;
  const R k1, k2;

  ExpressionAverage(const ExpressionVirtual & fh1, double kk1, double kk2)
  : fun1(fh1), k1(kk1), k2(kk2){
  }

  ~ExpressionAverage(){}
};
ExpressionAverage average(const ExpressionVirtual & fh1, const double kk1=0.5, const double kk2=0.5);
ExpressionAverage jump(const ExpressionVirtual & fh1, const double kk1=1, const double kk2=-1);
ExpressionAverage operator* (double c, const ExpressionAverage& fh);
ExpressionAverage operator* ( const ExpressionAverage& fh, double c);


template<typename M>
ExpressionFunFEM<M> dx(const ExpressionFunFEM<M>& u) {
  return ExpressionFunFEM<M>(u.fun, u.cu, op_dx, u.opt, u.domain);
}
template<typename M>
ExpressionFunFEM<M> dy(const ExpressionFunFEM<M>& u) {
  return ExpressionFunFEM<M>(u.fun, u.cu, op_dy, u.opt, u.domain);
}
template<typename M>
ExpressionFunFEM<M> dz(const ExpressionFunFEM<M>& u) {
  return ExpressionFunFEM<M>(u.fun, u.cu, op_dz, u.opt, u.domain);
}
template<typename M>
ExpressionFunFEM<M> dt(const ExpressionFunFEM<M>& u) {
  return ExpressionFunFEM<M>(u.fun, u.cu, u.op, op_dx, u.domain);
}



class ExpressionBurgerFlux : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
public:
  ExpressionBurgerFlux(const ExpressionVirtual & fh1)
  : fun1(fh1) {
  }

  R operator()(long i) const {return fabs(fun1(i));}

  R eval(const int k, const R* x, const R* normal)const  {
    double val = fun1.eval(k,x, normal);
    return 0.5*val*val;
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    double val = fun1.eval(k,x,t, normal);
    return 0.5*val*val;
  }

  R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
    double val = fun1.evalOnBackMesh(k,dom,x,normal);
    return 0.5*val*val;
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    double val = fun1.evalOnBackMesh(k,dom,x,t,normal);
    return 0.5*val*val;
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionBurgerFlux(){}
};
class ExpressionNormalBurgerFlux : public ExpressionVirtual {
  const ExpressionVirtual & fun1;
public:
  ExpressionNormalBurgerFlux(const ExpressionVirtual & fh1)
  : fun1(fh1) {
  }

  R operator()(long i) const {return fabs(fun1(i));}

  R eval(const int k, const R* x, const R* normal)const  {
    assert(normal);
    double val = fun1.eval(k,x, normal);
    return 0.5*val*val*(normal[0] + normal[1]);
  }
  R eval(const int k, const R* x, const R t, const R* normal)const  {
    assert(normal);
    double val = fun1.eval(k,x,t, normal);
    return 0.5*val*val*(normal[0] + normal[1]);
  }

  R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
    assert(normal);
    double val = fun1.evalOnBackMesh(k,dom,x,normal);
    return 0.5*val*val*(normal[0] + normal[1]);
  }
  R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
    assert(normal);
    double val = fun1.evalOnBackMesh(k,dom,x,t,normal);
    return 0.5*val*val*(normal[0] + normal[1]);
  }
  int idxElementFromBackMesh(int kb, int dd=0) const {
    return fun1.idxElementFromBackMesh(kb, dd);
  }
  ~ExpressionNormalBurgerFlux(){}
};
ExpressionBurgerFlux burgerFlux(const ExpressionVirtual & f1);
ExpressionNormalBurgerFlux burgerFlux(const ExpressionVirtual & f1, const Normal& n ) ;


// class ExpressionNormalBurgerFlux : public ExpressionVirtual {
//   const ExpressionVirtual & fun1;
// public:
//   ExpressionNormalBurgerFlux(const ExpressionVirtual & fh1)
//   : fun1(fh1) {
//   }
//
//   R operator()(long i) const {return fabs(fun1(i));}
//
//   R eval(const int k, const R* x, const R* normal)const  {
//     double val = fun1.eval(k,x, normal)
//     return 0.5*val*val;
//   }
//   R eval(const int k, const R* x, const R t, const R* normal)const  {
//     double val = fun1.eval(k,x,t, normal);
//     return 0.5*val*val;
//   }
//
//   R evalOnBackMesh(const int k, const int dom,const R* x, const R* normal)const  {
//     double val = fun1.evalOnBackMesh(k,dom,x,normal);
//     return 0.5*val*val;
//   }
//   R evalOnBackMesh(const int k, const int dom, const R* x, const R t, const R* normal)const  {
//     double val = fun1.evalOnBackMesh(k,dom,x,t,normal);
//     return 0.5*val*val;
//   }
//   int idxElementFromBackMesh(int kb, int dd=0) const {
//     return fun1.idxElementFromBackMesh(kb, dd);
//   }
//   ~ExpressionNormalBurgerFlux(){}
// };


#include "expression.tpp"


typedef FunFEM<Mesh2> Fun2_h;
typedef ExpressionFunFEM<Mesh2> Expression2;
typedef FunFEM<Mesh3> Fun3_h;
typedef ExpressionFunFEM<Mesh3> Expression3;
#endif
