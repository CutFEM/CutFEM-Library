#ifndef _ITEMVF_HPP
#define _ITEMVF_HPP
#include "testFunction.hpp"
#include <list>

template<int N = 2>
struct ItemVF {
  typedef typename typeMesh<N>::Mesh Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename typeRd<N>::Rd Rd;
  double c;
  int cu,du,cv,dv;
  KN<int> ar_nu, ar_nv;
  int domu, domv;
  vector<string> coefu, coefv ;
  int dtu, dtv;
  const ExpressionVirtual* expru=nullptr;
  const ExpressionVirtual* exprv=nullptr;

  FESpace const * fespaceU = nullptr;
  FESpace const * fespaceV = nullptr;

  void(*pfunU)(RNMK_&,int,int) = f_id;
  void(*pfunV)(RNMK_&,int,int) = f_id;

  CutFEM_ParameterList parameterList;


  ItemVF()
  : c(0.), cu(-1),du(-1),cv(-1),dv(-1),domu(-1),domv(-1),dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l)
  : c(cc), cu(i),du(j),cv(k),dv(l),domu(-1),domv(-1),dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm)
  : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),domu(-1),domv(-1),dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm, int dou, int dov)
  : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),domu(dou),domv(dov),dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm, int dou, int dov,const vector<string>& weu, const vector<string>& wev, int tu, int tv)
  : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),domu(dou),domv(dov),dtu(tu),dtv(tv){
    for(int i=0;i<weu.size();++i) coefu.push_back(weu[i]);
    for(int i=0;i<wev.size();++i) coefv.push_back(wev[i]);
  }
  ItemVF(const ItemVF & U)
  :ItemVF(U.c,U.cu,U.du,U.cv,U.dv,U.ar_nu,U.ar_nv,U.domu,U.domv){
    for(int i=0;i<U.coefu.size();++i) coefu.push_back(U.coefu[i]);
    for(int i=0;i<U.coefv.size();++i) coefv.push_back(U.coefv[i]);
    dtu = U.dtu;
    dtv = U.dtv;
    expru = U.expru;
    exprv = U.exprv;
    fespaceU = U.fespaceU;
    fespaceV = U.fespaceV;
    pfunU = U.pfunU;
    pfunV = U.pfunV;
  }

  ItemVF(const ItemTestFunction<N>& U, const ItemTestFunction<N>& V)
    :ItemVF(U.c*V.c,U.cu,U.du,V.cu,V.du,U.ar_nu,V.ar_nu,
            U.dom,V.dom, U.coefu, V.coefu,  U.dtu, V.dtu)
  {
    expru = U.expru;
    exprv = V.expru;

    fespaceU = U.fespace;
    fespaceV = V.fespace;

    pfunU = U.pfun;
    pfunV = V.pfun;
  }

  bool on(int d) const {
    return ((domu == domv) && (domu == -1 || domu == d));
  }

  bool same() const {
    return (fespaceU == fespaceV) && (pfunU == pfunV); }

  bool operator==(const ItemVF& F){
    if(cu == F.cu && cv == F.cv && du == F.du && dv == F.dv
        && F.domu == domu && domv == F.domv && dtu == F.dtu && dtv == F.dtv
        && fespaceU == F.fespaceU && fespaceV == F.fespaceV
        && expru == F.expru && exprv == F.exprv) {

      if(ar_nu.size() == F.ar_nu.size() && ar_nv.size() == F.ar_nv.size()) {
        for(int i=0;i<ar_nu.size();++i) {
          if(ar_nu(i) != F.ar_nu(i)) return false;
        }
        for(int i=0;i<ar_nv.size();++i) {
          if(ar_nv(i) != F.ar_nv(i)) return false;
        }
        return true;
      }
    }
    return false;
  }


  double fxU(int k, Rd mip, const R* normal = nullptr) const {
    return ((expru)? expru->eval(k, mip, normal) : 1);
  }
  double fxV(int k, Rd mip, const R* normal = nullptr) const {
    return ((exprv)? exprv->eval(k, mip, normal) : 1);
  }
  double fxu(int k, Rd mip) const {
    return ((expru)? expru->eval(k, mip) : 1)*((exprv)? exprv->eval(k, mip) : 1);
  }
  double fxu(int k, Rd mip, double t) const {
    return ((expru)? expru->eval(k, mip, t) : 1)*((exprv)? exprv->eval(k, mip,t) : 1);
  }

  double fxu_backMesh(int k, int dom, Rd mip, const R* normal = nullptr) const {
    return ((expru)? expru->GevalOnBackMesh(k, dom, mip, normal) : 1)
          *((exprv)? exprv->GevalOnBackMesh(k, dom, mip, normal) : 1);
  }
  double fxu_backMesh(int k, int dom, Rd mip, double t, const R* normal = nullptr) const {
    return ((expru)? expru->GevalOnBackMesh(k, dom, mip, t, normal) : 1)
           *((exprv)? exprv->GevalOnBackMesh(k, dom, mip, t, normal) : 1);
  }

  double fx_backMesh_U(int k, int dom, Rd mip, const R* normal = nullptr) const {
    return ((expru)? expru->GevalOnBackMesh(k, dom, mip, normal) : 1);
  }
  double fx_backMesh_V(int k, int dom, Rd mip, const R* normal = nullptr) const {
    return ((exprv)? exprv->GevalOnBackMesh(k, dom, mip, normal) : 1);
  }
public:
  R getCoefU(const R* normal) const {
    R val = 1;
    for(int i=0;i<ar_nu.size();++i) val *= normal[ar_nu(i)];
    return val;
  }
  R getCoefV(const R* normal) const {
    R val = 1;
    for(int i=0;i<ar_nv.size();++i) val *= normal[ar_nv(i)];
    return val;
  }

public:
  R getCoef(const R* normal) const {
    return getCoefU(normal) * getCoefV(normal);
  }
  R computeNormal(const R* normal) const {
    return getCoef(normal);
  }

  R computeCoef(double h, double meas, double measK, int domain) const {
    R val = 1;
    for(int l=0;l<2;++l) {
      const vector<string>& listCoef = (l==0)?coefu : coefv;
      for(int i=0;i<listCoef.size();++i) {
        string coef = listCoef[i];

        if(parameterList.find(coef)) {
          CutFEM_Parameter& p(*parameterList.listParameter[coef]);
          val *= p(domain, h, meas, measK);
        }
      }
    }
    return val;
  }
  R computeCoefInterface( double h, double meas, double measK) const {
    R val = 1;
    for(int l=0;l<2;++l) {
      const vector<string>& listCoef = (l==0)?coefu : coefv;
      int domCoef = (l==0)?domu : domv;
      if(domCoef == -1) domCoef = 0;
      assert(domCoef == 0 || domCoef == 1);
      for(int i=0;i<listCoef.size();++i) {
        string coef = listCoef[i];

        if(this->parameterList.find(coef)) {
          CutFEM_Parameter& p(*this->parameterList.listParameter[coef]);
          val *= p(domCoef, h, meas, measK);
        }
      }
    }
    return val;
  }
  R computeCoef(int domain, const typename Mesh::Partition& cutK) const {
    R val = 1;
    double h = cutK.getEdgeLength();
    double meas = cutK.mesure(domain);
    double measK = cutK.T.mesure();
    for(int l=0;l<2;++l) {
      const vector<string>& listCoef = (l==0)?coefu : coefv;
      int domCoef = domain;
      for(int i=0;i<listCoef.size();++i) {
        string coef = listCoef[i];

        if(this->parameterList.find(coef)) {
          CutFEM_Parameter& p(*this->parameterList.listParameter[coef]);
          val *= p(domCoef, h, meas, measK);
        }
      }
    }
    return val;
  }

  void applyFunNL(RNMK_& bfu, RNMK_& bfv) const {
    pfunU(bfu, cu, du);
    pfunV(bfv, cv, dv);
  }


  friend std::ostream& operator <<(std::ostream& f, const ItemVF & u )
  {
    string n[3] = {"nx", "ny", "nz"};
    f << " FESpaces => " << u.fespaceU << " and " << u.fespaceV << "\t";
    f << u.c << "\t" << whichOperator( u.dtu) << whichOperator( u.du,u.cu);
    for(int i=0;i<u.ar_nu.size();++i) f << " * " << n[u.ar_nu(i)];
    for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
    f << " * " << whichOperator( u.dtv) << whichOperatorV( u.dv,u.cv);
    for(int i=0;i<u.ar_nv.size();++i) f << " * " << n[u.ar_nv(i)];
    for(int i=0;i<u.coefv.size();++i) f << " * " << u.coefv[i];
    if(u.domu == u.domv && u.domu != -1) f << "\t in Omega_" << u.domu+1;

    f << std::endl;
    return f;
  }
};

template<int N=2>
ItemVF<N>& operator*=(R cc, ItemVF<N>& F) {
  F.c *= cc;
  return F;
}
template<int N=2>
ItemVF<N>& operator*=(ItemVF<N>& F, R cc) {
  F.c *= cc;
  return F;
}
template<int N=2>
ItemVF<N> operator-(const ItemVF<N>& X) {
  ItemVF<N> F(X);
  F.c *= -1;
  return F;
}

template<int N = 2>
class  ListItemVF {

  public :
  KN<ItemVF<N>> VF;
  ListItemVF(int l) : VF(l) {};
  const ItemVF<N>& operator()(int i) const {return VF(i);}
  const ItemVF<N>& operator[](int i) const {return VF[i];}
  ItemVF<N>& operator()(int i){return VF(i);}
  ItemVF<N>& operator[](int i){return VF[i];}
  int size() const {return VF.size();}

  ListItemVF& operator* (const double cc) {
    for(int i=0;i<VF.size();++i) (*this)(i).c *= cc;
    return *this;
  }


  ListItemVF& operator+ (const ListItemVF& L) {
    int n0 = VF.size();
    int n = L.size() + this->size();
    this->VF.resize(n);
    for(int i=n0,j=0;i<n;++i,++j) (*this)(i) = L(j);
    return *this;
  }

  ListItemVF& operator- (const ListItemVF& L) {
    int n0 = VF.size();
    int n = L.size() + this->size();
    this->VF.resize(n);
    for(int i=n0,j=0;i<n;++i,++j) (*this)(i) = -L(j);
    return *this;
  }

  ListItemVF& operator- () {
    for(int i=0;i<this->size();++i) (*this)(i) = -(*this)(i);
    return *this;
  }
  ListItemVF& operator+ () {
    return *this;
  }

  void reduce() {
    // get size new list
    int l = VF.size();
    KN<int> s(VF.size(),-1);
    KN<int> s2k(VF.size(),-1);

    for (int i=0; i<VF.size();++i) {
      if(VF(i).c == 0){l-= 1; continue;}
      if (s(i) != -1) continue;
      for(int j=i+1; j<VF.size();++j) {
        if(VF(i) == VF(j)) {
          s(j) = i;
          l -= 1;
        }
      }
    }
    if(l == VF.size()) return;
    KN_<ItemVF<N>> v(VF);
    ItemVF<N> * u = new ItemVF<N>[l];
    VF.set(u,l);
    int k=0;

    for(int i=0;i<v.size();++i ) {
      if(v(i).c == 0){continue;}
      if(s(i) == -1) {
        s2k[i] = k;
        u[k++] = v(i);
      }
      else {
        u[s2k(s(i))].c += v(i).c;

      }
    }
  }


  friend std::ostream& operator <<(std::ostream& f, const ListItemVF & u )
  {
    f << u.VF;
    return f;
  }
};


template <int d>
ListItemVF<d> jump(const ListItemVF<d>& L) {
  int n0 = L.VF.size();
  int n = 2*n0;

  ListItemVF<d> item(n);
  for(int i=0;i<n0;++i) {
    item(i) = L(i);
    item(i).domu = 0; item(i).domv = 0;
  }
  for(int i=n0,j=0;i<n;++i,++j) {
    item(i) = L(j);
    item(i).c *= -1;
    item(i).domu = 1; item(i).domv = 1;
  }
  return item;
}

template <int d>
ListItemVF<d> average(const ListItemVF<d>& L, int v1, int v2) {
  int n0 = L.VF.size();
  int n = 2*n0;

  ListItemVF<d> item(n);
  for(int i=0;i<n0;++i) {
    item(i) = L(i);
    item(i).c *= v1;
    item(i).domu = 0; item(i).domv = 0;
  }
  for(int i=n0,j=0;i<n;++i,++j) {
    item(i) = L(j);
    item(i).c *= v2;
    item(i).domu = 1; item(i).domv = 1;
  }
  return item;
}


// only for vectors
template <int d>
ListItemVF<d> operator,(const TestFunction<d>& uu, const TestFunction<d>& vv) {
  assert(uu.A.N() == vv.A.N());
  assert(uu.A.M() == vv.A.M());// && A.M()==1);
  int l = 0;
  for(int i=0;i<uu.A.N();++i) {
    for(int j=0;j<uu.A.M();++j) {
      l += uu.A(i,j)->size() * vv.A(i,j)->size() ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  for(int i=0;i<uu.A.N();++i){
    for(int j=0;j<uu.A.M();++j){
      for(int ui=0;ui<uu.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& u(uu.A(i,j)->getItem(ui));
        for(int vi=0;vi<vv.A(i,j)->size();++vi) {
          const ItemTestFunction<d>& v(vv.A(i,j)->getItem(vi));
          item(k) = ItemVF<d>( u, v);
          k++;
        }
      }
    }
  }
  item.reduce();
  return item;
}

// only for vectors
template <int d>
ListItemVF<d> operator,(const R c, const TestFunction<d>& F) {
  int l = 0;
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      l += F.A(i,j)->size() ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  for(int i=0;i<F.A.N();++i){
    for(int j=0;j<F.A.M();++j){
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& v(F.A(i,j)->getItem(ui));
        item(k) = ItemVF<d>( v.c*c,v.cu,0,v.cu,v.du,0,v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,-1,v.dtu);
        item(k).expru = v.expru;
        item(k).fespaceV = v.fespace;
        k++;
      }
    }
  }
  item.reduce();

  return item;
}

// only for vectors
template <int d>
ListItemVF<d> operator,(const Rnm& c, const TestFunction<d>& F) {
  assert(c.N() == F.A.N() &&  c.M() == F.A.M());
  int l = 0;
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      l += F.A(i,j)->size() ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  for(int i=0;i<F.A.N();++i){
    for(int j=0;j<F.A.M();++j){
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& v(F.A(i,j)->getItem(ui));
        item(k) = ItemVF<d>( v.c*c(i,j),v.cu,0,v.cu,v.du,0,v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,-1,v.dtu);
        item(k).exprv = v.expru;
        item(k).fespaceV = v.fespace;

        k++;
      }
    }
  }
  item.reduce();

  return item;
}

// only for vectors
template <int d>
ListItemVF<d> operator,(const Projection& c, const TestFunction<d>& F) {
  int l = 0;
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      l += F.A(i,j)->size() * (1 + (i == j)) ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  for(int i=0;i<F.A.N();++i){
    for(int j=0;j<F.A.M();++j){
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& v(F.A(i,j)->getItem(ui));

        if(i==j) {
          item(k) = ItemVF<d>( v.c,v.cu,0,v.cu,v.du,0,v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,-1,v.dtu);
          item(k).exprv = v.expru;
          item(k).fespaceV = v.fespace;
          k++;
        }
        item(k) = ItemVF<d>( v.c*(-1),v.cu,0,v.cu,v.du,c(i,j),v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,-1,v.dtu);
        item(k).exprv = v.expru;
        item(k).fespaceV = v.fespace;
        k++;

      }
    }
  }
  item.reduce();

  return item;
}

template <int d>
ListItemVF<d> operator,(const ExpressionVirtual& fh, const TestFunction<d>& F) {
  int l = 0;
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      l += F.A(i,j)->size() ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  for(int i=0;i<F.A.N();++i){
    for(int j=0;j<F.A.M();++j){
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& v(F.A(i,j)->getItem(ui));
        item(k) = ItemVF<d>( v.c,0,0,v.cu,v.du,0,v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,0,v.dtu);
        item(k).expru = &fh;
        item(k).exprv = v.expru;
        item(k).fespaceV = v.fespace;

        k++;
      }
    }
  }

  item.reduce();
  return item;
}

template <int d>
ListItemVF<d> operator,(std::list<ExpressionFunFEM<typename typeMesh<d>::Mesh>*> fh, const TestFunction<d>& F) {
  if(F.A.N() != fh.size()){
    std::cout << "size expression \t" << fh.size() << std::endl;
    std::cout << "size test function \t" << F.A.N() << std::endl;
  }
  assert(F.A.N() == fh.size());
  int l = 0;
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      l += F.A(i,j)->size() ;
    }
  }

  ListItemVF<d> item(l);
  int k=0, kloc=0;
  auto it = fh.begin();
  for(int i=0;i<F.A.N();++i, ++it){
    for(int j=0;j<F.A.M();++j){
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<d>& v(F.A(i,j)->getItem(ui));
        item(k) = ItemVF<d>( v.c,0,0,v.cu,v.du,0,v.ar_nu,v.dom,v.dom,vector<string>(), v.coefu,0,v.dtu);
        item(k).expru = *it;
        item(k).exprv = v.expru;
        item(k).fespaceV = v.fespace;

        k++;
      }
    }
  }

  item.reduce();
  return item;
}


template <int d>
ListItemVF<d> innerProduct(double c, const TestFunction<d>& F) {
  return operator,(c, F);
}

template <int d>
ListItemVF<d> innerProduct(const ExpressionVirtual& fh, const TestFunction<d>& F) {
  return operator,(fh, F);
}

template <int d>
ListItemVF<d> innerProduct(std::list<ExpressionFunFEM<typename typeMesh<d>::Mesh>> fh, const TestFunction<d>& F) {
  std::list<ExpressionFunFEM<typename typeMesh<d>::Mesh>*> pointer2fh;
  for(  auto it = fh.begin();it!=fh.end();++it) {
    pointer2fh.push_back(&(*it));
  }
  return operator,(pointer2fh, F);

}

template <int d>
ListItemVF<d> innerProduct(const TestFunction<d>& F1, const TestFunction<d>& F2) {
  return (F1, F2);
}

template <int d>
ListItemVF<d> contractProduct(const TestFunction<d>& F1, const TestFunction<d>& F2) {
  return (F1, F2);
}

template <int d>
ListItemVF<d> contractProduct(const Rnm& F1, const TestFunction<d>& F2) {
  return (F1, F2);
}


template <int d>
ListItemVF<d> contractProduct(const Projection& F1, const TestFunction<d>& F2) {
  return (F1, F2);
}


#endif
