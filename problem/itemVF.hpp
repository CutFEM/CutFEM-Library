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
  KN<int> conormalU_, conormalV_;
  int face_sideU_, face_sideV_;
  int domainU_id_, domainV_id_;

  std::vector<const Virtual_Parameter*> coefu, coefv ;
  int dtu, dtv;
  const ExpressionVirtual* expru=nullptr;
  const ExpressionVirtual* exprv=nullptr;

  FESpace const * fespaceU = nullptr;
  FESpace const * fespaceV = nullptr;

  void(*pfunU)(RNMK_&,int,int) = f_id;
  void(*pfunV)(RNMK_&,int,int) = f_id;

  // CutFEM_ParameterList parameterList;


  ItemVF()
  : c(0.), cu(-1),du(-1),cv(-1),dv(-1),face_sideU_(-1),face_sideV_(-1),domainU_id_(-1), domainV_id_(-1), dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l)
  : c(cc), cu(i),du(j),cv(k),dv(l),face_sideU_(-1),face_sideV_(-1),domainU_id_(-1), domainV_id_(-1),dtu(-1),dtv(-1){}
  ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm)
  : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),face_sideU_(-1),face_sideV_(-1),domainU_id_(-1), domainV_id_(-1), dtu(-1),dtv(-1){}
  // ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm, int dou, int dov)
  // : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),face_sideU_(dou),face_sideV_(dov),dtu(-1),dtv(-1){}
  // ItemVF(double cc,int i,int j,int k,int l,const KN<int>& nn, const KN<int>& mm, int dou, int dov,
  //   const vector<const Virtual_Parameter*>& weu, const vector<const Virtual_Parameter*>& wev, int tu, int tv)
  // : c(cc), cu(i),du(j),cv(k),dv(l),ar_nu(nn),ar_nv(mm),face_sideU_(dou),face_sideV_(dov),dtu(tu),dtv(tv){
  //   for(int i=0;i<weu.size();++i) coefu.push_back(weu[i]);
  //   for(int i=0;i<wev.size();++i) coefv.push_back(wev[i]);
  // }
  ItemVF(const ItemVF & U)
  :ItemVF(U.c,U.cu,U.du,U.cv,U.dv,U.ar_nu,U.ar_nv){
    conormalU_ = U.conormalU_;
    conormalV_ = U.conormalV_;

    face_sideU_ = U.face_sideU_;
    face_sideV_ = U.face_sideV_;

    domainU_id_ = U.domainU_id_;
    domainV_id_ = U.domainV_id_;

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
    :ItemVF(U.c*V.c,U.cu,U.du,V.cu,V.du,U.ar_nu,V.ar_nu) {
    conormalU_ = U.conormal;
    conormalV_ = V.conormal;

    face_sideU_ = U.face_side_;
    face_sideV_ = V.face_side_;

    domainU_id_ = U.domain_id_;
    domainV_id_ = V.domain_id_,

    coefu = U.coefu;
    coefv = V.coefu;
    dtu = U.dtu;
    dtv = V.dtu;

    expru = U.expru;
    exprv = V.expru;

    fespaceU = U.fespace;
    fespaceV = V.fespace;

    pfunU = U.pfun;
    pfunV = V.pfun;
  }


  // // bool same() const {
  //   return (fespaceU == fespaceV) && (pfunU == pfunV); }

  bool operator==(const ItemVF& F){
    if(cu == F.cu && cv == F.cv && du == F.du && dv == F.dv
        && F.face_sideU_ == face_sideU_ && face_sideV_ == F.face_sideV_ && dtu == F.dtu && dtv == F.dtv
        && domainU_id_ == F.domainU_id_ && domainV_id_ == F.domainV_id_
        && fespaceU == F.fespaceU && fespaceV == F.fespaceV
        && expru == F.expru && exprv == F.exprv) {
      if(ar_nu.size() == F.ar_nu.size() && ar_nv.size() == F.ar_nv.size()) {
        for(int i=0;i<ar_nu.size();++i) {
          if(ar_nu(i) != F.ar_nu(i)) return false;
        }
        for(int i=0;i<ar_nv.size();++i) {
          if(ar_nv(i) != F.ar_nv(i)) return false;
        }
      }
      if(conormalU_.size() == F.conormalU_.size() && conormalV_.size() == F.conormalV_.size()) {
        for(int i=0;i<conormalU_.size();++i) {
          if(conormalU_(i) != F.conormalU_(i)) return false;
        }
        for(int i=0;i<conormalV_.size();++i) {
          if(conormalV_(i) != F.conormalV_(i)) return false;
        }
      }
      return true;
    }
    return false;
  }


  // double fxU(int k, Rd mip, const R* normal = nullptr) const {
  //   return ((expru)? expru->eval(k, mip, normal) : 1);
  // }
  // double fxV(int k, Rd mip, const R* normal = nullptr) const {
  //   return ((exprv)? exprv->eval(k, mip, normal) : 1);
  // }
  // double fxu(int k, Rd mip) const {
  //   return ((expru)? expru->eval(k, mip) : 1)*((exprv)? exprv->eval(k, mip) : 1);
  // }
  // double fxu(int k, Rd mip, double t) const {
  //   return ((expru)? expru->eval(k, mip, t) : 1)*((exprv)? exprv->eval(k, mip,t) : 1);
  // }
  //
  // double fxu_backMesh(int k, int dom, Rd mip, const R* normal = nullptr) const {
  //   return ((expru)? expru->GevalOnBackMesh(k, dom, mip, normal) : 1)
  //         *((exprv)? exprv->GevalOnBackMesh(k, dom, mip, normal) : 1);
  // }
  // double fxu_backMesh(int k, int dom, Rd mip, double t, const R* normal = nullptr) const {
  //   return ((expru)? expru->GevalOnBackMesh(k, dom, mip, t, normal) : 1)
  //          *((exprv)? exprv->GevalOnBackMesh(k, dom, mip, t, normal) : 1);
  // }
  //
  // double fx_backMesh_U(int k, int dom, Rd mip, const R* normal = nullptr) const {
  //   return ((expru)? expru->GevalOnBackMesh(k, dom, mip, normal) : 1);
  // }
  // double fx_backMesh_V(int k, int dom, Rd mip, const R* normal = nullptr) const {
  //   return ((exprv)? exprv->GevalOnBackMesh(k, dom, mip, normal) : 1);
  // }


public:
  void applyFunNL(RNMK_& bfu, RNMK_& bfv) const {
    pfunU(bfu, cu, du);
    pfunV(bfv, cv, dv);
  }

  // FOR NEW VERSION
  double computeCoefElement(double h, double meas, double measK, double measCut, int domain) const {
    R val = 1;
    for(int l=0;l<2;++l) {
      const vector<const Virtual_Parameter*>& listCoef = (l==0)?coefu : coefv;
      for(int i=0;i<listCoef.size();++i) {
        val *= listCoef[i]->evaluate(domain, h, meas, measK, measCut);
      }
    }
    return val;
  }
  double computeCoefInterface(double h, double meas, double measK, double measCut, int domi, int domj) const {
    R val = 1;
    for(int l=0;l<2;++l) {
      const vector<const Virtual_Parameter*>& listCoef = (l==0)?coefu : coefv;
      int dom = (l==0)?domi:domj;
      for(int i=0;i<listCoef.size();++i) {
        val *= listCoef[i]->evaluate(dom, h, meas, measK, measCut);
      }
    }
    return val;
  }
  double evaluateFunctionOnBackgroundMesh(int k, int dom, Rd mip, const R* normal = nullptr) const {
    return ((expru)? expru->GevalOnBackMesh(k, dom, mip, normal) : 1)
    *((exprv)? exprv->GevalOnBackMesh(k, dom, mip, normal) : 1);
  }
  double evaluateFunctionOnBackgroundMesh(const std::pair<int,int>& k, const std::pair<int,int>& dom, Rd mip, const R* normal = nullptr) const {
    return ((expru)? expru->GevalOnBackMesh(k.first , dom.first , mip, normal) : 1)
    *((exprv)? exprv->GevalOnBackMesh(k.second, dom.second, mip, normal) : 1);
  }
  double evaluateFunctionOnBackgroundMesh(const std::pair<int,int>& k, const std::pair<int,int>& dom, Rd mip, const std::pair<Rd, Rd>normal) const {
    return ((expru)? expru->GevalOnBackMesh(k.first , dom.first , mip, normal.first) : 1)
    *((exprv)? exprv->GevalOnBackMesh(k.second, dom.second, mip, normal.second) : 1);
  }
  double computeCoefFromNormal(const Rd normal) const {
    return computeCoeffFromArray(ar_nu, normal) * computeCoeffFromArray(ar_nv,normal);;
  }
  double computeCoefFromConormal(const Rd conormal) const {
    return computeCoeffFromArray(conormalU_, conormal) * computeCoeffFromArray(conormalV_, conormal);
  }
  double computeCoefFromNormal(const Rd normal0,const Rd normal1) const {
    return computeCoeffFromArray(ar_nu, normal0) * computeCoeffFromArray(ar_nv,normal1);;
  }
  double computeCoefFromConormal(const Rd conormalu, const Rd conormalv) const {
    return computeCoeffFromArray(conormalU_, conormalu) * computeCoeffFromArray(conormalV_, conormalv);
  }
  int onWhatElementIsTrialFunction (std::vector<int> k) const {
    if(k.size()==1) {return k[0];}
    else {assert(face_sideU_ != -1); return (face_sideU_==0)? k[0] : k[1];}
  }
  int onWhatElementIsTestFunction  (std::vector<int> k) const {
    if(k.size()==1) { return k[0];}
    else {assert(face_sideV_ != -1); return (face_sideV_==0)? k[0] : k[1];}
  }
  template<typename T>
  T onWhatElementIsTrialFunction (T ki, T kj) const {
    assert(face_sideU_ != -1); return (face_sideU_==0)? ki : kj;
  }
  template<typename T>
  T onWhatElementIsTestFunction  (T ki, T kj) const {
    assert(face_sideV_ != -1); return (face_sideV_==0)? ki : kj;
  }
  int get_domain_test_function() const {return domainV_id_;}
  int get_domain_trial_function() const {return domainU_id_;}


  bool on(int d) const {
    return ((domainU_id_ == domainV_id_) && (domainU_id_ == -1 || domainU_id_ == d));
  }

private:
  double computeCoeffFromArray(const KN<int>& array_idx, const R* v) const {
    R val = 1;
    for(int i=0;i<array_idx.size();++i) val *= v[array_idx(i)];
    return val;
  }

public:



// ------------------------------------------

friend std::ostream& operator <<(std::ostream& f, const ItemVF & u )
{
  string n[3] = {"nx", "ny", "nz"};
  f << " FESpaces => " << u.fespaceU << " and " << u.fespaceV << "\t";
  f << u.c << "\t" << whichOperator( u.dtu) << whichOperator( u.du,u.cu);
  for(int i=0;i<u.ar_nu.size();++i) f << " * " << n[u.ar_nu(i)];
  // for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
  f << " * " << whichOperator( u.dtv) << whichOperatorV( u.dv,u.cv);
  for(int i=0;i<u.ar_nv.size();++i) f << " * " << n[u.ar_nv(i)];
  // for(int i=0;i<u.coefv.size();++i) f << " * " << u.coefv[i];
  if(u.face_sideU_ == u.face_sideV_ && u.face_sideU_ != -1) f << "\t in Omega_" << u.face_sideU_+1;

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

  typedef typename typeMesh<N>::Mesh Mesh;
  typedef GFESpace<Mesh> FESpace;
  public :
  KN<ItemVF<N>> VF;
  bool isRHS_ = true;
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

  int get_lastOp() const {
    int n = 0;
    for(int i=0; i<this->size();++i) n = Max(n, getLastop(VF[i].du, VF[i].dv));
    return n;
  }
  bool isRHS() const { return isRHS_;}
  const FESpace& get_spaceU(int i) const {assert(VF(i).fespaceU||isRHS_);return (VF(i).fespaceU)?*VF(i).fespaceU:*VF(i).fespaceV;}
  const FESpace& get_spaceV(int i) const {assert(VF(i).fespaceV);        return *VF(i).fespaceV;}


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
    item(i).face_sideU_ = 0; item(i).face_sideV_ = 0;
  }
  for(int i=n0,j=0;i<n;++i,++j) {
    item(i) = L(j);
    item(i).c *= -1;
    item(i).face_sideU_ = 1; item(i).face_sideV_ = 1;
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
    item(i).face_sideU_ = 0; item(i).face_sideV_ = 0;
  }
  for(int i=n0,j=0;i<n;++i,++j) {
    item(i) = L(j);
    item(i).c *= v2;
    item(i).face_sideU_ = 1; item(i).face_sideV_ = 1;
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
  item.isRHS_ = false;
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
        item(k) = ItemVF<d>( v.c*c,v.cu,0,v.cu,v.du,0,v.ar_nu);
        item(k).face_sideU_ = v.face_side_;
        item(k).face_sideV_ = v.face_side_;
        item(k).coefv = v.coefu;
        item(k).dtu = -1;
        item(k).dtv = v.dtu;
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
        item(k) = ItemVF<d>( v.c*c(i,j),v.cu,0,v.cu,v.du,0,v.ar_nu);
        item(k).face_sideU_ = v.face_side_;
        item(k).face_sideV_ = v.face_side_;
        item(k).coefv = v.coefu;
        item(k).dtu = -1;
        item(k).dtv = v.dtu;
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
          item(k) = ItemVF<d>( v.c,v.cu,0,v.cu,v.du,0,v.ar_nu,v.face_side_,v.face_side_,vector<const Virtual_Parameter*>(), v.coefu,-1,v.dtu);
          item(k).exprv = v.expru;
          item(k).fespaceV = v.fespace;
          k++;
        }
        item(k) = ItemVF<d>( v.c*(-1),v.cu,0,v.cu,v.du,c(i,j),v.ar_nu);
        item(k).face_sideU_ = v.face_side_;
        item(k).face_sideV_ = v.face_side_;
        item(k).coefv = v.coefu;
        item(k).dtu = -1;
        item(k).dtv = v.dtu;

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
        item(k) = ItemVF<d>( v.c,0,-1,v.cu,v.du,0,v.ar_nu);
        item(k).face_sideU_ = v.face_side_;
        item(k).face_sideV_ = v.face_side_;
        item(k).coefv = v.coefu;
        item(k).dtu = 0;
        item(k).dtv = v.dtu;
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
        item(k) = ItemVF<d>( v.c,0,-1,v.cu,v.du,0,v.ar_nu);
        item(k).face_sideU_ = v.face_side_;
        item(k).face_sideV_ = v.face_side_;
        item(k).coefv = v.coefu;
        item(k).dtu = 0;
        item(k).dtv = v.dtu;
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
