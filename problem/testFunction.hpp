
#ifndef _TEST_FUNCTION_HPP
#define _TEST_FUNCTION_HPP

#include "../FESpace/expression.hpp"
#include "CutFEM_parameter.hpp"


string whichOperator(int op);
string whichOperator(int op, int cu);
string whichOperatorV(int op, int cu);

static int D(int i){
  int op[3] = {op_dx, op_dy,op_dz};
  return op[i];
}
static int D2(int i, int j){
  int op[9] = {op_dxx, op_dxy,op_dxz,
               op_dyx, op_dyy,op_dyz,
               op_dzx, op_dzy,op_dzz};
  return op[i + 3*j];
}


void f_id(RNMK_&  x, int cu, int du) ;
void f_ln(RNMK_&  x, int cu, int du) ;


template<int N = 2>
struct ItemTestFunction {
  typedef typename typeMesh<N>::Mesh Mesh;

  double c;
  int cu,du,dtu;
  KN<int> ar_nu;
  std::vector<const Virtual_Parameter*> coefu;
  int domain_id_, face_side_;
  const ExpressionVirtual * expru = nullptr;
  GFESpace<Mesh> const * fespace = nullptr;
  void (*pfun)(RNMK_&, int , int) = f_id;
  bool alloc_expr = false;


  ItemTestFunction() : c(0.), cu(-1),du(-1),dtu(-1),domain_id_(-1),face_side_(-1){}
  ItemTestFunction(double cc,int i,int j, int tu, int dd, vector<const Virtual_Parameter*> cou)
  : c(cc), cu(i),du(j),dtu(tu), domain_id_(dd), face_side_(-1){ coefu =cou;}
  ItemTestFunction(const ItemTestFunction& F)
  : c(F.c), cu(F.cu),du(F.du),dtu(F.dtu),ar_nu(F.ar_nu), domain_id_(F.domain_id_), face_side_(F.face_side_),expru(F.expru), fespace(F.fespace) {
    for(int i=0;i<F.coefu.size();++i) coefu.push_back(F.coefu[i]);
  }
  ItemTestFunction(const ItemTestFunction& F, void(*f)(RNMK_&,int,int))
  : c(F.c), cu(F.cu),du(F.du),dtu(F.dtu),ar_nu(F.ar_nu), domain_id_(F.domain_id_), face_side_(F.face_side_),expru(F.expru), fespace(F.fespace), pfun(f) {
    for(int i=0;i<F.coefu.size();++i) coefu.push_back(F.coefu[i]);
  }
  // ItemTestFunction(const FunFEM<Mesh>& ff, int ic) : c(1.), cu(ic),du(op_id),dtu(0),domain_id_(-1),fespace(ff.Vh),alloc_expr(true){
  //   expru = new ExpressionFunFEM<Mesh>(ff, ic, op_id);
  // }
  // ItemTestFunction(const GFESpace<Mesh>& vh,  const ExpressionVirtual& ff) : c(1.), cu(0),du(-1),dtu(0),domain_id_(-1),fespace(&vh){
  //   expru = &ff;
  // }

  ItemTestFunction& operator= (const ItemTestFunction& L) {
    c = L.c;
    cu = L.cu;
    du = L.du;
    dtu = L.dtu;
    ar_nu.init(L.ar_nu);
    domain_id_ = L.domain_id_;
    face_side_ = L.face_side_;
    coefu = L.coefu;
    expru = L.expru;
    fespace = L.fespace;
    assert(!L.alloc_expr);
    return *this;
  }

  void addNormal(int i) {
    int l = ar_nu.size();
    ar_nu.resize(l+1);
    ar_nu(l) = i;
  }
  void addTangent(int i) {
    int Ni = (i==0); // which normal component
    int l = ar_nu.size();
    ar_nu.resize(l+1);
    ar_nu(l) = Ni;
    if(Ni == 1) c*=-1;  //(-b,a) with (a,b) normal
  }

  void addParameter(const Virtual_Parameter& x) {
    coefu.push_back(&x);
  }
  void addParameter(const Virtual_Parameter* x) {
    coefu.push_back(x);
  }

  ItemTestFunction operator*=(R cc){
    this->c *= cc;
    return *this;
  }

  ItemTestFunction operator*=(const Normal_Component& cc){
    this->addNormal(cc.component());
    return *this;
  }


  friend std::ostream& operator << (std::ostream& f, const ItemTestFunction& u ) {
    string n[3] = {"nx", "ny", "nz"};
    f << " FESpace => " << u.fespace << "\t";
    f << to_string(u.c) << " * "
    << whichOperator(u.du, u.cu);
    for(int i=0;i<u.ar_nu.size();++i) f << " * " << n[u.ar_nu(i)];

    // for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
    f << "\t in Omega_" << u.domain_id_ ;
    f << "\n";
    return f;
  }

  ~ItemTestFunction(){
    if(alloc_expr) {
      delete expru;
    }
  }

};



template<int N = 2>
class ItemList {
  // typedef ItemTestFunction<N> ItemTestFunction;
  typedef typename typeMesh<N>::Mesh Mesh;
  typedef GFESpace<Mesh> FESpace;

public:

  KN<ItemTestFunction<N>*> U;

  ItemList() : U(1) {
    U(0) = nullptr;//new ItemTestFunction<N>();
  }
  ItemList(double cc,int i,int j, int dd=-1) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,0,dd,vector<const Virtual_Parameter*>());
  }
  ItemList(double cc,int i,int j,int tu, int dd, const vector<const Virtual_Parameter*>& cuu) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,tu,dd,cuu);
  }

  ItemList(int l) : U(l) {
    for(int i=0;i<l;++i) U(i) = new ItemTestFunction<N>();
  }

  ItemList(const ItemList& L) : U(L.U.size()) {
    for(int i=0;i<U.size();++i) U(i) = new ItemTestFunction<N>(*L.U(i));
  }

  ItemList(const ItemList& L, void(*f)(RNMK_&,int,int)) : U(L.U.size()) {
    for(int i=0;i<U.size();++i) U(i) = new ItemTestFunction<N>(*L.U(i), f);
  }

  // ItemList(const FunFEM<Mesh>& ff, int ic) : U(1) {
  //   U(0) = new ItemTestFunction<N>(ff, ic);
  // }
  // ItemList(const FESpace& Vh,const ExpressionVirtual& ff) : U(1) {
  //   U(0) = new ItemTestFunction<N>(Vh, ff);
  // }

  ItemList(const FESpace& Vh, double cc,int i,int j, int dd=-1) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,0,dd,vector<const Virtual_Parameter*>());
    U(0)->fespace = &Vh;
  }

  int size() const {return U.size();}
  const ItemTestFunction<N>& operator()(int i)const {return *U(i);}
  const ItemTestFunction<N>& getItem(int i)const {return *U(i);}
  ItemTestFunction<N>& getItem(int i) {return *U(i);}


  ItemList operator*=(R cc){
    for(int i=0;i<U.size();++i) *U(i) *= cc;
    return *this;
  }
  ItemList operator*=(const Normal_Component& cc){
    for(int i=0;i<U.size();++i) *U(i) *= cc;
    return *this;
  }


  friend std::ostream& operator <<(std::ostream& f, const ItemList & l ){
    for(int i=0;i<l.size();++i) {
      f << l(i);
    }
    return f;
  }

  ~ItemList() {
    for(int i=0;i<U.size();++i) delete U(i);
  }

};



template<int dim = 2>
class TestFunction {
  typedef typename typeRd<dim>::Rd Rd;
  typedef typename typeMesh<dim>::Mesh Mesh;
  typedef GFESpace<Mesh> FESpace;

  public :
  KNM<ItemList<dim> *> A;
private :
  TestFunction() {
    A.init(1,1);
    A(0,0) = nullptr;//new ItemList<dim> ();
  }
  TestFunction(int d) {
    A.init(d,1);
    A = nullptr;
    // for(int i=0;i<d;++i)  A(i,0) = new ItemList<dim> (1,i,0);
  }
  TestFunction(int d, int l) {
   A.init(d,l);
   A = nullptr;
 }

 TestFunction(const FESpace& Vh,int d, int comp0, int domm) {
   A.init(d,1);
   for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (Vh,1,comp0+i,0, domm);
 }


public:

  TestFunction(const FESpace& Vh,int d) {
    A.init(d,1);
    for(int i=0;i<d;++i)  A(i,0) = new ItemList<dim> (Vh,1,i,0);
  }
  TestFunction(const FESpace& Vh, int d, int comp0) {
    A.init(d,1);
    for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (Vh, 1,comp0+i,0);
  }


  TestFunction(const TestFunction& U) {
    A.init(U.A.N(), U.A.M());
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        A(i,j) = new ItemList<dim> (*U.A(i,j));
      }
    }
  }
  TestFunction(const TestFunction& U, void(*f)(RNMK_&,int,int)) {
    A.init(U.A.N(), U.A.M());
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        A(i,j) = new ItemList<dim> (*U.A(i,j), f);
      }
    }
  }
  // TestFunction(const FESpace& Vh, const ExpressionVirtual& ff) {
  //   A.init(1,1);
  //   for(int i=0;i<1;++i)  A(i,0) = new ItemList<dim> (Vh, ff);
  // }
  // TestFunction(const FESpace& Vh, const ExpressionVirtual& ff, const ExpressionVirtual& gg) {
  //   A.init(2,1);
  //   A(0,0) = new ItemList<dim> (Vh,ff);
  //   A(1,0) = new ItemList<dim> (Vh,gg);
  // }
  // TestFunction(const FESpace& Vh, const ExpressionVirtual& f1, const ExpressionVirtual& f2, const ExpressionVirtual& f3, const ExpressionVirtual& f4) {
  //   A.init(4,1);
  //   A(0,0) = new ItemList<dim> (Vh,f1);
  //   A(1,0) = new ItemList<dim> (Vh,f2);
  //   A(2,0) = new ItemList<dim> (Vh,f3);
  //   A(3,0) = new ItemList<dim> (Vh,f4);
  // }
  // TestFunction(const FunFEM<Mesh>& ff, int i0, int nb_comp) {
  //   int n = nb_comp;
  //   A.init(n,1);
  //   for(int i=i0,j=0;i<i0+nb_comp;++i,++j)  A(j,0) = new ItemList<dim> (ff, i);
  // }

  TestFunction t() const {
    TestFunction Ut(A.M(), A.N()); //Ut.init(A.M(), A.N());
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        Ut.A(j,i) = new ItemList<dim> (A(i,j)->size());
        for(int ui=0;ui<A(i,j)->size();++ui) {
          ItemTestFunction<dim>& v(A(i,j)->getItem(ui));
          ItemTestFunction<dim>& u(Ut.A(j,i)->getItem(ui));
          u = v;
        }
      }
    }
    return Ut;
  }

  // We need a line vector because N is column
  TestFunction operator*(const Normal& N) {
    // assert(!(A.M() == 1 && A.N() == dim ));// no column accepted
    assert(A.M() == dim || A.M() == 1);
    bool scalar (A.M() == 1 && A.N() == 1);
    bool line (A.M() == dim && A.N() == 1);
    bool column (A.M() == 1 && A.N() == dim);

    int s = (scalar)? dim : ((column)? 1 : A.N());

    TestFunction Un(s);
    if(scalar) {
      int ksum=0;
      ksum += A(0,0)->size();
      for(int i=0;i<dim;++i){
        int k = 0;
        Un.A(i,0) = new ItemList<dim> (ksum);
        for(int ui=0;ui<A(0,0)->size();++ui,++k) {
          ItemTestFunction<dim>& v(A(0,0)->getItem(ui));
          ItemTestFunction<dim>& u(Un.A(i,0)->getItem(k));
          u = v;
          u.addNormal(i);
        }
      }
    }
    else if(column){
      int ksum=0, k=0;
      for(int i=0;i<A.N();++i) ksum += A(i,0)->size();
      Un.A(0,0) = new ItemList<dim> (ksum);
      for(int i=0;i<A.N();++i){
        for(int ui=0;ui<A(i,0)->size();++ui,++k) {
          ItemTestFunction<dim>& v(A(i,0)->getItem(ui));
          ItemTestFunction<dim>& u(Un.A(0,0)->getItem(k));
          u = v;
          u.addNormal(i);
        }
      }
    }
    else {
      for(int i=0;i<A.N();++i){
        int ksum=0, k=0;
        for(int j=0;j<A.M();++j) ksum += A(i,j)->size();
        Un.A(i,0) = new ItemList<dim> (ksum);
        for(int j=0;j<A.M();++j){
          for(int ui=0;ui<A(i,j)->size();++ui,++k) {
            ItemTestFunction<dim>& v(A(i,j)->getItem(ui));
            ItemTestFunction<dim>& u(Un.A(i,0)->getItem(k));
            u = v;
            u.addNormal(j);
          }
        }
      }
    }
    return Un;
  }

  // We need a line vector because N is column
  TestFunction operator*(const Tangent& N) {
    // assert(!(A.M() == 1 && A.N() == dim ));// no column accepted
    assert(A.M() == dim || A.M() == 1);
    bool scalar (A.M() == 1 && A.N() == 1);
    bool line (A.M() == dim && A.N() == 1);
    bool column (A.M() == 1 && A.N() == dim);

    int s = (scalar)? dim : ((column)? 1 : A.N());

    TestFunction Un(s);
    if(scalar) {
      int ksum=0;
      ksum += A(0,0)->size();
      for(int i=0;i<dim;++i){
        int k = 0;
        Un.A(i,0) = new ItemList<dim> (ksum);
        for(int ui=0;ui<A(0,0)->size();++ui,++k) {
          ItemTestFunction<dim>& v(A(0,0)->getItem(ui));
          ItemTestFunction<dim>& u(Un.A(i,0)->getItem(k));
          u = v;
          u.addTangent(i);
        }
      }
    }
    else if(column){
      int ksum=0, k=0;
      for(int i=0;i<A.N();++i) ksum += A(i,0)->size();
      Un.A(0,0) = new ItemList<dim> (ksum);
      for(int i=0;i<A.N();++i){
        for(int ui=0;ui<A(i,0)->size();++ui,++k) {
          ItemTestFunction<dim>& v(A(i,0)->getItem(ui));
          ItemTestFunction<dim>& u(Un.A(0,0)->getItem(k));
          u = v;
          u.addTangent(i);
        }
      }
    }
    else {
      for(int i=0;i<A.N();++i){
        int ksum=0, k=0;
        for(int j=0;j<A.M();++j) ksum += A(i,j)->size();
        Un.A(i,0) = new ItemList<dim> (ksum);
        for(int j=0;j<A.M();++j){
          for(int ui=0;ui<A(i,j)->size();++ui,++k) {
            ItemTestFunction<dim>& v(A(i,j)->getItem(ui));
            ItemTestFunction<dim>& u(Un.A(i,0)->getItem(k));
            u = v;
            u.addTangent(j);
          }
        }
      }
    }
    return Un;
  }

  TestFunction operator * (const Projection& Pg ) {
    // Only for scalr right now
      assert(A.N() == 1 && A.M() == 1);
      int N = dim;
      int s = 1;
      TestFunction Un(N, N);
      // Un.init(N, N);
      const int ksum = A(0,0)->size();

      for(int i=0;i<N;++i){
        for(int j=0;j<N;++j) {
          int k = 0;
          int cst = 1+(i==j);
          Un.A(i,j) = new ItemList<dim> (cst*ksum);

          for(int ui=0;ui<ksum;++ui,++k) {
            ItemTestFunction<dim>& v(A(0,0)->getItem(ui));

            if(i==j) {
              ItemTestFunction<dim>& u(Un.A(i,j)->getItem(k));
              u = v;
              k++;
            }
            ItemTestFunction<dim>& u(Un.A(i,j)->getItem(k));
            u = v;
            u.c *= -1;
            u.addNormal(i);
            u.addNormal(j);
            k++;
          }
        }
      }
      return Un;
    }

  friend std::ostream& operator <<(std::ostream& f, const TestFunction & u ) {
    f << u.A.N() << " * " << u.A.M() << std::endl;
    for(int i=0;i<u.A.N();++i) {
      for(int j=0;j<u.A.M();++j) {
        f << *u.A(i,j);
      }
    }
    return f;
  }

  void destroy() {
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        if(A(i,j)) delete A(i,j);
      }
    }
  }

  ~TestFunction(){
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        if(A(i,j)) delete A(i,j);
      }
    }
  }


template<int N> friend TestFunction<N> grad(const TestFunction<N> & T);
template<int N> friend TestFunction<N> gradS(const TestFunction<N> & T);
template<int N> friend TestFunction<N> div(const TestFunction<N> & T);
template<int N> friend TestFunction<N> divS(const TestFunction<N> & T);
template<int N> friend TestFunction<N> divT(const TestFunction<N> & T);

template<int N> friend TestFunction<N> Eps(const TestFunction<N> & T);
template<int N> friend TestFunction<N> grad2(const TestFunction<N> & T);

template<int N> friend TestFunction<N> jump(const TestFunction<N> & T);
template<int N> friend TestFunction<N> jump(const TestFunction<N> & U, const TestFunction<N> & V );
template<int N> friend TestFunction<N> jump(const TestFunction<N> & T, int c1, int c2);
template<int N> friend TestFunction<N> jump(const TestFunction<N> & U, const TestFunction<N> & V , int c1, int c2);


template<int N> friend TestFunction<N> average(const TestFunction<N> & T, double v1, double v2);
template<int N> friend TestFunction<N> average(const TestFunction<N> & T, const Virtual_Parameter& para, const Virtual_Parameter& para2);

template<int N> friend TestFunction<N> dx(const TestFunction<N> & T);
template<int N> friend TestFunction<N> dy(const TestFunction<N> & T);
template<int N> friend TestFunction<N> dz(const TestFunction<N> & T);

template<int N> friend TestFunction<N> dxS(const TestFunction<N> & T);
template<int N> friend TestFunction<N> dyS(const TestFunction<N> & T);
template<int N> friend TestFunction<N> dzS(const TestFunction<N> & T);

template<int N> friend TestFunction<N> dt(const TestFunction<N> & T);

template<int N> friend TestFunction<N> operator - (const TestFunction<N>& F1, const TestFunction<N>& F2);
template<int N> friend TestFunction<N> operator + (const TestFunction<N>& F1, const TestFunction<N>& F2);
template<int N> friend TestFunction<N> operator , (std::list<ExpressionFunFEM<typename typeMesh<N>::Mesh>> fh, const TestFunction<N>& F2);
// template<int N> friend TestFunction<N> operator * (const CutFEM_R2& cc, const TestFunction<N>& T);

};



template <int N>
TestFunction<N> operator + (const TestFunction<N>& F1, const TestFunction<N>& F2) {
  int row = max(F1.A.N(),F2.A.N());
  int col = max(F1.A.M(),F2.A.M());

  TestFunction<N> sumU(row,col);
  for(int i=0;i<row;++i) {
    for(int j=0;j<col;++j) {
      int s = 0;
      if(i < F1.A.N() && j < F1.A.M()) s += F1.A(i,j)->size();
      if(i < F2.A.N() && j < F2.A.M()) s += F2.A(i,j)->size();
      sumU.A(i,j) = new ItemList<N>(s);
      int k = 0;
      if(i < F1.A.N() && j < F1.A.M()){
        for(int ui=0;ui<F1.A(i,j)->size();++ui,++k) {
          ItemTestFunction<N>& v(F1.A(i,j)->getItem(ui));
          ItemTestFunction<N>& u(sumU.A(i,j)->getItem(k));
          u = v;
        }
      }
      if(i < F2.A.N() && j < F2.A.M()){
        for(int ui=0;ui<F2.A(i,j)->size();++ui,++k) {
          ItemTestFunction<N>& v(F2.A(i,j)->getItem(ui));
          ItemTestFunction<N>& u(sumU.A(i,j)->getItem(k));
          u = v;
        }
      }
    }
  }
  return sumU;
}

template <int N>
TestFunction<N> operator - (const TestFunction<N>& F1, const TestFunction<N>& F2) {
  int row = max(F1.A.N(),F2.A.N());
  int col = max(F1.A.M(),F2.A.M());


  TestFunction<N> sumU(row,col); //sumU.init(row,col);
  for(int i=0;i<row;++i) {
    for(int j=0;j<col;++j) {
      int s = 0;
      if(i < F1.A.N() && j < F1.A.M()) s += F1.A(i,j)->size();
      if(i < F2.A.N() && j < F2.A.M()) s += F2.A(i,j)->size();
      sumU.A(i,j) = new ItemList<N>(s);
      int k = 0;
      if(i < F1.A.N() && j < F1.A.M()){
        for(int ui=0;ui<F1.A(i,j)->size();++ui,++k) {
          ItemTestFunction<N>& v(F1.A(i,j)->getItem(ui));
          ItemTestFunction<N>& u(sumU.A(i,j)->getItem(k));
          u = v;
        }
      }
      if(i < F2.A.N() && j < F2.A.M()){
        for(int ui=0;ui<F2.A(i,j)->size();++ui,++k) {
          ItemTestFunction<N>& v(F2.A(i,j)->getItem(ui));
          ItemTestFunction<N>& u(sumU.A(i,j)->getItem(k));
          u = v;
          u.c *= -1;
        }
      }
    }
  }
  return sumU;
}

template <int N>
TestFunction<N> operator * (const ExpressionVirtual& expr, const TestFunction<N>& F) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.expru = &expr;
      }
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator * (const TestFunction<N>& F, const ExpressionVirtual& expr) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.expru = &expr;
      }
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator, (std::list<ExpressionFunFEM<typename typeMesh<N>::Mesh>> fh, const TestFunction<N>& F) {
  assert(F.A.M() == 1);
  assert(F.A.N() == N);
  TestFunction<N> multU(1);
  if(F.A.N() != fh.size()){
    std::cout << "size expression \t" << fh.size() << std::endl;
    std::cout << "size test function \t" << F.A.N() << std::endl;
  }
  assert(F.A.N() == fh.size());

  int k = 0,ksum=0;
  for(int i=0;i<N;++i) ksum += F.A(i,0)->size();
  multU.A(0,0) = new ItemList<N>(ksum);


  auto it = fh.begin();
  for(int i=0;i<F.A.N();++i, ++it) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<F.A(i,j)->size();++ui) {
        const ItemTestFunction<N>& v(F.A(i,j)->getItem(ui));
        ItemTestFunction<N>& u(multU.A(0,0)->getItem(k));
        u = v;
        u.expru = &(*it);
        k++;
      }
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator, (const FunFEM<typename typeMesh<N>::Mesh>& fh, const TestFunction<N>& F) {
  return (fh.expression(), F);
}


template <int N>
TestFunction<N> operator * (const TestFunction<N>& F, double cc) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      *multU.A(i,j) *= cc ;
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator * (double cc, const TestFunction<N>& F) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      *multU.A(i,j) *= cc ;
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator * (const TestFunction<N>& F, const Normal_Component& c) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      *multU.A(i,j) *= c ;
    }
  }
  return multU;
}



template <int N>
TestFunction<N> ln(const TestFunction<N>& F) {
  return TestFunction<N>(F, f_ln);
}



// Need to modify that to make it inner product
// template <int N>
// TestFunction<N> operator * (const TestFunction<N>& F, typename typeRd<N>::Rd cc) {
//   assert(F.A.N() == N);
//   TestFunction<N> multU(F);
//   for(int i=0;i<F.A.N();++i) {
//     for(int j=0;j<F.A.M();++j) {
//       *multU.A(i,j) *= cc[i] ;
//     }
//   }
//   return multU;
// }
//
// template <int N>
// TestFunction<N> operator * (typename typeRd<N>::Rd cc, const TestFunction<N>& F) {
//   assert(F.A.N() == N);
//   TestFunction<N> multU(F);
//   for(int i=0;i<F.A.N();++i) {
//     for(int j=0;j<F.A.M();++j) {
//       *multU.A(i,j) *= cc[i] ;
//     }
//   }
//   return multU;
// }

template <int N>
TestFunction<N> operator * (const TestFunction<N>& F, const Virtual_Parameter& cc ) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<multU.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.addParameter(cc);
      }
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator * (const Virtual_Parameter& cc, const TestFunction<N>& F) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<multU.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.addParameter(cc);
      }
    }
  }
  return multU;
}

template <int d>
TestFunction<d> operator * (const CutFEM_Rd<d>& cc, const TestFunction<d>& T) {
  assert(T.A.M() == 1);
  int N = T.A.N();

  bool scalar = (N==1);
  int r = (scalar)?d:1;
  TestFunction<d> resU(r,1);
  if(scalar){
    int nitem = T.A(0,0)->size();

    for(int j=0;j<d;++j) {
      resU.A(j,0) = new ItemList<d>(nitem);
      for(int ui=0;ui<nitem;++ui) {
        const ItemTestFunction<d>& v(T.A(0,0)->getItem(ui));
        ItemTestFunction<d>& u(resU.A(j,0)->getItem(ui));
        u = v;
        u.addParameter(cc.get_parameter(j));
      }
    }
  }
  else {
    assert(N == d);  // column
    int nitem = 0;
    for(int i=0;i<d;++i) {
      nitem += T.A(i,0)->size();
    }
    resU.A(0,0) = new ItemList<d>(nitem);
    int k = 0;
    for(int j=0;j<d;++j) {
      int nloc = T.A(j,0)->size();
      for(int ui=0;ui<nloc;++ui) {
        const ItemTestFunction<d>& v(T.A(j,0)->getItem(ui));
        ItemTestFunction<d>& u(resU.A(0,0)->getItem(k++));
        u = v;
        u.addParameter(cc.get_parameter(j));
       }
      }
    }
    return resU;
  }






static int nextDerivative(int c, int du) {
  assert(du <= op_dz);
  if(du == op_id) return D(c);
  else return D2(du-1, c);
}

template <int d>
TestFunction<d> grad(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  bool scalar = (N==1);

  int row = scalar ? d : N;
  int col = scalar ? 1 : d;

  TestFunction<d> gradU(row, col); //gradU.init(row, col);
  for(int i=0;i<N;++i) {
    int nitem = T.A(i,0)->size();
    for(int j=0;j<d;++j) {
      int irow = scalar ? j: i;
      int jrow = scalar ? 0: j;
      gradU.A(irow,jrow) = new ItemList<d>(nitem);
      for(int ui=0;ui<nitem;++ui) {
        const ItemTestFunction<d>& v(T.A(i,0)->getItem(ui));
        ItemTestFunction<d>& u(gradU.A(irow,jrow)->getItem(ui));
        u = v;
        u.du = nextDerivative(j, v.du);//D(j);
      }
    }
  }
  return gradU;
}

template <int d>
TestFunction<d> gradS(const TestFunction<d> & T){
  assert(T.A.M() == 1);

  int N = T.A.N();
  bool scalar = (N==1);
  int col = (scalar)?1 : d;
  TestFunction<d> gradU(d, col);

  if(scalar) {
    assert(T.A(0,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(0,0)->getItem(0));

    for(int i=0;i<d;++i) {
      gradU.A(i,0) = new ItemList<d>(d+1);

      int kk=0;
      for(int j=0;j<d;++j) {
        if(i==j) {
          ItemTestFunction<d>& u(gradU.A(i,0)->getItem(kk++));
          u=v;
          // u.c = 1;
          // u.cu = i;
          u.du = D(i);
        }
        ItemTestFunction<d>& u(gradU.A(i,0)->getItem(kk++));
        u=v;
        u.c *= -1;
        // u.cu = i;
        u.du = D(j);
        u.addNormal(i);
        u.addNormal(j);
      }
    }
  }
  else {
    int n = T.A(0, 0)->size();

    for(int row=0;row<d;++row) {
      for(int col=0;col<d;++col) {
        gradU.A(row,col) = new ItemList<d>(n*(d+1));
      }
    }

    for(int row=0;row<d;++row) {
      assert(T.A(row, 0)->size() == n);

      // const ItemTestFunction<d>& v(T.A(row,0)->getItem(0));

      for(int col=0;col<d;++col) {
        int kk=0;

        for(int l=0;l<n; ++l){
          const ItemTestFunction<d>& v(T.A(row,0)->getItem(l));


          for(int j=0;j<d;++j) {
            if(col==j) {
              ItemTestFunction<d>& u(gradU.A(row,col)->getItem(kk++));
              u=v;
              u.du = D(col);
            }
            ItemTestFunction<d>& u(gradU.A(row, col)->getItem(kk++));
            u=v;
            u.c *= -1;
            u.du = D(j);
            u.addNormal(col);
            u.addNormal(j);
          }
        }
      }
    }
  }

  return gradU;
}


template <int d>
TestFunction<d> div(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  assert(T.A.N() == d);
  int N = T.A.N();

  TestFunction<d> divU(1);
  int k = 0,ksum=0;
  for(int i=0;i<N;++i) ksum += T.A(i,0)->size();
  divU.A(0,0) = new ItemList<d>(ksum);

  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    ItemTestFunction<d>& u(divU.A(0,0)->getItem(k));
    u = v;
    u.du = D(i);
    k++;
  }
  return divU;
}

template <int d>
TestFunction<d> divS(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  assert(T.A.N() == d);
  int N = T.A.N();

  TestFunction<d> divU(1);
  int k = 0,ksum=0;
  for(int i=0;i<N;++i) ksum += T.A(i,0)->size();
  divU.A(0,0) = new ItemList<d>(ksum*(d+1));

  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    ItemTestFunction<d>& u(divU.A(0,0)->getItem(k));
    u=v;
    u.du = D(i);
    k++;

    for(int j=0;j<d;++j) {
      ItemTestFunction<d>& u(divU.A(0,0)->getItem(k));
      u=v;
      u.c *= -1;
      u.du = D(j);
      u.addNormal(i);
      u.addNormal(j);
      k++;
    }
  }
  return divU;
}

template <int d>
TestFunction<d> divT(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  assert(T.A.N() == d);
  assert(d == 2);
  int N = T.A.N();

  TestFunction<d> divU(1);
  int k = 0,ksum=0;
  for(int i=0;i<N;++i) ksum += T.A(i,0)->size();
  divU.A(0,0) = new ItemList<d>(ksum*(d));

  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));


    for(int j=0;j<d;++j) {
      int ii = (i==0);
      int jj = (j==0);
      ItemTestFunction<d>& u(divU.A(0,0)->getItem(k));
      u=v;
      u.c *= (ii==jj)?1 : -1;
      u.du = D(j);
      u.addNormal(jj);
      u.addNormal(ii);
      k++;
    }
  }
  return divU;
}

// time gradient
template <int d>
TestFunction<d> dt(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  int row = N;
  int col = 1;

  TestFunction<d> gradU(row, col); //gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);
    ItemTestFunction<d>& u(gradU.A(irow,jrow)->getItem(0));
    u = v;
    u . dtu = D(0);
  }
  return gradU;
}

template <int d>
TestFunction<d> dx(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  int row = N;
  int col = 1;

  TestFunction<d> gradU(row, col);// gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);
    ItemTestFunction<d>& u(gradU.A(irow,jrow)->getItem(0));
    u = v;
    u . du = D(0);
  }
  return gradU;
}

template <int d>
TestFunction<d> dy(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  int row = N;
  int col = 1;

  TestFunction<d> gradU(row, col); //gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);
    ItemTestFunction<d>& u(gradU.A(irow,jrow)->getItem(0));
    u = v;
    u . du = D(1);
  }
  return gradU;
}

template <int d>
TestFunction<d> dz(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  int row = N;
  int col = 1;

  TestFunction<d> gradU(row, col); //gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);
    ItemTestFunction<d>& u(gradU.A(irow,jrow)->getItem(0));
    u = v;
    u . du = D(2);
  }
  return gradU;
}

template <int d>
TestFunction<d> dxS(const TestFunction<d> & T){
  assert(T.A.M() == 1 && T.A.N() == 1);
  int N = T.A.N();

  TestFunction<d> gradU(N,1); // gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);
    {
      ItemTestFunction<d>& u1(gradU.A(i,0)->getItem(0));
      u1 = v;
      u1.du = D(0);
    }
    int kk=1;
    for(int j=0;j<d;++j) {
      ItemTestFunction<d>& u(gradU.A(i,0)->getItem(kk++));
      u=v;
      u.c *= -1;
      u.du = D(j);
      u.addNormal(0);
      u.addNormal(j);
    }
  }
  return gradU;
}

template <int d>
TestFunction<d> dyS(const TestFunction<d> & T){
  assert(T.A.M() == 1 && T.A.N() == 1);
  int N = T.A.N();

  TestFunction<d> gradU(N,1); // gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);
    {
      ItemTestFunction<d>& u1(gradU.A(i,0)->getItem(0));
      u1 = v;
      u1.du = D(1);
    }
    int kk=1;
    for(int j=0;j<d;++j) {
      ItemTestFunction<d>& u(gradU.A(i,0)->getItem(kk++));
      u=v;
      u.c *= -1;
      u.du = D(j);
      u.addNormal(1);
      u.addNormal(j);
    }
  }
  return gradU;
}

template <int d>
TestFunction<d> dzS(const TestFunction<d> & T){
  assert(T.A.M() == 1 && T.A.N() == 1);
  int N = T.A.N();

  TestFunction<d> gradU(N,1); //gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);
    {
      ItemTestFunction<d>& u1(gradU.A(i,0)->getItem(0));
      u1 = v;
      u1.du = D(2);
    }
    int kk=1;
    for(int j=0;j<d;++j) {
      ItemTestFunction<d>& u(gradU.A(i,0)->getItem(kk++));
      u=v;
      u.c *= -1;
      u.du = D(j);
      u.addNormal(2);
      u.addNormal(j);
    }
  }
  return gradU;
}



template <int d>
TestFunction<d> Eps(const TestFunction<d> & T){
  assert(T.A.N() == d && T.A.M() == 1);

  TestFunction<d> epsU(d, d); //epsU.init(d, d);
  for(int i=0;i<d;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    for(int j=0;j<d;++j) {
      if(i==j){
        epsU.A(i,j) = new ItemList<d>(1);
        ItemTestFunction<d>& u(epsU.A(i,j)->getItem(0));
        u = v;
        u.du = D(i);

      }
      else  {
        epsU.A(i,j) = new ItemList<d>(2);
        {
          ItemTestFunction<d>& u(epsU.A(i,j)->getItem(0));
          u=v;
          u.c = 0.5*v.c;
          u.cu = i;
          u.du = D(j);
        }
        {
          ItemTestFunction<d>& u(epsU.A(i,j)->getItem(1));
          u=v;
          u.c = 0.5*v.c;
          u.cu = j;
          u.du = D(i);
        }
      }
    }
  }
  return epsU;
}
//
//
template <int d>
TestFunction<d> grad2(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  assert(T.A.N() == d || T.A.N() == 1);
  int N = T.A.N();

  TestFunction<d> DDU(T.A.N(), T.A.M()); //DDU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) assert(T.A(i,0)->size()==1);


  for(int i=0;i<N;++i) {
    DDU.A(i,0) = new ItemList<d>(3);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int k = 0;
    for(int i1=0;i1<d;++i1) {
      for(int i2=i1;i2<d;++i2){
      // 0,0-> dxx, 0,1->2dxy  1,1->dyy
        ItemTestFunction<d>& u(DDU.A(i,0)->getItem(k));
        u = v;
        u.du = D2(i1, i2);
        u.c = (i1 != i2) + 1;
        k++;
      }
    }
  }
  return DDU;
}


template <int d>
TestFunction<d> jump(const TestFunction<d> & T, int c1, int c2){
  // assert(T.A.M() == 1);
  int N = T.A.N();
  int M = T.A.M();
  TestFunction<d> jumpU(T.A.N(), T.A.M()); //jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    for(int j=0;j<M;++j) {
      int l = T.A(i,j)->size();
      jumpU.A(i,j) = new ItemList<d>(2*l);
      for(int e=0;e<l;++e) {
        const ItemTestFunction<d>& v(T.A(i,j)->getItem(e));
        {
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e));
          u = v;
          u.c *= c1;
          u.face_side_ = 0;
        }
        {
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e+1));
          u=v;
          u.c *= c2;
          u.face_side_ = 1;
        }
      }
    }
  }
  return jumpU;
}
template <int d>
TestFunction<d> jump(const TestFunction<d> & T){
  return jump(T, 1, -1);
}

// NEED DO FIXE THIS FUNCTION
// HAS TO WORK FOR GENERAL DOMAIN
// template <int d>
// TestFunction<d> jump(const TestFunction<d> & U, const TestFunction<d> & V, int c1, int c2){
//   assert(U.A.M() == 1 && U.A.N() == 1);
//   assert(V.A.M() == 1 && V.A.N() == 1);
//   int N = U.A.N();
//   int M = U.A.M();
//   TestFunction<d> jumpU(U.A.N(), U.A.M()); //jumpU.init(U.A.N(), U.A.M());
//   for(int i=0;i<N;++i) {
//     for(int j=0;j<M;++j) {
//       assert(U.A(i,j)->size() == V.A(i,j)->size());
//       int l = U.A(i,j)->size();
//       jumpU.A(i,j) = new ItemList<d>(2*l);
//       for(int e=0;e<l;++e) {
//         {
//           const ItemTestFunction<d>& v(U.A(i,j)->getItem(e));
//           ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e));
//           u = v;
//           u.c *= c1;
//           // u.dom = 0;
//         }
//         {
//           const ItemTestFunction<d>& v(V.A(i,j)->getItem(e));
//           ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e+1));
//           u=v;
//           u.c *= c2;
//           // u.dom = 1;
//         }
//       }
//     }
//   }
//   return jumpU;
// }

template <int d>
TestFunction<d> jump(const TestFunction<d> & U, const TestFunction<d> & V){
  return jump(U, V, 1, -1);
}



// template <int d>
// TestFunction<d> average1(const TestFunction<d> & T){
//   assert(T.A.M() == 1);
//   int N = T.A.N();
//   TestFunction<d> jumpU(T.A.N(), T.A.M()); //jumpU.init(T.A.N(), T.A.M());
//   for(int i=0;i<N;++i) {
//     int l = T.A(i,0)->size();
//     jumpU.A(i,0) = new ItemList<d>(2*l);
//     for(int e=0;e<l;++e) {
//       const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
//       {
//         ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
//         u = v;
//         u.face_side_ = 0;
//         u.addParameter("kappa1");
//       }
//       {
//         ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
//         u=v;
//         u.face_side_ = 1;
//         u.addParameter("kappa2");
//       }
//     }
//   }
//   return jumpU;
// }

// template <int d>
// TestFunction<d> average2(const TestFunction<d> & T){
//   assert(T.A.M() == 1);
//   int N = T.A.N();
//   TestFunction<d> jumpU(T.A.N(), T.A.M()); //jumpU.init(T.A.N(), T.A.M());
//   for(int i=0;i<N;++i) {
//     int l = T.A(i,0)->size();
//     jumpU.A(i,0) = new ItemList<d>(2*l);
//     for(int e=0;e<l;++e) {
//       const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
//       {
//         ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
//         u = v;
//         u.face_side_ = 0;
//         u.addParameter("kappa2");
//       }
//       {
//         ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
//         u=v;
//         u.face_side_ = 1;
//         u.addParameter("kappa1");
//       }
//     }
//   }
//   return jumpU;
// }

template <int d>
TestFunction<d> average(const TestFunction<d> & T, double v1=0.5, double v2=0.5){
  assert(T.A.M() == 1);
  int N = T.A.N();
  TestFunction<d> jumpU(T.A.N(), T.A.M()); //jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    int l = T.A(i,0)->size();
    jumpU.A(i,0) = new ItemList<d>(2*l);
    for(int e=0;e<l;++e) {
      const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
        u = v;
        u.c *= v1;
        u.face_side_ = 0;
      }
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
        u=v;
        u.c *= v2;
        u.face_side_ = 1;
      }
    }
  }
  return jumpU;
}

template <int d>
TestFunction<d> average(const TestFunction<d> & T, const Virtual_Parameter& para1, const Virtual_Parameter& para2){
  assert(T.A.M() == 1);
  int N = T.A.N();
  TestFunction<d> jumpU(T.A.N(), T.A.M()); //jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    int l = T.A(i,0)->size();
    jumpU.A(i,0) = new ItemList<d>(2*l);
    for(int e=0;e<l;++e) {
      const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
        u = v;
        u.face_side_ = 0;
        u.addParameter(para1);
      }
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
        u=v;
        u.face_side_ = 1;
        u.addParameter(para2);
      }
    }
  }
  return jumpU;
}


typedef TestFunction<2> TestFunction2;
typedef TestFunction<3> TestFunction3;

#endif
