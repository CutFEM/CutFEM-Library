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





template<int N = 2>
struct ItemTestFunction {
  typedef typename typeMesh<N>::Mesh Mesh;

  double c;
  int cu,du,dtu;
  KN<int> ar_nu;
  vector<string> coefu;
  int dom;
  const ExpressionVirtual * expru = nullptr;
  GFESpace<Mesh> const * fespace = nullptr;;


  ItemTestFunction() : c(0.), cu(-1),du(-1),dtu(-1),dom(-1){}
  ItemTestFunction(double cc,int i,int j, int tu, int dd, vector<string> cou)
  : c(cc), cu(i),du(j),dtu(tu), dom(dd){ coefu =cou;}
  // ItemTestFunction(double cc,int i,int j, int dd)
  // : c(cc), cu(i),du(j),dtu(-1), dom(dd){}
  ItemTestFunction(const ItemTestFunction& F)
  : c(F.c), cu(F.cu),du(F.du),dtu(F.dtu),ar_nu(F.ar_nu), dom(F.dom), expru(F.expru), fespace(F.fespace) {
    for(int i=0;i<F.coefu.size();++i) coefu.push_back(F.coefu[i]);
  }

  ItemTestFunction& operator = (const ItemTestFunction& L) {
    c = L.c;
    cu = L.cu;
    du = L.du;
    dtu = L.dtu;
    ar_nu.init(L.ar_nu);
    dom = L.dom;
    coefu = L.coefu;
    expru = L.expru;
    fespace = L.fespace;
    return *this;

  }

  void addNormal(int i) {
    int l = ar_nu.size();
    ar_nu.resize(l+1);
    ar_nu(l) = i;
  }
  void addTangent(int Ni) {
    int l = ar_nu.size();
    ar_nu.resize(l+1);
    ar_nu(l) = Ni;
    if(Ni == 1) c*=-1;  //(-b,a) with (a,b) normal
  }

  void addParameter(string x) {
    coefu.push_back(x);
  }

  ItemTestFunction operator*=(R cc){
    this->c *= cc;
    return *this;
  }


  friend std::ostream& operator <<(std::ostream& f, const ItemTestFunction & u )
  {

    string n[3] = {"nx", "ny", "nz"};
    f << " FESpace => " << u.fespace << "\t";
    f << to_string(u.c) << " * "
      << whichOperator(u.du, u.cu);
    for(int i=0;i<u.ar_nu.size();++i) f << " * " << n[u.ar_nu(i)];
    // f << "\n";
    for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
    if(u.dom != -1) f << "\t in Omega_" << u.dom+1;
    else f << "\t in Omega";
    f << "\n";
    return f;
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
    U(0) = new ItemTestFunction<N>();
  }
  ItemList(double cc,int i,int j, int dd=-1) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,0,dd,vector<string>());
  }
  ItemList(double cc,int i,int j,int tu, int dd, const vector<string>& cuu) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,tu,dd,cuu);
  }

  ItemList(int l) : U(l) {
    for(int i=0;i<l;++i) U(i) = new ItemTestFunction<N>();
  }

  ItemList(const ItemList& L) : U(L.U.size()) {
    for(int i=0;i<U.size();++i) U(i) = new ItemTestFunction<N>(*L.U(i));
  }

  ItemList(const FESpace& Vh, double cc,int i,int j, int dd=-1) : U(1) {
    U(0) = new ItemTestFunction<N>(cc,i,j,0,dd,vector<string>());
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
  // typedef ItemTestFunction<dim> ItemTestFunction;
  // typedef  ItemList<dim> ItemList;
  typedef GFESpace<Mesh> FESpace;

  public :
  KNM<ItemList<dim> *> A;
private :
  TestFunction() {
    A.init(1,1);
    A(0,0) = new ItemList<dim> ();
  }
  TestFunction(int d) {
    A.init(d,1);
    for(int i=0;i<d;++i)  A(i,0) = new ItemList<dim> (1,i,0);
  }

  // TestFunction(int d, int comp0) {
  //   A.init(d,1);
  //   for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (1,comp0+i,0);
  // }
  // TestFunction(int d, int comp0, int domm) {
  //   A.init(d,1);
  //   for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (1,comp0+i,0, domm);
  // }

public:

  TestFunction(const FESpace& Vh,int d) {
    A.init(d,1);
    for(int i=0;i<d;++i)  A(i,0) = new ItemList<dim> (Vh,1,i,0);
  }
  TestFunction(const FESpace& Vh, int d, int comp0) {
    A.init(d,1);
    for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (Vh, 1,comp0+i,0);
  }
  TestFunction(const FESpace& Vh,int d, int comp0, int domm) {
    A.init(d,1);
    for(int i=0;i<d;++i) A(i,0) = new ItemList<dim> (Vh,1,comp0+i,0, domm);
  }


  TestFunction(const TestFunction& U) {
    A.init(U.A.N(), U.A.M());
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        A(i,j) = new ItemList<dim> (*U.A(i,j));
        // A(i,j)->fespace = U.A(i,j)->fespace;
      }
    }
    // A.fespace = U.fespace;
  }

  void init(int d, int l) {
    destroy();
    A.init(d,l);
    for(int i=0;i<d;++i) {
      for(int j=0;j<l;++j) {
        A(i,j) = new ItemList<dim> (1,i,0);
      }
    }
  }

  TestFunction t() const {
    TestFunction Ut; Ut.init(A.M(), A.N());
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
          u.addNormal(N[i]);
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
          u.addNormal(N[i]);
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
            u.addNormal(N[j]);
          }
        }
      }
    }
    return Un;
  }



  // We need a line vector because N is column
  TestFunction operator*(const Tangent& T) {
    assert(A.M() == dim || A.M() == 1);
    bool scalar (A.M() == 1 && A.N() == 1);
    int s = (scalar)? dim : A.N();
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
          u=v;
          u.addTangent(T[i]);
          // u.fespace = v.fespace;
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
            u=v;
            u.addTangent(T[j]);
            // u.fespace = v.fespace;
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
      TestFunction Un;
      Un.init(N, N);
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


  friend std::ostream& operator <<(std::ostream& f, const TestFunction & u )
  {
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
        delete A(i,j);
      }
    }
  }

  ~TestFunction(){
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
        delete A(i,j);
      }
    }
  }



// friend A foo<T>(A& a);

template<int N> friend TestFunction<N> grad(const TestFunction<N> & T);
template<int N> friend TestFunction<N> gradS(const TestFunction<N> & T);
template<int N> friend TestFunction<N> div(const TestFunction<N> & T);
template<int N> friend TestFunction<N> Eps(const TestFunction<N> & T);
template<int N> friend TestFunction<N> grad2(const TestFunction<N> & T);

template<int N> friend TestFunction<N> jump(const TestFunction<N> & T);
template<int N> friend TestFunction<N> jump(const TestFunction<N> & U, const TestFunction<N> & V );

template<int N> friend TestFunction<N> average(const TestFunction<N> & T);
template<int N> friend TestFunction<N> average1(const TestFunction<N> & T);
template<int N> friend TestFunction<N> average2(const TestFunction<N> & T);

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
template<int N> friend TestFunction<N> operator * (const CutFEM_R2& cc, const TestFunction<N>& T);

};

















template <int N>
TestFunction<N> operator + (const TestFunction<N>& F1, const TestFunction<N>& F2) {
  int row = max(F1.A.N(),F2.A.N());
  int col = max(F1.A.M(),F2.A.M());


  TestFunction<N> sumU; sumU.init(row,col);
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


  TestFunction<N> sumU; sumU.init(row,col);
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
TestFunction<N> operator , (std::list<ExpressionFunFEM<typename typeMesh<N>::Mesh>> fh, const TestFunction<N>& F) {
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
TestFunction<N> operator * (const TestFunction<N>& F, const CutFEM_Parameter& cc ) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<multU.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.addParameter(cc.getName());
      }
    }
  }
  return multU;
}

template <int N>
TestFunction<N> operator * (const CutFEM_Parameter& cc, const TestFunction<N>& F) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<multU.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        v.addParameter(cc.getName());
      }
    }
  }
  return multU;
}

template <int d>
TestFunction<d> operator * (const CutFEM_R2& cc, const TestFunction<d>& T) {
  assert(T.A.M() == 1);
  int N = T.A.N();

  bool scalar = (N==1);
  TestFunction<d> resU;
  if(scalar){
    resU.init(d,1);
    int nitem = T.A(0,0)->size();

    for(int j=0;j<d;++j) {
      resU.A(j,0) = new ItemList<d>(nitem);
      for(int ui=0;ui<nitem;++ui) {
        const ItemTestFunction<d>& v(T.A(0,0)->getItem(ui));
        ItemTestFunction<d>& u(resU.A(j,0)->getItem(ui));
        u = v;
        u.addParameter(cc.getName(j));
      }
    }
  }
  else {
    assert(N == d);  // column
    resU.init(1,1);
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
        u.addParameter(cc.getName(j));
       }
      }
    }
    return resU;
  }


template <int N>
TestFunction<N> operator * (Pow_Par cc, const TestFunction<N>& F) {
  TestFunction<N> multU(F);
  for(int i=0;i<F.A.N();++i) {
    for(int j=0;j<F.A.M();++j) {
      for(int ui=0;ui<multU.A(i,j)->size();++ui) {
        ItemTestFunction<N>& v(multU.A(i,j)->getItem(ui));
        for(int ll=0;ll<cc.list_name.size();++ll){
          v.addParameter(cc.list_name[ll]);
        }
        for(int ll=0;ll<cc.list_cst.size();++ll){
          v.c *= cc.list_cst[ll];
        }
      }
    }
  }
  return multU;
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

  TestFunction<d> gradU; gradU.init(row, col);
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
  TestFunction<d> gradU;

  if(scalar) {
    gradU.init(d, 1);
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
    gradU.init(d, d);
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

// time gradient
template <int d>
TestFunction<d> dt(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();

  int row = N;
  int col = 1;

  TestFunction<d> gradU; gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);//(v.c,i,v.du,D(0),v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);//(v.c,i,D(0),v.dtu,v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);//(v.c,i,D(1),v.dtu,v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(row, col);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    int irow = i;
    int jrow = 0;
    gradU.A(irow,jrow) = new ItemList<d>(1);//(v.c,i,D(2),v.dtu,v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);//(v.c,i,D(0),v.dtu,v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);//(v.c,i,D(0),v.dtu,v.dom,v.coefu);
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

  TestFunction<d> gradU; gradU.init(N,1);
  for(int i=0;i<N;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));

    gradU.A(i,0) = new ItemList<d>(d+1);//(v.c,i,D(0),v.dtu,v.dom,v.coefu);
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
TestFunction<d> Eps(const TestFunction<d> & T){
  assert(T.A.N() == d && T.A.M() == 1);

  TestFunction<d> epsU; epsU.init(d, d);
  for(int i=0;i<d;++i) {
    assert(T.A(i,0)->size() == 1);
    const ItemTestFunction<d>& v(T.A(i,0)->getItem(0));
    for(int j=0;j<d;++j) {
      if(i==j){
        epsU.A(i,j) = new ItemList<d>(1);
        ItemTestFunction<d>& u(epsU.A(i,j)->getItem(0));
        u = v;
        // u.cu = i;
        u.du = D(i);
        //(v.c,i,D(i),v.dtu,v.dom,v.coefu);

      }
      else  {
        epsU.A(i,j) = new ItemList<d>(2);
        {
          ItemTestFunction<d>& u(epsU.A(i,j)->getItem(0));
          u=v;
          u.c = 0.5*v.c;
          u.cu = i;
          u.du = D(j);
          // u.dom = v.dom;
          // u.dtu = v.dtu;
          // u.coefu =v.coefu;
        }
        {
          ItemTestFunction<d>& u(epsU.A(i,j)->getItem(1));
          u=v;
          u.c = 0.5*v.c;
          u.cu = j;
          u.du = D(i);
          // u.dom = v.dom;
          // u.dtu = v.dtu;
          // u.coefu =v.coefu;
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

  TestFunction<d> DDU; DDU.init(T.A.N(), T.A.M());
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
TestFunction<d> jump(const TestFunction<d> & T){
  // assert(T.A.M() == 1);
  int N = T.A.N();
  int M = T.A.M();
  TestFunction<d> jumpU; jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    for(int j=0;j<M;++j) {
      int l = T.A(i,j)->size();
      jumpU.A(i,j) = new ItemList<d>(2*l);
      for(int e=0;e<l;++e) {
        const ItemTestFunction<d>& v(T.A(i,j)->getItem(e));
        {
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e));
          u = v;
          u.dom = 0;
        }
        {
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e+1));
          u=v;
          u.c *= -1;
          u.dom = 1;
        }
      }
    }
  }
  return jumpU;
}

template <int d>
TestFunction<d> jump(const TestFunction<d> & U, const TestFunction<d> & V){
  assert(U.A.M() == 1 && U.A.N() == 1);
  assert(V.A.M() == 1 && V.A.N() == 1);
  int N = U.A.N();
  int M = U.A.M();
  TestFunction<d> jumpU; jumpU.init(U.A.N(), U.A.M());
  for(int i=0;i<N;++i) {
    for(int j=0;j<M;++j) {
      assert(U.A(i,j)->size() == V.A(i,j)->size());
      int l = U.A(i,j)->size();
      jumpU.A(i,j) = new ItemList<d>(2*l);
      for(int e=0;e<l;++e) {
        {
          const ItemTestFunction<d>& v(U.A(i,j)->getItem(e));
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e));
          u = v;
          // u.dom = 0;
        }
        {
          const ItemTestFunction<d>& v(V.A(i,j)->getItem(e));
          ItemTestFunction<d>& u(jumpU.A(i,j)->getItem(2*e+1));
          u=v;
          u.c *= -1;
          // u.dom = 1;
        }
      }
    }
  }
  return jumpU;
}


template <int d>
TestFunction<d> average1(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();
  TestFunction<d> jumpU; jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    int l = T.A(i,0)->size();
    jumpU.A(i,0) = new ItemList<d>(2*l);
    for(int e=0;e<l;++e) {
      const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
        u = v;
        // u.c *= kap1;
        u.dom = 0;
        u.addParameter("kappa1");
      }
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
        u=v;
        // u.c *= kap2;
        u.dom = 1;
        u.addParameter("kappa2");
      }
    }
  }
  return jumpU;
}

template <int d>
TestFunction<d> average2(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();
  TestFunction<d> jumpU; jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    int l = T.A(i,0)->size();
    jumpU.A(i,0) = new ItemList<d>(2*l);
    for(int e=0;e<l;++e) {
      const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
        u = v;
        // u.c *= kap2;
        u.dom = 0;
        u.addParameter("kappa2");
      }
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
        u=v;
        // u.c *= kap1;
        u.dom = 1;
        u.addParameter("kappa1");
      }
    }
  }
  return jumpU;
}

template <int d>
TestFunction<d> average(const TestFunction<d> & T){
  assert(T.A.M() == 1);
  int N = T.A.N();
  TestFunction<d> jumpU; jumpU.init(T.A.N(), T.A.M());
  for(int i=0;i<N;++i) {
    int l = T.A(i,0)->size();
    jumpU.A(i,0) = new ItemList<d>(2*l);
    for(int e=0;e<l;++e) {
      const ItemTestFunction<d>& v(T.A(i,0)->getItem(e));
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e));
        u = v;
        u.c *= 0.5;
        u.dom = 0;
      }
      {
        ItemTestFunction<d>& u(jumpU.A(i,0)->getItem(2*e+1));
        u=v;
        u.c *= 0.5;
        u.dom = 1;
      }
    }
  }
  return jumpU;
}


typedef TestFunction<2> TestFunction2;
typedef TestFunction<3> TestFunction3;

#endif
