
#ifndef _TEST_FUNCTION_HPP
#define _TEST_FUNCTION_HPP

#include "../FESpace/expression.hpp"
#include "CutFEM_parameter.hpp"

std::string whichOperator(int op);
std::string whichOperator(int op, int cu);
std::string whichOperatorV(int op, int cu);

static int D(int i) {
   int op[3] = {op_dx, op_dy, op_dz};
   return op[i];
}
static int D2(int i, int j) {
   int op[9] = {op_dxx, op_dxy, op_dxz, op_dyx, op_dyy,
                op_dyz, op_dzx, op_dzy, op_dzz};
   return op[i + 3 * j];
}
static int nextDerivative(int c, int du) {
   assert(du <= op_dz);
   if (du == op_id)
      return D(c);
   else
      return D2(du - 1, c);
}

void f_id(RNMK_ &x, int cu, int du);
void f_ln(RNMK_ &x, int cu, int du);

template <int D> struct TestFunction;

template <int N = 2> struct ItemTestFunction {
   typedef typename typeMesh<N>::Mesh Mesh;
   typedef TestFunction<N> testFun_t;

   double c;
   int cu, du, dtu;
   std::vector<int> ar_nu, conormal;
   std::vector<const VirtualParameter *> coefu;
   int domain_id_, face_side_;
   const ExpressionVirtual *expru = nullptr;
   const GFESpace<Mesh> *fespace  = nullptr;
   const testFun_t *root_fun_p    = nullptr;

   void (*pfun)(RNMK_ &, int, int) = f_id;

   ItemTestFunction()
       : c(0.), cu(-1), du(-1), dtu(-1), domain_id_(-1), face_side_(-1) {}
   ItemTestFunction(const GFESpace<Mesh> *fes, double cc, int i, int j, int tu,
                    int dd, testFun_t *ff)
       : c(cc), cu(i), du(j), dtu(tu), coefu(0), domain_id_(dd), face_side_(-1),
         fespace(fes), root_fun_p(ff){};
   ItemTestFunction(const ItemTestFunction &F)
       : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu),
         conormal(F.conormal), domain_id_(F.domain_id_),
         face_side_(F.face_side_), expru(F.expru), fespace(F.fespace),
         root_fun_p(F.root_fun_p) {
      for (int i = 0; i < F.coefu.size(); ++i)
         coefu.push_back(F.coefu[i]);
   }
   ItemTestFunction(const ItemTestFunction &F, void (*f)(RNMK_ &, int, int))
       : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu),
         conormal(F.conormal), domain_id_(F.domain_id_),
         face_side_(F.face_side_), expru(F.expru), fespace(F.fespace),
         root_fun_p(F.root_fun_p), pfun(f) {
      for (int i = 0; i < F.coefu.size(); ++i)
         coefu.push_back(F.coefu[i]);
   }

   ItemTestFunction &operator=(const ItemTestFunction &L) {
      c          = L.c;
      cu         = L.cu;
      du         = L.du;
      dtu        = L.dtu;
      ar_nu      = L.ar_nu;
      conormal   = L.conormal;
      domain_id_ = L.domain_id_;
      face_side_ = L.face_side_;
      coefu      = L.coefu;
      expru      = L.expru;
      fespace    = L.fespace;
      root_fun_p = L.root_fun_p;
      return *this;
   }
   ItemTestFunction(ItemTestFunction &&v) = default;

   void addNormal(int i) { ar_nu.push_back(i); }
   void addTangent(int i) {
      int Ni = (i == 0); // which normal component
      ar_nu.push_back(Ni);
      if (Ni == 1)
         c *= -1; //(-b,a) with (a,b) normal
   }
   void addCoNormal(int i) { conormal.push_back(i); }
   void addParameter(const VirtualParameter &x) { coefu.push_back(&x); }
   void addParameter(const VirtualParameter *x) { coefu.push_back(x); }

   bool operator==(const ItemTestFunction &F) const {
      if (root_fun_p == F.root_fun_p && coefu.size() == 0 &&
          F.coefu.size() == 0 && cu == F.cu && du == F.du && dtu == F.dtu &&
          face_side_ == F.face_side_ && domain_id_ == F.domain_id_ &&
          ar_nu == F.ar_nu && conormal == F.conormal && expru == F.expru &&
          fespace == F.fespace)
         return true;
      else
         return false;
   }

   friend std::ostream &operator<<(std::ostream &f, const ItemTestFunction &u) {
      std::string n[3] = {"nx", "ny", "nz"};
      f << " FESpace => " << u.fespace << "\t";
      f << std::to_string(u.c) << " * " << whichOperator(u.du, u.cu);
      for (int i = 0; i < u.ar_nu.size(); ++i)
         f << " * " << n[u.ar_nu[i]];

      // for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
      f << "\t in Omega_" << u.domain_id_;
      f << "\n";
      return f;
   }
};

template <int N = 2> struct ItemList {

   typedef typename typeMesh<N>::Mesh Mesh;
   typedef GFESpace<Mesh> FESpace;
   typedef ItemTestFunction<N> item_t;
   typedef TestFunction<N> testFun_t;

   std::vector<item_t> U;

   ItemList() {}
   ItemList(const ItemList &L) {
      for (const auto &item : L.U)
         U.push_back(item);
   }
   ItemList(ItemList &&v) = default;

   // ItemList(const ItemList &L, void (*f)(RNMK_ &, int, int)) {
   //    for (int i = 0; i < U.size(); ++i)
   //       U.push_back(item_t(*L.U(i), f));
   // }

   ItemList(const FESpace &Vh, double cc, int i, int j, int dd, testFun_t *ff) {
      U.push_back(item_t(&Vh, cc, i, j, 0, dd, ff));
   }

   int size() const { return U.size(); }
   const ItemTestFunction<N> &operator()(int i) const { return U[i]; }
   const ItemTestFunction<N> &getItem(int i) const { return U[i]; }
   ItemTestFunction<N> &getItem(int i) { return U[i]; }
   void push(const item_t &u) { U.push_back(u); }

   friend std::ostream &operator<<(std::ostream &f, const ItemList &l) {
      for (int i = 0; i < l.size(); ++i) {
         f << l(i);
      }
      return f;
   }
};

template <int dim = 2> struct TestFunction {
   typedef typename typeRd<dim>::Rd Rd;
   typedef typename typeMesh<dim>::Mesh Mesh;
   typedef GFESpace<Mesh> FESpace;
   typedef ItemTestFunction<dim> item_t;
   typedef ItemList<dim> itemList_t;

   std::vector<std::vector<itemList_t>> A;

   TestFunction() {}
   TestFunction(const FESpace &Vh, int d) {
      for (int i = 0; i < d; ++i)
         this->push(itemList_t(Vh, 1, i, 0, -1, this));
   }
   TestFunction(const FESpace &Vh, int d, int comp0) {
      for (int i = 0; i < d; ++i)
         this->push(itemList_t(Vh, 1, comp0 + i, 0, -1, this));
   }
   TestFunction(const FESpace &Vh, int d, int comp0, int domm) {
      for (int i = 0; i < d; ++i)
         this->push(itemList_t(Vh, 1, comp0 + i, 0, domm, this));
   }
   TestFunction(TestFunction &&U) = default;
   // TestFunction(const TestFunction &U, void (*f)(RNMK_ &, int, int)) {
   //    A.init(U.nbRow(), U.nbCol());
   //    for (int i = 0; i < nbRow(); ++i) {
   //       for (int j = 0; j < nbCol(); ++j) {
   //          A(i, j) = new ItemList<dim>(*U.A(i, j), f);
   //       }
   //    }
   // }

   itemList_t &operator()(int i, int j) { return A[i][j]; }
   const itemList_t &operator()(int i, int j) const { return A[i][j]; }
   const itemList_t &getList(int i, int j) const { return A[i][j]; }
   const item_t &getItem(std::pair<int, int> ij, int k) const {
      return getList(ij.first, ij.second).getItem(k);
   }
   int sizeItemList(int i, int j) const { return A[i][j].size(); }

   void push(const itemList_t &u) {
      std::vector<itemList_t> l_u(1, u);
      A.push_back(l_u);
   }
   void push(const std::vector<itemList_t> &l_u) { A.push_back(l_u); }
   int nbRow() const { return A.size(); }
   int nbCol(int r = 0) const { return A[r].size(); }
   std::pair<int, int> size() const { return {nbRow(), nbCol()}; }
   bool isScalar() const { return (nbRow() == 1 && nbCol() == 1); }

   const std::vector<itemList_t> &getRow(int i) const { return A[i]; }

   void push(std::pair<int, int> ij, const itemList_t &u) {
      auto [nr, nc] = size();
      auto [i, j]   = ij;
      if (nr <= i)
         A.resize(i + 1);
      if (nbCol(i) <= j)
         A[i].resize(j + 1);
      for (const auto &ui : u.U)
         A[i][j].push(ui);
   }

   // TestFunction t() const {
   //    TestFunction Ut(this->nbCol(),
   //                    this->nbRow()); // Ut.init(nbCol(), nbRow());
   //    for (int i = 0; i < this->nbRow(); ++i) {
   //       for (int j = 0; j < this->nbCol(); ++j) {
   //          Ut.A(j, i) = new ItemList<dim>(A(i, j)->size());
   //          for (int ui = 0; ui < A(i, j)->size(); ++ui) {
   //             ItemTestFunction<dim> &v(A(i, j)->getItem(ui));
   //             ItemTestFunction<dim> &u(Ut.A(j, i)->getItem(ui));
   //             u = v;
   //          }
   //       }
   //    }
   //    return Ut;
   // }

   // TestFunction operator*(const Normal &N) {
   //    // assert(!(nbCol() == 1 && nbRow() == dim ));// no column accepted
   //    assert(nbCol() == dim || nbCol() == 1);
   //    bool scalar(nbCol() == 1 && nbRow() == 1);
   //    bool line(nbCol() == dim && nbRow() == 1);
   //    bool column(nbCol() == 1 && nbRow() == dim);

   //    int s = (scalar) ? dim : ((column) ? 1 : nbRow());

   //    TestFunction Un(s);
   //    if (scalar) {
   //       int ksum = 0;
   //       ksum += A(0, 0)->size();
   //       for (int i = 0; i < dim; ++i) {
   //          int k      = 0;
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int ui = 0; ui < A(0, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(0, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //             u = v;
   //             u.addNormal(i);
   //          }
   //       }
   //    } else if (column) {
   //       int ksum = 0, k = 0;
   //       for (int i = 0; i < nbRow(); ++i)
   //          ksum += A(i, 0)->size();
   //       Un.A(0, 0) = new ItemList<dim>(ksum);
   //       for (int i = 0; i < nbRow(); ++i) {
   //          for (int ui = 0; ui < A(i, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(i, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(0, 0)->getItem(k));
   //             u = v;
   //             u.addNormal(i);
   //          }
   //       }
   //    } else {
   //       for (int i = 0; i < nbRow(); ++i) {
   //          int ksum = 0, k = 0;
   //          for (int j = 0; j < nbCol(); ++j)
   //             ksum += A(i, j)->size();
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int j = 0; j < nbCol(); ++j) {
   //             for (int ui = 0; ui < A(i, j)->size(); ++ui, ++k) {
   //                ItemTestFunction<dim> &v(A(i, j)->getItem(ui));
   //                ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //                u = v;
   //                u.addNormal(j);
   //             }
   //          }
   //       }
   //    }
   //    return Un;
   // }
   // TestFunction operator*(const Tangent &N) {
   //    // assert(!(nbCol() == 1 && nbRow() == dim ));// no column accepted
   //    assert(nbCol() == dim || nbCol() == 1);
   //    bool scalar(nbCol() == 1 && nbRow() == 1);
   //    bool line(nbCol() == dim && nbRow() == 1);
   //    bool column(nbCol() == 1 && nbRow() == dim);

   //    int s = (scalar) ? dim : ((column) ? 1 : nbRow());

   //    TestFunction Un(s);
   //    if (scalar) {
   //       int ksum = 0;
   //       ksum += A(0, 0)->size();
   //       for (int i = 0; i < dim; ++i) {
   //          int k      = 0;
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int ui = 0; ui < A(0, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(0, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //             u = v;
   //             u.addTangent(i);
   //          }
   //       }
   //    } else if (column) {
   //       int ksum = 0, k = 0;
   //       for (int i = 0; i < nbRow(); ++i)
   //          ksum += A(i, 0)->size();
   //       Un.A(0, 0) = new ItemList<dim>(ksum);
   //       for (int i = 0; i < nbRow(); ++i) {
   //          for (int ui = 0; ui < A(i, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(i, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(0, 0)->getItem(k));
   //             u = v;
   //             u.addTangent(i);
   //          }
   //       }
   //    } else {
   //       for (int i = 0; i < nbRow(); ++i) {
   //          int ksum = 0, k = 0;
   //          for (int j = 0; j < nbCol(); ++j)
   //             ksum += A(i, j)->size();
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int j = 0; j < nbCol(); ++j) {
   //             for (int ui = 0; ui < A(i, j)->size(); ++ui, ++k) {
   //                ItemTestFunction<dim> &v(A(i, j)->getItem(ui));
   //                ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //                u = v;
   //                u.addTangent(j);
   //             }
   //          }
   //       }
   //    }
   //    return Un;
   // }
   // TestFunction operator*(const Conormal &N) {
   //    // assert(!(nbCol() == 1 && nbRow() == dim ));// no column accepted
   //    assert(nbCol() == dim || nbCol() == 1);
   //    bool scalar(nbCol() == 1 && nbRow() == 1);
   //    bool line(nbCol() == dim && nbRow() == 1);
   //    bool column(nbCol() == 1 && nbRow() == dim);

   //    int s = (scalar) ? dim : ((column) ? 1 : nbRow());

   //    TestFunction Un(s);
   //    if (scalar) {
   //       int ksum = 0;
   //       ksum += A(0, 0)->size();
   //       for (int i = 0; i < dim; ++i) {
   //          int k      = 0;
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int ui = 0; ui < A(0, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(0, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //             u = v;
   //             u.addCoNormal(i);
   //          }
   //       }
   //    } else if (column) {
   //       int ksum = 0, k = 0;
   //       for (int i = 0; i < nbRow(); ++i)
   //          ksum += A(i, 0)->size();
   //       Un.A(0, 0) = new ItemList<dim>(ksum);
   //       for (int i = 0; i < nbRow(); ++i) {
   //          for (int ui = 0; ui < A(i, 0)->size(); ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(i, 0)->getItem(ui));
   //             ItemTestFunction<dim> &u(Un.A(0, 0)->getItem(k));
   //             u = v;
   //             u.addCoNormal(i);
   //          }
   //       }
   //    } else {
   //       for (int i = 0; i < nbRow(); ++i) {
   //          int ksum = 0, k = 0;
   //          for (int j = 0; j < nbCol(); ++j)
   //             ksum += A(i, j)->size();
   //          Un.A(i, 0) = new ItemList<dim>(ksum);
   //          for (int j = 0; j < nbCol(); ++j) {
   //             for (int ui = 0; ui < A(i, j)->size(); ++ui, ++k) {
   //                ItemTestFunction<dim> &v(A(i, j)->getItem(ui));
   //                ItemTestFunction<dim> &u(Un.A(i, 0)->getItem(k));
   //                u = v;
   //                u.addCoNormal(j);
   //             }
   //          }
   //       }
   //    }
   //    return Un;
   // }
   // TestFunction operator*(const Projection &Pg) {
   //    // Only for scalr right now
   //    assert(nbRow() == 1 && nbCol() == 1);
   //    int N = dim;
   //    int s = 1;
   //    TestFunction Un(N, N);
   //    // Un.init(N, N);
   //    const int ksum = A(0, 0)->size();

   //    for (int i = 0; i < N; ++i) {
   //       for (int j = 0; j < N; ++j) {
   //          int k      = 0;
   //          int cst    = 1 + (i == j);
   //          Un.A(i, j) = new ItemList<dim>(cst * ksum);

   //          for (int ui = 0; ui < ksum; ++ui, ++k) {
   //             ItemTestFunction<dim> &v(A(0, 0)->getItem(ui));

   //             if (i == j) {
   //                ItemTestFunction<dim> &u(Un.A(i, j)->getItem(k));
   //                u = v;
   //                k++;
   //             }
   //             ItemTestFunction<dim> &u(Un.A(i, j)->getItem(k));
   //             u = v;
   //             u.c *= -1;
   //             u.addNormal(i);
   //             u.addNormal(j);
   //             k++;
   //          }
   //       }
   //    }
   //    return Un;
   // }

   friend std::ostream &operator<<(std::ostream &f, const TestFunction &u) {
      f << u.nbRow() << " * " << u.nbCol() << std::endl;
      for (int i = 0; i < u.nbRow(); ++i) {
         for (int j = 0; j < u.nbCol(); ++j) {
            f << u(i, j);
         }
      }
      return f;
   }
};

template <int D> void simplify(TestFunction<D> &T) {

   for (auto &row_i : T.A) {
      for (auto &list_item : row_i) {
         for (int i = 0; i < list_item.size(); ++i) {
            for (int j = i + 1; j <= list_item.size(); ++j) {
               if (list_item(i) == list_item(j)) {
                  list_item.U[i].c += list_item.U[j].c;
                  list_item.U[j].c = 0;
               }
            }
         }
      }
   }

   auto is_nul = [](typename TestFunction<D>::item_t item) -> bool {
      return (fabs(item.c) < globalVariable::Epsilon);
   };
   // remove all zero
   for (auto &row_i : T.A) {
      for (auto &list_item : row_i) {
         auto iterator =
             std::remove_if(list_item.U.begin(), list_item.U.end(), is_nul);
         list_item.U.erase(iterator, list_item.U.end());
      }
   }
}
template <int D>
TestFunction<D> operator+(const TestFunction<D> &f1,
                          const TestFunction<D> &f2) {
   assert(f1.size() == f2.size());
   TestFunction<D> sum_fi;
   for (int i = 0; i < f1.nbRow(); ++i) {
      for (int j = 0; j < f1.nbCol(); ++j) {
         const auto &item_list1(f1.getList(i, j));
         sum_fi.push({i, j}, item_list1);
         const auto &item_list2(f2.getList(i, j));
         sum_fi.push({i, j}, item_list2);
      }
   }
   simplify(sum_fi);
   return sum_fi;
}

template <int D>
TestFunction<D> operator-(const TestFunction<D> &f1,
                          const TestFunction<D> &f2) {
   assert(f1.size() == f2.size());
   TestFunction<D> sub_fi;
   for (int i = 0; i < f1.nbRow(); ++i) {
      for (int j = 0; j < f1.nbCol(); ++j) {
         const auto &item_list1(f1.getList(i, j));
         sub_fi.push({i, j}, item_list1);
         auto new_list(f2.getList(i, j));
         for (auto &item : new_list.U)
            item.c *= -1;
         sub_fi.push({i, j}, new_list);
      }
   }
   simplify(sub_fi);
   return sub_fi;
}

template <int D>
TestFunction<D> operator*(const TestFunction<D> &T, double cc) {
   TestFunction<D> mulT;
   for (auto row_i : T.A) {
      for (auto &list_item : row_i) {
         for (auto &item : list_item.U) {
            item.c *= cc;
         }
      }
      mulT.push(row_i);
   }
   return mulT;
}

template <int D>
TestFunction<D> operator*(double cc, const TestFunction<D> &T) {
   return T * cc;
}

template <int D> TestFunction<D> dt(const TestFunction<D> &T) {
   assert(T.nbCol() == 1);

   TestFunction<D> dt_U;
   for (auto comp_i : T.A) {
      for (auto &item_list : comp_i) {
         for (auto &item : item_list.U)
            item.dtu = op_dx;
      }
      dt_U.push(comp_i);
   }
   return dt_U;
}

template <int D> TestFunction<D> dx(const TestFunction<D> &T) {
   assert(T.nbCol() == 1);

   TestFunction<D> dx_U;
   for (auto comp_i : T.A) {
      for (auto &item_list : comp_i) {
         for (auto &item : item_list.U)
            item.du = op_dx;
      }
      dx_U.push(comp_i);
   }
   return dx_U;
}

template <int D> TestFunction<D> dy(const TestFunction<D> &T) {
   assert(T.nbCol() == 1);

   TestFunction<D> dx_U;
   for (auto comp_i : T.A) {
      for (auto &item_list : comp_i) {
         for (auto &item : item_list.U)
            item.du = op_dy;
      }
      dx_U.push(comp_i);
   }
   return dx_U;
}

template <int D> TestFunction<D> dz(const TestFunction<D> &T) {
   assert(T.nbCol() == 1);

   TestFunction<D> dx_U;
   for (auto comp_i : T.A) {
      for (auto &item_list : comp_i) {
         for (auto &item : item_list.U)
            item.du = op_dz;
      }
      dx_U.push(comp_i);
   }
   return dx_U;
}

template <int D> TestFunction<D> grad(const TestFunction<D> &T) {

   auto [N, M] = T.size();
   TestFunction<D> gradU;
   for (int i = 0; i < N; ++i) {
      assert(T.nbCol(i) == 1);
      for (int d = 0; d < D; ++d) {
         auto new_list = T.getList(i, 0);
         for (auto &item : new_list.U) {
            item.du = nextDerivative(d, item.du);
         }
         int irow = T.isScalar() ? d : i;
         int jrow = T.isScalar() ? 0 : d;
         gradU.push({irow, jrow}, new_list);
      }
   }
   return gradU;
}

// template <int N>
// TestFunction<N>
// operator*(std::list<ExpressionFunFEM<typename typeMesh<N>::Mesh>> fh,
//           const TestFunction<N> &F) {
//    assert(F.nbCol() == 1);
//    assert(fh.size() == 1 || fh.size() == N);
//    assert(F.nbRow() == 1 || F.nbRow() == N);
//    int dim_vh = F.nbRow();
//    int dim_fh = fh.size();
//    bool scalar(fh.size() == F.nbRow());
//    int s = (scalar) ? 1 : N;
//    TestFunction<N> Un(s); // initialize a scalar or vector
//    if (scalar) {
//       int ksum = 0;
//       for (int i = 0; i < dim_vh; ++i)
//          ksum += F.A(i, 0)->size();
//       Un.A(0, 0) = new ItemList<N>(ksum);
//       auto it    = fh.begin();
//       if (dim_fh == 1) {
//          int k = 0;
//          for (int i = 0; i < F.nbRow(); ++i) {
//             for (int j = 0; j < F.nbCol(); ++j) {
//                for (int ui = 0; ui < F.A(i, j)->size(); ++ui) {
//                   const ItemTestFunction<N> &v(F.A(i, j)->getItem(ui));
//                   ItemTestFunction<N> &u(Un.A(0, 0)->getItem(k));
//                   u       = v;
//                   u.expru = &(*it);
//                   k++;
//                }
//             }
//          }
//       } else {
//          int k = 0;
//          for (int i = 0; i < F.nbRow(); ++i, ++it) {
//             for (int j = 0; j < F.nbCol(); ++j) {
//                for (int ui = 0; ui < F.A(i, j)->size(); ++ui) {
//                   const ItemTestFunction<N> &v(F.A(i, j)->getItem(ui));
//                   ItemTestFunction<N> &u(Un.A(0, 0)->getItem(k));
//                   u       = v;
//                   u.expru = &(*it);
//                   k++;
//                }
//             }
//          }
//       }
//    } else {              // result is a vector
//       if (dim_fh == 1) { // v_h is a vector , fh is a scalar
//          auto it = fh.begin();
//          for (int i = 0; i < F.nbRow(); ++i) {
//             int k      = 0;
//             int ksum   = F.A(i, 0)->size();
//             Un.A(i, 0) = new ItemList<N>(ksum);
//             for (int ui = 0; ui < ksum; ++ui) {
//                const ItemTestFunction<N> &v(F.A(i, 0)->getItem(ui));
//                ItemTestFunction<N> &u(Un.A(i, 0)->getItem(k));
//                u       = v;
//                u.expru = &(*it);
//                k++;
//             }
//          }
//       } else { // v_h is a scalar , fh is a vector
//          int ksum = F.A(0, 0)->size();
//          auto it  = fh.begin();
//          for (int i = 0; i < dim_fh; ++i, ++it) {
//             int k      = 0;
//             Un.A(i, 0) = new ItemList<N>(ksum);
//             for (int ui = 0; ui < ksum; ++ui) {
//                const ItemTestFunction<N> &v(F.A(0, 0)->getItem(ui));
//                ItemTestFunction<N> &u(Un.A(i, 0)->getItem(k));
//                u       = v;
//                u.expru = &(*it);
//                k++;
//             }
//          }
//       }
//    }
//    // int k = 0,ksum=0;
//    // for(int i=0;i<N;++i) ksum += F.A(i,0)->size();
//    // multU.A(0,0) = new ItemList<N>(ksum);
//    //
//    // auto it = fh.begin();
//    // for(int i=0;i<F.nbRow();++i, ++it) {
//    //   for(int j=0;j<F.nbCol();++j) {
//    //     for(int ui=0;ui<F.A(i,j)->size();++ui) {
//    //       const ItemTestFunction<N>& v(F.A(i,j)->getItem(ui));
//    //       ItemTestFunction<N>& u(multU.A(0,0)->getItem(k));
//    //       u = v;
//    //       u.expru = &(*it);
//    //       k++;
//    //     }
//    //   }
//    // }
//    return Un;
// }

// template <int N>
// TestFunction<N> operator*(const ExpressionVirtual &expr,
//                           const TestFunction<N> &F) {
//    TestFunction<N> multU(F);
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < F.A(i, j)->size(); ++ui) {
//             ItemTestFunction<N> &v(multU.A(i, j)->getItem(ui));
//             v.expru = &expr;
//          }
//       }
//    }
//    return multU;
// }

// template <int N>
// TestFunction<N> operator*(const TestFunction<N> &F,
//                           const ExpressionVirtual &expr) {
//    TestFunction<N> multU(F);
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < F.A(i, j)->size(); ++ui) {
//             ItemTestFunction<N> &v(multU.A(i, j)->getItem(ui));
//             v.expru = &expr;
//          }
//       }
//    }
//    return multU;
// }

// template <int N>
// TestFunction<N> operator*(const TestFunction<N> &F, const
// Normal_Component &c) {
//    TestFunction<N> multU(F);
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          *multU.A(i, j) *= c;
//       }
//    }
//    return multU;
// }

// template <int N>
// TestFunction<N> operator*(const TestFunction<N> &F,
//                           const VirtualParameter &cc) {
//    TestFunction<N> multU(F);
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < multU.A(i, j)->size(); ++ui) {
//             ItemTestFunction<N> &v(multU.A(i, j)->getItem(ui));
//             v.addParameter(cc);
//          }
//       }
//    }
//    return multU;
// }

// template <int N>
// TestFunction<N> operator*(const VirtualParameter &cc,
//                           const TestFunction<N> &F) {
//    TestFunction<N> multU(F);
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < multU.A(i, j)->size(); ++ui) {
//             ItemTestFunction<N> &v(multU.A(i, j)->getItem(ui));
//             v.addParameter(cc);
//          }
//       }
//    }
//    return multU;
// }

// template <int d>
// TestFunction<d> operator*(const CutFEM_Rd<d> &cc, const TestFunction<d>
// &T) {
//    assert(T.nbCol() == 1);
//    int N = T.nbRow();

//    bool scalar = (N == 1);
//    int r       = (scalar) ? d : 1;
//    TestFunction<d> resU(r, 1);
//    if (scalar) {
//       int nitem = T.A(0, 0)->size();

//       for (int j = 0; j < d; ++j) {
//          resU.A(j, 0) = new ItemList<d>(nitem);
//          for (int ui = 0; ui < nitem; ++ui) {
//             const ItemTestFunction<d> &v(T.A(0, 0)->getItem(ui));
//             ItemTestFunction<d> &u(resU.A(j, 0)->getItem(ui));
//             u = v;
//             u.addParameter(cc.get_parameter(j));
//          }
//       }
//    } else {
//       assert(N == d); // column
//       int nitem = 0;
//       for (int i = 0; i < d; ++i) {
//          nitem += T.A(i, 0)->size();
//       }
//       resU.A(0, 0) = new ItemList<d>(nitem);
//       int k        = 0;
//       for (int j = 0; j < d; ++j) {
//          int nloc = T.A(j, 0)->size();
//          for (int ui = 0; ui < nloc; ++ui) {
//             const ItemTestFunction<d> &v(T.A(j, 0)->getItem(ui));
//             ItemTestFunction<d> &u(resU.A(0, 0)->getItem(k++));
//             u = v;
//             u.addParameter(cc.get_parameter(j));
//          }
//       }
//    }
//    return resU;
// }

// template <int N> TestFunction<N> ln(const TestFunction<N> &F) {
//    return TestFunction<N>(F, f_ln);
// }

// template <int d>
// TestFunction<d> operator*(const typename typeRd<d>::Rd &cc,
//                           const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);
//    int N       = T.nbRow();
//    bool scalar = (N == 1);
//    int r       = (scalar) ? d : 1; // dim of the result
//    TestFunction<d> resU(r, 1);
//    if (scalar) {
//       int nitem = T.A(0, 0)->size();
//       for (int j = 0; j < d; ++j) {
//          resU.A(j, 0) = new ItemList<d>(nitem);
//          for (int ui = 0; ui < nitem; ++ui) {
//             const ItemTestFunction<d> &v(T.A(0, 0)->getItem(ui));
//             ItemTestFunction<d> &u(resU.A(j, 0)->getItem(ui));
//             u = v;
//             u.c *= cc[j];
//          }
//       }
//    } else {
//       assert(N == d); // column
//       int nitem = 0;
//       for (int i = 0; i < d; ++i) {
//          nitem += T.A(i, 0)->size();
//       }
//       resU.A(0, 0) = new ItemList<d>(nitem);
//       int k        = 0;
//       for (int j = 0; j < d; ++j) {
//          int nloc = T.A(j, 0)->size();
//          for (int ui = 0; ui < nloc; ++ui) {
//             const ItemTestFunction<d> &v(T.A(j, 0)->getItem(ui));
//             ItemTestFunction<d> &u(resU.A(0, 0)->getItem(k++));
//             u = v;
//             u.c *= cc[j];
//          }
//       }
//    }
//    return resU;
// }

// template <int d> TestFunction<d> gradS(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);

//    int N       = T.nbRow();
//    bool scalar = (N == 1);
//    int col     = (scalar) ? 1 : d;
//    TestFunction<d> gradU(d, col);

//    if (scalar) {
//       assert(T.A(0, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(0, 0)->getItem(0));

//       for (int i = 0; i < d; ++i) {
//          gradU.A(i, 0) = new ItemList<d>(d + 1);

//          int kk = 0;
//          for (int j = 0; j < d; ++j) {
//             if (i == j) {
//                ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//                u    = v;
//                // u.c = 1;
//                // u.cu = i;
//                u.du = D(i);
//             }
//             ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//             u = v;
//             u.c *= -1;
//             // u.cu = i;
//             u.du = D(j);
//             u.addNormal(i);
//             u.addNormal(j);
//          }
//       }
//    } else {
//       int n = T.A(0, 0)->size();

//       for (int row = 0; row < d; ++row) {
//          for (int col = 0; col < d; ++col) {
//             gradU.A(row, col) = new ItemList<d>(n * (d + 1));
//          }
//       }

//       for (int row = 0; row < d; ++row) {
//          assert(T.A(row, 0)->size() == n);

//          // const ItemTestFunction<d>& v(T.A(row,0)->getItem(0));

//          for (int col = 0; col < d; ++col) {
//             int kk = 0;

//             for (int l = 0; l < n; ++l) {
//                const ItemTestFunction<d> &v(T.A(row, 0)->getItem(l));

//                for (int j = 0; j < d; ++j) {
//                   if (col == j) {
//                      ItemTestFunction<d> &u(gradU.A(row,
//                      col)->getItem(kk++)); u    = v; u.du = D(col);
//                   }
//                   ItemTestFunction<d> &u(gradU.A(row,
//                   col)->getItem(kk++)); u = v; u.c *= -1; u.du = D(j);
//                   u.addNormal(col);
//                   u.addNormal(j);
//                }
//             }
//          }
//       }
//    }

//    return gradU;
// }

// template <int d> TestFunction<d> div(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);
//    assert(T.nbRow() == d);
//    int N = T.nbRow();

//    TestFunction<d> divU(1);
//    int k = 0, ksum = 0;
//    for (int i = 0; i < N; ++i)
//       ksum += T.A(i, 0)->size();
//    divU.A(0, 0) = new ItemList<d>(ksum);

//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));
//       ItemTestFunction<d> &u(divU.A(0, 0)->getItem(k));
//       u    = v;
//       u.du = D(i);
//       k++;
//    }
//    return divU;
// }

// template <int d> TestFunction<d> divS(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);
//    assert(T.nbRow() == d);
//    int N = T.nbRow();

//    TestFunction<d> divU(1);
//    int k = 0, ksum = 0;
//    for (int i = 0; i < N; ++i)
//       ksum += T.A(i, 0)->size();
//    divU.A(0, 0) = new ItemList<d>(ksum * (d + 1));

//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));

//       ItemTestFunction<d> &u(divU.A(0, 0)->getItem(k));
//       u    = v;
//       u.du = D(i);
//       k++;

//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(divU.A(0, 0)->getItem(k));
//          u = v;
//          u.c *= -1;
//          u.du = D(j);
//          u.addNormal(i);
//          u.addNormal(j);
//          k++;
//       }
//    }
//    return divU;
// }

// template <int d> TestFunction<d> divT(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);
//    assert(T.nbRow() == d);
//    assert(d == 2);
//    int N = T.nbRow();

//    TestFunction<d> divU(1);
//    int k = 0, ksum = 0;
//    for (int i = 0; i < N; ++i)
//       ksum += T.A(i, 0)->size();
//    divU.A(0, 0) = new ItemList<d>(ksum * (d));

//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));

//       for (int j = 0; j < d; ++j) {
//          int ii = (i == 0);
//          int jj = (j == 0);
//          ItemTestFunction<d> &u(divU.A(0, 0)->getItem(k));
//          u = v;
//          u.c *= (ii == jj) ? 1 : -1;
//          u.du = D(j);
//          u.addNormal(jj);
//          u.addNormal(ii);
//          k++;
//       }
//    }
//    return divU;
// }

// template <int d> TestFunction<d> dxS(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1 && T.nbRow() == 1);
//    int N = T.nbRow();

//    TestFunction<d> gradU(N, 1); // gradU.init(N,1);
//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));

//       gradU.A(i, 0) = new ItemList<d>(d + 1);
//       {
//          ItemTestFunction<d> &u1(gradU.A(i, 0)->getItem(0));
//          u1    = v;
//          u1.du = D(0);
//       }
//       int kk = 1;
//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//          u = v;
//          u.c *= -1;
//          u.du = D(j);
//          u.addNormal(0);
//          u.addNormal(j);
//       }
//    }
//    return gradU;
// }

// template <int d> TestFunction<d> dyS(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1 && T.nbRow() == 1);
//    int N = T.nbRow();

//    TestFunction<d> gradU(N, 1); // gradU.init(N,1);
//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));

//       gradU.A(i, 0) = new ItemList<d>(d + 1);
//       {
//          ItemTestFunction<d> &u1(gradU.A(i, 0)->getItem(0));
//          u1    = v;
//          u1.du = D(1);
//       }
//       int kk = 1;
//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//          u = v;
//          u.c *= -1;
//          u.du = D(j);
//          u.addNormal(1);
//          u.addNormal(j);
//       }
//    }
//    return gradU;
// }

// template <int d> TestFunction<d> dzS(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1 && T.nbRow() == 1);
//    int N = T.nbRow();

//    TestFunction<d> gradU(N, 1); // gradU.init(N,1);
//    for (int i = 0; i < N; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));

//       gradU.A(i, 0) = new ItemList<d>(d + 1);
//       {
//          ItemTestFunction<d> &u1(gradU.A(i, 0)->getItem(0));
//          u1    = v;
//          u1.du = D(2);
//       }
//       int kk = 1;
//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//          u = v;
//          u.c *= -1;
//          u.du = D(j);
//          u.addNormal(2);
//          u.addNormal(j);
//       }
//    }
//    return gradU;
// }

// template <int d> TestFunction<d> Eps(const TestFunction<d> &T) {
//    assert(T.nbRow() == d && T.nbCol() == 1);

//    TestFunction<d> epsU(d, d); // epsU.init(d, d);
//    for (int i = 0; i < d; ++i) {
//       assert(T.A(i, 0)->size() == 1);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));
//       for (int j = 0; j < d; ++j) {
//          if (i == j) {
//             epsU.A(i, j) = new ItemList<d>(1);
//             ItemTestFunction<d> &u(epsU.A(i, j)->getItem(0));
//             u    = v;
//             u.du = D(i);

//          } else {
//             epsU.A(i, j) = new ItemList<d>(2);
//             {
//                ItemTestFunction<d> &u(epsU.A(i, j)->getItem(0));
//                u    = v;
//                u.c  = 0.5 * v.c;
//                u.cu = i;
//                u.du = D(j);
//             }
//             {
//                ItemTestFunction<d> &u(epsU.A(i, j)->getItem(1));
//                u    = v;
//                u.c  = 0.5 * v.c;
//                u.cu = j;
//                u.du = D(i);
//             }
//          }
//       }
//    }
//    return epsU;
// }
// //
// //
// template <int d> TestFunction<d> grad2(const TestFunction<d> &T) {
//    assert(T.nbCol() == 1);
//    assert(T.nbRow() == d || T.nbRow() == 1);
//    int N = T.nbRow();

//    TestFunction<d> DDU(T.nbRow(),
//                        T.nbCol()); // DDU.init(T.nbRow(), T.nbCol());
//    for (int i = 0; i < N; ++i)
//       assert(T.A(i, 0)->size() == 1);

//    for (int i = 0; i < N; ++i) {
//       DDU.A(i, 0) = new ItemList<d>(3);
//       const ItemTestFunction<d> &v(T.A(i, 0)->getItem(0));
//       int k = 0;
//       for (int i1 = 0; i1 < d; ++i1) {
//          for (int i2 = i1; i2 < d; ++i2) {
//             // 0,0-> dxx, 0,1->2dxy  1,1->dyy
//             ItemTestFunction<d> &u(DDU.A(i, 0)->getItem(k));
//             u    = v;
//             u.du = D2(i1, i2);
//             u.c  = (i1 != i2) + 1;
//             k++;
//          }
//       }
//    }
//    return DDU;
// }

// template <int d>
// TestFunction<d> jump(const TestFunction<d> &T, int c1, int c2) {
//    // assert(T.nbCol() == 1);
//    int N = T.nbRow();
//    int M = T.nbCol();
//    TestFunction<d> jumpU(T.nbRow(),
//                          T.nbCol()); // jumpU.init(T.nbRow(),
//                          T.nbCol());
//    for (int i = 0; i < N; ++i) {
//       for (int j = 0; j < M; ++j) {
//          int l         = T.A(i, j)->size();
//          jumpU.A(i, j) = new ItemList<d>(2 * l);
//          for (int e = 0; e < l; ++e) {
//             const ItemTestFunction<d> &v(T.A(i, j)->getItem(e));
//             {
//                ItemTestFunction<d> &u(jumpU.A(i, j)->getItem(2 * e));
//                u = v;
//                u.c *= c1;
//                u.face_side_ = 0;
//             }
//             {
//                ItemTestFunction<d> &u(jumpU.A(i, j)->getItem(2 * e +
//                1)); u = v; u.c *= c2; u.face_side_ = 1;
//             }
//          }
//       }
//    }
//    return jumpU;
// }
// template <int d> TestFunction<d> jump(const TestFunction<d> &T) {
//    return jump(T, 1, -1);
// }

// // NEED DO FIXE THIS FUNCTION
// // HAS TO WORK FOR GENERAL DOMAIN
// template <int d>
// TestFunction<d> jump(const TestFunction<d> &U, const TestFunction<d>
// &V, int c1,
//                      int c2) {
//    assert(U.nbCol() == 1 && U.nbRow() == 1);
//    assert(V.nbCol() == 1 && V.nbRow() == 1);
//    int N = U.nbRow();
//    int M = U.nbCol();
//    TestFunction<d> jumpU(U.nbRow(),
//                          U.nbCol()); // jumpU.init(U.nbRow(),
//                          U.nbCol());
//    for (int i = 0; i < N; ++i) {
//       for (int j = 0; j < M; ++j) {
//          assert(U.A(i, j)->size() == V.A(i, j)->size());
//          int l         = U.A(i, j)->size();
//          jumpU.A(i, j) = new ItemList<d>(2 * l);
//          for (int e = 0; e < l; ++e) {
//             {
//                const ItemTestFunction<d> &v(U.A(i, j)->getItem(e));
//                ItemTestFunction<d> &u(jumpU.A(i, j)->getItem(2 * e));
//                u = v;
//                u.c *= c1;
//                u.face_side_ = 0;
//             }
//             {
//                const ItemTestFunction<d> &v(V.A(i, j)->getItem(e));
//                ItemTestFunction<d> &u(jumpU.A(i, j)->getItem(2 * e +
//                1)); u = v; u.c *= c2; u.face_side_ = 1;
//             }
//          }
//       }
//    }
//    return jumpU;
// }

// template <int d>
// TestFunction<d> jump(const TestFunction<d> &U, const TestFunction<d>
// &V) {
//    return jump(U, V, 1, -1);
// }

// template <int d>
// TestFunction<d> average(const TestFunction<d> &T, double v1 = 0.5,
//                         double v2 = 0.5) {
//    assert(T.nbCol() == 1);
//    int N = T.nbRow();
//    TestFunction<d> jumpU(T.nbRow(),
//                          T.nbCol()); // jumpU.init(T.nbRow(),
//                          T.nbCol());
//    for (int i = 0; i < N; ++i) {
//       int l         = T.A(i, 0)->size();
//       jumpU.A(i, 0) = new ItemList<d>(2 * l);
//       for (int e = 0; e < l; ++e) {
//          const ItemTestFunction<d> &v(T.A(i, 0)->getItem(e));
//          {
//             ItemTestFunction<d> &u(jumpU.A(i, 0)->getItem(2 * e));
//             u = v;
//             u.c *= v1;
//             u.face_side_ = 0;
//          }
//          {
//             ItemTestFunction<d> &u(jumpU.A(i, 0)->getItem(2 * e + 1));
//             u = v;
//             u.c *= v2;
//             u.face_side_ = 1;
//          }
//       }
//    }
//    return jumpU;
// }

// template <int d>
// TestFunction<d> average(const TestFunction<d> &T, const
// VirtualParameter &para1,
//                         const VirtualParameter &para2) {
//    assert(T.nbCol() == 1);
//    int N = T.nbRow();
//    TestFunction<d> jumpU(T.nbRow(),
//                          T.nbCol()); // jumpU.init(T.nbRow(),
//                          T.nbCol());
//    for (int i = 0; i < N; ++i) {
//       int l         = T.A(i, 0)->size();
//       jumpU.A(i, 0) = new ItemList<d>(2 * l);
//       for (int e = 0; e < l; ++e) {
//          const ItemTestFunction<d> &v(T.A(i, 0)->getItem(e));
//          {
//             ItemTestFunction<d> &u(jumpU.A(i, 0)->getItem(2 * e));
//             u            = v;
//             u.face_side_ = 0;
//             u.addParameter(para1);
//          }
//          {
//             ItemTestFunction<d> &u(jumpU.A(i, 0)->getItem(2 * e + 1));
//             u            = v;
//             u.face_side_ = 1;
//             u.addParameter(para2);
//          }
//       }
//    }
//    return jumpU;
// }

// template <int d>
// TestFunction<d> average(const TestFunction<d> &T,
//                         const VirtualParameter &para1) {
//    return average(T, para1, para1);
// }

typedef TestFunction<2> TestFunction2;
typedef TestFunction<3> TestFunction3;

#endif
