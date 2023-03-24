/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#ifndef _TEST_FUNCTION_HPP
#define _TEST_FUNCTION_HPP

#include "../FESpace/expression.hpp"
#include "CutFEM_parameter.hpp"

std::string whichOperator(int op);
std::string whichOperator(int op, int cu);
std::string whichOperatorV(int op, int cu);

static int D1(int i) {
    int op[3] = {op_dx, op_dy, op_dz};
    return op[i];
}
static int D2(int i, int j) {
    int op[9] = {op_dxx, op_dxy, op_dxz, op_dyx, op_dyy, op_dyz, op_dzx, op_dzy, op_dzz};
    return op[i + 3 * j];
}
static int nextDerivative(int c, int du) {
    assert(du <= op_dz);
    // c = index we want to differentiate with
    // du = current differentiation index of test function

    if (du == op_id)    
        return D1(c);   // if test function is not differentiated, return first derivative in direction c
    else
        return D2(du - 1, c);   // else, return second derivative in direction c
}

static int rotgradD(int i) {
    int op[2] = {op_dy, op_dx}; 
    return op[i];
}


// Given test function u with u=(ux, uy, uz) where c denotes the component of u, and 

void f_id(RNMK_ &x, int cu, int du);
void f_ln(RNMK_ &x, int cu, int du);

template <int D> struct TestFunction;

/**
 * @brief ItemTestFunction is the core class representing a test function
 * @tparam N Physical dimension //?
 * @note This class holds all necessary information about the test function, 
 * such as the coefficient, the component of the field, the derivative index,
 * the time derivative index, the normal and conormal components, the domain id
 * and the face side. It also holds the pointer to the root TestFunction object
 * and the pointer to the function that evaluates the test function.
 */
template <int N = 2> struct ItemTestFunction {
    typedef typename MeshType<N>::Mesh Mesh;
    typedef TestFunction<N> testFun_t; // why does ItemTestFunction have a TestFunction obejct, when TestFunction has an
                                       // ItemTestFunction object?

    double c; ///< c – scalar multiple of the test function
    int cu;   ///< cu – component of u
    int du;   ///< du – derivative index of u, i.e. dx, dy or dz
    int dtu;  ///< dtu – time derivative index

    std::vector<int> ar_nu, conormal; ///< arrays of normal and conormal components of test function //?
    std::vector<const VirtualParameter *> coefu;    
    int domain_id_, face_side_; ///< e.g. average(u) saves the same test function but with different domain_id/face_side_
    std::shared_ptr<ExpressionVirtual> expru = nullptr;     ///> if multiplied by a function
    const GFESpace<Mesh> *fespace            = nullptr;     ///< so that integration is performed on the right space
    const testFun_t *root_fun_p              = nullptr;

    void (*pfun)(RNMK_ &, int, int) = f_id;

    ItemTestFunction() : c(0.), cu(-1), du(-1), dtu(-1), domain_id_(-1), face_side_(-1) {}
    ItemTestFunction(const GFESpace<Mesh> *fes, double cc, int i, int j, int tu, int dd, testFun_t *ff)
        : c(cc), cu(i), du(j), dtu(tu), coefu(0), domain_id_(dd), face_side_(-1), fespace(fes), root_fun_p(ff){};
    ItemTestFunction(const ItemTestFunction &F)
        : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu), conormal(F.conormal), domain_id_(F.domain_id_),
          face_side_(F.face_side_), expru(F.expru), fespace(F.fespace), root_fun_p(F.root_fun_p) {
        for (int i = 0; i < F.coefu.size(); ++i)
            coefu.push_back(F.coefu[i]);
    }
    ItemTestFunction(const ItemTestFunction &F, void (*f)(RNMK_ &, int, int))
        : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu), conormal(F.conormal), domain_id_(F.domain_id_),
          face_side_(F.face_side_), expru(F.expru), fespace(F.fespace), root_fun_p(F.root_fun_p), pfun(f) {
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
    void addParameter(const VirtualParameter &x) { coefu.push_back(&x); }
    void addParameter(const VirtualParameter *x) { coefu.push_back(x); }

    bool operator==(const ItemTestFunction &F) const {
        if (root_fun_p == F.root_fun_p && coefu.size() == 0 && F.coefu.size() == 0 && cu == F.cu && du == F.du &&
            dtu == F.dtu && face_side_ == F.face_side_ && domain_id_ == F.domain_id_ && ar_nu == F.ar_nu &&
            conormal == F.conormal && expru.get() == F.expru.get() && fespace == F.fespace)
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

/**
 * @brief ItemList handles test functions related to each other through an operator
 * @note For example, if T1 and T2 are two test functions, then the test function T1+T2 
 * has an ItemList of size two, containing both T1 and T2.
 * @tparam N Physical dimension
 */
template <int N = 2> struct ItemList {

    typedef typename MeshType<N>::Mesh Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef ItemTestFunction<N> item_t;
    typedef TestFunction<N> testFun_t;

    std::vector<item_t> U;      ///< vector of Items

    //! When does the vector consist of more than one item?

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

    //! CONTINUE HERE
    /**
     * @brief Construct a new ItemList object
     * 
     * @param Vh Finite Element Space
     * @param cc Scalar coefficient that is multiplied with the testfunction
     * @param i  Component index
     * @param j  Derivative index (op_id = 0, op_dx = 1, op_dy = 2, op_dz = 3 etc)
     * @param dd Domain index: -1 equals main domain if only one is defined
     * @param ff Pointer to the testfunction
     */
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

/**
 * @brief Test function class
 * @note A test function object is just a matrix A of ItemLists. 
 * A[i] returns a vector of the ith components of the test function
 * If one takes grad in 2D of a a 2D vector, it becomes a 2x2 matrix
 * | dx(ux) |  dy(ux) |
 * | dx(uy) |  dy(uy) | 
 * @tparam dim Physical dimension of the problem
 */
template <int dim = 2> struct TestFunction {
    typedef typename typeRd<dim>::Rd Rd;
    typedef typename MeshType<dim>::Mesh Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef ItemTestFunction<dim> item_t;
    typedef ItemList<dim> itemList_t;

    std::vector<std::vector<itemList_t>> A;

    TestFunction() {} ///< Default constructor

    /**
     * @brief Construct a new Test Function object
     *
     * @param Vh Finite Element space
     * @param d Number of components of test function
     */
    TestFunction(const FESpace &Vh, int d) {
        for (int i = 0; i < d; ++i)
            // (FE space, coefficient 1, component i, derivative index identity, domain -1, pointer to testfunction)
            this->push(itemList_t(Vh, 1, i, 0, -1, this));  
    }

    /**
     * @brief Construct a new Test Function object
     *
     * @param Vh Finite Element space
     * @param d Number of total components of test function
     * @param comp0 Component number of first component //?
     */
    TestFunction(const FESpace &Vh, int d, int comp0) {
        for (int i = 0; i < d; ++i)
            this->push(itemList_t(Vh, 1, comp0 + i, 0, -1, this));
    }

    /**
     * @brief Construct a new Test Function object
     *
     * @param Vh Finite Element space
     * @param d Number of total components of test function
     * @param comp0 Component number of first component //?
     * @param domm Domain number
     */
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
    const item_t &getItem(std::pair<int, int> ij, int k) const { return getList(ij.first, ij.second).getItem(k); }
    int sizeItemList(int i, int j) const { return A[i][j].size(); }

    // pushes a new row of ItemList
    void push(const itemList_t &u) {
        std::vector<itemList_t> l_u(1, u);
        A.push_back(l_u);
    }
    void push(const std::vector<itemList_t> &l_u) { A.push_back(l_u); } // push a row of itemList_t 
    int nbRow() const { return A.size(); }      
    int nbCol(int r = 0) const { return A[r].size(); }
    std::pair<int, int> size() const { return {nbRow(), nbCol()}; } // return the size of the matrix A
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

    friend std::ostream &operator<<(std::ostream &f, const TestFunction &u) {
        f << "Number of rows (components): " << u.nbRow() << ", Number of columns: " << u.nbCol() << std::endl;
        for (int i = 0; i < u.nbRow(); ++i) {
            for (int j = 0; j < u.nbCol(); ++j) {
                f << "( " << i << " , " << j << " ) : \t" << u(i, j);
            }
        }
        return f;
    }
};

template <int D> void simplify(TestFunction<D> &T) {

    for (auto &row_i : T.A) {
        for (auto &list_item : row_i) {
            for (int i = 0; i < list_item.size(); ++i) {
                for (int j = i + 1; j < list_item.size(); ++j) {
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
            auto iterator = std::remove_if(list_item.U.begin(), list_item.U.end(), is_nul);
            list_item.U.erase(iterator, list_item.U.end());
        }
    }
}

template <int D> TestFunction<D> transpose(TestFunction<D> &T) {
    TestFunction<D> Ut;
    for (int i = 0; i < T.nbRow(); ++i) {
        for (int j = 0; j < T.nbCol(); ++j) {
            auto list = T.getList(i, j);
            Ut.push({j, i}, list);
        }
    }
    return Ut;
}

template <int D> TestFunction<D> operator+(const TestFunction<D> &f1, const TestFunction<D> &f2) {
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

template <int D> TestFunction<D> operator-(const TestFunction<D> &f1, const TestFunction<D> &f2) {
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

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, double cc) {
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

template <int D> TestFunction<D> operator*(double cc, const TestFunction<D> &T) { return T * cc; }

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, const BaseVector &vec) {
    auto [N, M] = T.size();

    TestFunction<D> Un;

    if (T.isScalar()) {
        for (int i = 0; i < D; ++i) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.addNormal(vec[i]);
                item.c *= vec.cst(i);
            }
            Un.push(new_list);
        }
    } else if (N == D && M == 1) {
        for (int i = 0; i < D; ++i) {
            auto new_list = T.getList(i, 0);
            for (auto &item : new_list.U) {
                item.addNormal(vec[i]);
                item.c *= vec.cst(i);
            }
            Un.push({0, 0}, new_list);
        }
    } else if (N == D && M == D) {
        for (int i = 0; i < T.nbRow(); ++i) {
            for (int j = 0; j < T.nbCol(); ++j) {
                auto new_list = T.getList(i, j);
                for (auto &item : new_list.U) {
                    item.addNormal(vec[j]);
                    item.c *= vec.cst(j);
                }
                Un.push({i, 0}, new_list);
            }
        }
    } else {
        std::cout << " Cannot operate the multiplication witht the normal vector" << std::endl;
        assert(0);
    }
    return Un;
}

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, const Conormal &vec) {
    auto [N, M] = T.size();

    TestFunction<D> Un;

    if (T.isScalar()) {
        for (int i = 0; i < D; ++i) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.addCoNormal(vec[i]);
                item.c *= vec.cst(i);
            }
            Un.push(new_list);
        }
    } else if (N == D && M == 1) {
        for (int i = 0; i < D; ++i) {
            auto new_list = T.getList(i, 0);
            for (auto &item : new_list.U) {
                item.addCoNormal(vec[i]);
                item.c *= vec.cst(i);
            }
            Un.push({0, 0}, new_list);
        }
    } else if (N == D && M == D) {
        for (int i = 0; i < T.nbRow(); ++i) {
            for (int j = 0; j < T.nbCol(); ++j) {
                auto new_list = T.getList(i, j);
                for (auto &item : new_list.U) {
                    item.addCoNormal(vec[j]);
                    item.c *= vec.cst(j);
                }
                Un.push({i, 0}, new_list);
            }
        }
    } else {
        std::cout << " Cannot operate the multiplication witht the normal vector" << std::endl;
        assert(0);
    }
    return Un;
}

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, const Normal_Component &cc) {
    TestFunction<D> multU;
    auto [N, M] = T.size();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            auto new_list = T.getList(i, j);
            for (auto &item : new_list) {
                item.addNormal(cc.component());
            }
            multU.push({i, j}, new_list);
        }
    }
    return multU;
}

template <int D, typename Expr>
TestFunction<D> operator*(const std::list<std::shared_ptr<Expr>> &fh, const TestFunction<D> &T) {
    auto [N, M] = T.size();
    assert(M == 1);

    TestFunction<D> Un;
    if (N == fh.size()) {
        auto it = fh.begin();
        for (int i = 0; i < N; i++, it++) {
            auto new_list = T.getList(i, 0);
            for (auto &item : new_list.U) {
                // if (item.expru.get() == nullptr)
                item.expru = (*it);
                // else {
                //    auto temp_p = item.expru;
                //    item.expru  = temp_p * (*it);
                // }
            }
            Un.push({0, 0}, new_list);
        }
    } else if (N == D && fh.size() == 1) {
        auto it = fh.begin();
        for (int i = 0; i < N; i++) {
            auto new_list = T.getList(i, 0);
            for (auto &item : new_list.U) {
                // if (item.expru.get() == nullptr)
                item.expru = (*it);
                // else {
                //    auto temp_p = item.expru;
                //    item.expru  = temp_p * (*it);
                // }
            }
            Un.push({i, 0}, new_list);
        }
    } else if (N == 1 && fh.size() == D) {
        for (const auto &ff : fh) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                // if (item.expru.get() == nullptr)
                item.expru = ff;
                // else {
                //    auto temp_p = item.expru;
                //    item.expru  = temp_p * ff;
                // }
            }
            Un.push(new_list);
        }
    } else {
        assert(0);
    }
    return Un;
}

template <int D>
TestFunction<D> operator*(const TestFunction<D> &T, const std::list<std::shared_ptr<ExpressionVirtual>> &fh) {
    return fh * T;
}

template <int D> TestFunction<D> operator*(ExpressionFunFEM<typename MeshType<D>::Mesh> &fh, const TestFunction<D> &T) {
    auto fh_p = std::make_shared<ExpressionFunFEM<typename MeshType<D>::Mesh>>(fh);
    std::list<std::shared_ptr<ExpressionFunFEM<typename MeshType<D>::Mesh>>> ff = {fh_p};
    return ff * T;
}

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, ExpressionFunFEM<typename MeshType<D>::Mesh> &fh) {
    return fh * T;
}

template <int D> TestFunction<D> operator*(const std::shared_ptr<ExpressionVirtual> &fh, const TestFunction<D> &T) {
    std::list<std::shared_ptr<ExpressionVirtual>> ff = {fh};
    return ff * T;
}

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, const std::shared_ptr<ExpressionVirtual> &fh) {
    std::list<std::shared_ptr<ExpressionVirtual>> ff = {fh};
    return ff * T;
}

template <int D> TestFunction<D> operator*(const VirtualParameter &cc, const TestFunction<D> &T) {
    TestFunction<D> multU;
    for (auto row_i : T.A) {
        for (auto &list : row_i) {
            for (auto &item : list.U) {
                item.addParameter(cc);
            }
        }
        multU.push(row_i);
    }
    return multU;
}

template <int D> TestFunction<D> operator*(const TestFunction<D> &T, const VirtualParameter &cc) { return cc * T; }

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
    auto [N, M] = T.size(); // N = number of components of T, M = number of test functions in each component

    TestFunction<D> gradU;

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);    // cannot take gradient of a matrix
        
        // Iterate over physical dimensions
        for (int d = 0; d < D; ++d) {

            // Get the ItemList corresponding to component i
            auto new_list = T.getList(i, 0);   // e.g. if T is a sum of two test functions, new_list will have two components

            // Iterate over the test functions in the item list
            for (auto &item : new_list.U) {
                item.du = nextDerivative(d, item.du);   // take first derivative if test function not already differentiated, otherwise take second derivative
            }
            int irow = T.isScalar() ? d : i;            // if T is a scalar, then grad(T) is a vector and the derivatives of T are placed in the rows, 
            int jrow = T.isScalar() ? 0 : d;            // otherwise grad(T) is a matrix and the components of T are placed in the rows, and the derivatives in the columns
            gradU.push({irow, jrow}, new_list);         // push the differentiated item list into the gradU matrix
        }
    }
    return gradU;
}


template <int D> TestFunction<D> div(const TestFunction<D> &T) {
    auto [N, M] = T.size();

    assert(M == 1);
    assert(N == D);
    TestFunction<D> divU;

    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);
        auto new_list = T.getList(i, 0);
        for (auto &item : new_list.U) {
            item.du = D1(i);
        }
        divU.push({0, 0}, new_list);
    }
    return divU;
}


template<int D> TestFunction<D> rot(const TestFunction<D> &T) {
    // T = [u0, u1, u2]

    auto [N, M] = T.size();
    assert(M == 1);             // we do not want to take rot of matrices
    
    // ONLY FOR 3D
    assert(N == 3);
    assert(D == 3);

    // Temporary test functions
    TestFunction<D> rotU1;
    TestFunction<D> rotU2;

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);

        // In component i, we want to return duj/dxk - duk/dxj 
        int j      = (i + 1) % 3;        // (if i = x then j = y, k = z, if i=y then j=z, k=x, if i=z then j=x, k=y)
        int k      = (i + 2) % 3;

        auto uj = T.getList(j, 0);       // component j of T
        auto uk = T.getList(k, 0);       // component k of T

        // if uj = uj1 + uj2 + ... + ukn, we need to differentiate each term separately
        assert(uj.U.size() == uk.U.size());
        for (int l = 0; l < uj.U.size(); ++l) {

            auto &itemj = uj.U[l];       // ujl
            auto &itemk = uk.U[l];       // ukl

            itemj.du = nextDerivative(k, itemj.du);     // take derivative dxk of ujl
            itemk.du = nextDerivative(j, itemk.du);     // take derivative dxj of ukl
        }

        // result should go into row i
        int irow = i;
        int jrow = 0;

        rotU1.push({irow, jrow}, uj);   // push duj/dxk 
        rotU2.push({irow, jrow}, uk);   // push duk/dxj 
    }
    
    return rotU1 - rotU2;       // return duj/dxk - duk/dxj
}

template <int D> TestFunction<D> rotgrad(const TestFunction<D> &T) {
    auto [N, M] = T.size();
    assert(M == 1);
    TestFunction<D> gradU;

    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);
        for (int d = 0; d < D; ++d) {
            auto new_list = T.getList(i, 0);
            for (auto &item : new_list.U) {
                item.du = rotgradD(d);
                if (d == 0)
                    item.c *= -1;
            }
            int irow = T.isScalar() ? d : i;
            int jrow = T.isScalar() ? 0 : d;
            gradU.push({irow, jrow}, new_list);
        }
    }
    return gradU;
}




template <int D> TestFunction<D> Eps(const TestFunction<D> &T) {
    auto [N, M] = T.size();
    assert(N == D && M == 1);
    TestFunction<D> gradU = grad(T);
    TestFunction<D> epsU  = 0.5 * gradU + 0.5 * transpose(gradU);
    return epsU;
}

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

template <int D> TestFunction<D> average(const TestFunction<D> &T, double v1 = 0.5, double v2 = 0.5) {
    double vv[2] = {v1, v2};
    int N        = T.nbRow();
    TestFunction<D> jumpU;
    for (int row = 0; row < N; ++row) {
        int M = T.nbCol();
        for (int col = 0; col < M; ++col) {
            for (int i = 0; i <= 1; ++i) {
                auto new_list = T.getList(row, col);
                for (auto &item : new_list.U) {
                    item.c *= vv[i];
                    item.face_side_ = i;
                }
                jumpU.push({row, col}, new_list);
            }
        }
    }

    return jumpU;
}

template <int D> TestFunction<D> jump(const TestFunction<D> &T) { return average(T, 1, -1); }

template <int d> TestFunction<d> operator*(const CutFEM_Rd<d> &cc, const TestFunction<d> &T) {
    assert(T.nbCol() == 1);
    int N = T.nbRow();

    bool scalar = (N == 1);
    int r       = (scalar) ? d : 1;
    TestFunction<d> resU;
    if (scalar) {
        for (int j = 0; j < d; ++j) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.addParameter(cc.get_parameter(j));
            }
            resU.push({j, 0}, new_list);
        }
    } else {
        assert(N == d); // column
        for (int j = 0; j < d; ++j) {
            auto new_list = T.getList(j, 0);
            for (auto &item : new_list.U) {
                item.addParameter(cc.get_parameter(j));
            }
            resU.push({0, 0}, new_list);
        }
    }
    return resU;
}

template <int d> TestFunction<d> operator*(const typename typeRd<d>::Rd &cc, const TestFunction<d> &T) {
    assert(T.nbCol() == 1);
    int N       = T.nbRow();
    bool scalar = (N == 1);
    int r       = (scalar) ? d : 1; // dim of the result
    TestFunction<d> resU;
    if (scalar) {
        for (int j = 0; j < d; ++j) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.c *= cc[j];
            }
            resU.push({j, 0}, new_list);
        }
    } else {
        assert(N == d); // column
        for (int j = 0; j < d; ++j) {
            auto new_list = T.getList(j, 0);
            for (auto &item : new_list.U) {
                item.c *= cc[j];
            }
            resU.push({0, 0}, new_list);
        }
    }
    return resU;
}

template <int D> TestFunction<D> gradS(const TestFunction<D> &T) {
    assert(T.nbCol() == 1);

    auto [N, M] = T.size();
    TestFunction<D> gradSU;
    TestFunction<D> gradU = grad(T);

    if (T.nbRow() == 1) {
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
                auto new_list = gradU.getList(j, 0);
                if (i == j) {
                    gradSU.push({i, 0}, new_list);
                }
                for (auto &item : new_list.U) {
                    item.c *= -1.;
                    item.addNormal(i);
                    item.addNormal(j);
                }
                gradSU.push({i, 0}, new_list);
            }
        }
    } else {
        for (int i_line = 0; i_line < N; ++i_line) {
            for (int i = 0; i < D; ++i) {
                for (int j = 0; j < D; ++j) {
                    auto new_list = gradU.getList(i_line, j);
                    if (i == j) {
                        gradSU.push({i_line, i}, new_list);
                    }
                    for (auto &item : new_list.U) {
                        item.c *= -1.;
                        item.addNormal(i);
                        item.addNormal(j);
                    }
                    gradSU.push({i_line, i}, new_list);
                }
            }
        }
    }

    return gradSU;
}

template <int D> TestFunction<D> dxS(const TestFunction<D> &T) {
    assert(T.nbCol() == 1 && T.nbRow() == 1);
    auto [N, M] = T.size();

    TestFunction<D> gradSU;
    TestFunction<D> gradU = grad(T);

    if (T.nbRow() == 1) {
        for (int j = 0; j < D; ++j) {
            auto new_list = gradU.getList(j, 0);
            if (0 == j) {
                gradSU.push({0, 0}, new_list);
            }
            for (auto &item : new_list.U) {
                item.c *= -1.;
                item.addNormal(0);
                item.addNormal(j);
            }
            gradSU.push({0, 0}, new_list);
        }

    } else {
        assert(0);
    }

    return gradSU;
}

template <int D> TestFunction<D> dyS(const TestFunction<D> &T) {
    assert(T.nbCol() == 1 && T.nbRow() == 1);
    auto [N, M] = T.size();

    TestFunction<D> gradSU;
    TestFunction<D> gradU = grad(T);

    if (T.nbRow() == 1) {
        for (int j = 0; j < D; ++j) {
            auto new_list = gradU.getList(j, 0);
            if (1 == j) {
                gradSU.push({0, 0}, new_list);
            }
            for (auto &item : new_list.U) {
                item.c *= -1.;
                item.addNormal(1);
                item.addNormal(j);
            }
            gradSU.push({0, 0}, new_list);
        }
    } else {
        assert(0);
    }

    return gradSU;
}



template <int D> TestFunction<D> average(const TestFunction<D> &U, const TestFunction<D> &V, int c1, int c2) {
    assert(U.nbCol() == V.nbCol());
    assert(U.nbRow() == V.nbRow());
    auto [N, M] = U.size();
    TestFunction<D> jumpU;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            auto new_list_U = U.getList(i, j);
            auto new_list_V = V.getList(i, j);
            for (auto &item : new_list_U.U) {
                item.c *= c1;
                item.face_side_ = 0;
            }
            for (auto &item : new_list_V.U) {
                item.c *= c2;
                item.face_side_ = 1;
            }
            jumpU.push({i, j}, new_list_U);
            jumpU.push({i, j}, new_list_V);
        }
    }
    return jumpU;
}

template <int D> TestFunction<D> jump(const TestFunction<D> &U, const TestFunction<D> &V) {
    return average(U, V, 1, -1);
}

template <int D>
TestFunction<D> average(const TestFunction<D> &T, const VirtualParameter &para1, const VirtualParameter &para2) {

    int N = T.nbRow();
    TestFunction<D> jumpU;
    for (int row = 0; row < N; ++row) {
        int M = T.nbCol();
        for (int col = 0; col < M; ++col) {
            for (int i = 0; i <= 1; ++i) {
                auto new_list = T.getList(row, col);
                const auto &p = (i == 0) ? para1 : para2;
                for (auto &item : new_list.U) {
                    item.face_side_ = i;
                    item.addParameter(p);
                }
                jumpU.push({row, col}, new_list);
            }
        }
    }

    return jumpU;
}

template <int D> TestFunction<D> average(const TestFunction<D> &T, const VirtualParameter &para1) {
    return average(T, para1, para1);
}

// template <int N> TestFunction<N> ln(const TestFunction<N> &F) {
//    return TestFunction<N>(F, f_ln);
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
//       u.du = D1(i);
//       k++;
//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(divU.A(0, 0)->getItem(k));
//          u = v;
//          u.c *= -1;
//          u.du = D1(j);
//          u.addNormal(i);
//          u.addNormal(j);
//          k++;
//       }
//    }
//    return divU;
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
//          u1.du = D1(2);
//       }
//       int kk = 1;
//       for (int j = 0; j < d; ++j) {
//          ItemTestFunction<d> &u(gradU.A(i, 0)->getItem(kk++));
//          u = v;
//          u.c *= -1;
//          u.du = D1(j);
//          u.addNormal(2);
//          u.addNormal(j);
//       }
//    }
//    return gradU;
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

typedef TestFunction<2> TestFunction2;
typedef TestFunction<3> TestFunction3;

#endif
