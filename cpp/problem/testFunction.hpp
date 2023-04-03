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
        return D1(c); // if test function is not differentiated, return first derivative in direction c
    else
        return D2(du - 1, c); // else, return second derivative in direction c
}

static int rotgradD(int i) {
    int op[2] = {op_dy, op_dx};
    return op[i];
}

// Given test function u with u=(ux, uy, uz) where c denotes the component of u, and

void f_id(RNMK_ &x, int cu, int du);
void f_ln(RNMK_ &x, int cu, int du);

template <typeMesh M> struct TestFunction;

/**
 * @brief ItemTestFunction is the core class representing a test function
 * @tparam N Physical dimension
 * @note This class holds all necessary information about the test function, 
 * such as the coefficient, the component of the vector field, the derivative index,
 * the time derivative index, the normal and conormal components, the domain id
 * and the face side index. It also holds the pointer to the root TestFunction object
 * and the pointer to the function that evaluates the test function.
 */
template <typeMesh M> struct ItemTestFunction {
    using mesh_t         = M;
    using fespace_t      = GFESpace<mesh_t>;
    using testfunction_t = TestFunction<mesh_t>;

    static const int D = mesh_t::D;

    double c;       ///< c – scalar multiple of the test function
    int cu;         ///< cu – component of root test function root_fun_p //?
    int du;         ///< du – derivative index of root test function, e.g. dx, dy or dz
    int dtu;        ///< dtu – time derivative index

    std::vector<int> ar_nu, conormal;                       ///< arrays of normal and conormal components of test function //? is it which normal component this ItemTestFunction is from root test function? e.g. if ar_nu = {1}, then this ItemTestFunction is the normal component of root test function in direction y?
    std::vector<const VirtualParameter *> coefu;        ///< CutFEM parameter multiple of the test function
    int domain_id_,
        face_side_;                             ///< e.g. average(u) saves the same test function but with different domain_id/face_side_
    std::shared_ptr<ExpressionVirtual> expru = nullptr; ///< if multiplied by an expression function
    const fespace_t *fespace                 = nullptr; ///< space of integration
    const testfunction_t *root_fun_p         = nullptr; ///< the root test function object     ///< pointer to the root test function

    void (*pfun)(RNMK_ &, int, int) = f_id;                 ///< //? what does this do?

    /**
     * @brief Default constructor
     * @note Creates a test function with coefficient 0, no derivative, no time derivative,
     * no normal or conormal components, default domain id and default face side index (-1).
     */
    ItemTestFunction() : c(0.), cu(-1), du(-1), dtu(-1), domain_id_(-1), face_side_(-1) {}

    /**
     * @brief Construct an ItemTestFunction object based on data
     * 
     * @param fes Finite element space GFESpace<Mesh>
     * @param cc Scalar multiple
     * @param i Component of root test function 
     * @param j Derivative index of root test function
     * @param tu Time derivative index of root test function
     * @param dd Domain id
     * @param ff Root test function
     */
    ItemTestFunction(const fespace_t *fes, double cc, int i, int j, int tu, int dd, testfunction_t *ff)
        : c(cc), cu(i), du(j), dtu(tu), coefu(0), domain_id_(dd), face_side_(-1), fespace(fes), root_fun_p(ff){};
    
    /**
     * @brief Copy constructor
     * 
     * @param F ItemTestFunction object to copy
     */
    ItemTestFunction(const ItemTestFunction &F)
        : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu), conormal(F.conormal), domain_id_(F.domain_id_),
          face_side_(F.face_side_), expru(F.expru), fespace(F.fespace), root_fun_p(F.root_fun_p) {
        for (int i = 0; i < F.coefu.size(); ++i)
            coefu.push_back(F.coefu[i]);
    }
    /**
     * @brief //? What does this do?
     * 
     * @param F 
     * @param f 
     */
    ItemTestFunction(const ItemTestFunction &F, void (*f)(RNMK_ &, int, int))
        : c(F.c), cu(F.cu), du(F.du), dtu(F.dtu), ar_nu(F.ar_nu), conormal(F.conormal), domain_id_(F.domain_id_),
          face_side_(F.face_side_), expru(F.expru), fespace(F.fespace), root_fun_p(F.root_fun_p), pfun(f) {
        for (int i = 0; i < F.coefu.size(); ++i)
            coefu.push_back(F.coefu[i]);
    }

    /**
     * @brief Assignment operator
     * 
     * @param L ItemTestFunction object to copy
     * @return ItemTestFunction& 
     */
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

    void addNormal(int i) { ar_nu.push_back(i); }                           ///< Adds which normal component this ItemTestFunction corresponds to
    void addParameter(const VirtualParameter &x) { coefu.push_back(&x); }   ///< Adds a parameter multiple to the test function
    void addParameter(const VirtualParameter *x) { coefu.push_back(x); }    ///< Adds a parameter multiple to the test function

    /**
     * @brief Checks if two ItemTestFunction objects are equal
     * 
     * @param F 
     * @return true 
     * @return false 
     */
    bool operator==(const ItemTestFunction &F) const {
        if (root_fun_p == F.root_fun_p && coefu.size() == 0 && F.coefu.size() == 0 && cu == F.cu && du == F.du &&
            dtu == F.dtu && face_side_ == F.face_side_ && domain_id_ == F.domain_id_ && ar_nu == F.ar_nu &&
            conormal == F.conormal && expru.get() == F.expru.get() && fespace == F.fespace)
            return true;
        else
            return false;
    }

    /**
     * @brief Write ItemTestFunction to output stream
     * 
     * @param f ostream object
     * @param u ItemTestFunction object
     */
    friend std::ostream &operator<<(std::ostream &f, const ItemTestFunction &u) {
        std::string n[3] = {"nx", "ny", "nz"};
        f << " fespace_t => " << u.fespace << "\t";
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
template <typeMesh M> struct ItemList {

    using mesh_t         = M;
    using fespace_t      = GFESpace<mesh_t>;
    using item_t         = ItemTestFunction<mesh_t>;
    using testfunction_t = TestFunction<mesh_t>;
    static const int D   = mesh_t::D;

    std::vector<item_t> U;      ///< vector of Items

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

    /**
     * @brief Construct a new ItemList object
     *
     * @param Vh Finite Element Space
     * @param cc Scalar multiple of the test function
     * @param i  Component index
     * @param j  Derivative index (op_id = 0, op_dx = 1, op_dy = 2, op_dz = 3 etc)
     * @param dd Domain index: -1 equals main domain if only one is defined
     * @param ff Pointer to the test function
     */
    ItemList(const fespace_t &Vh, double cc, int i, int j, int dd, testfunction_t *ff) {
        U.push_back(item_t(&Vh, cc, i, j, 0, dd, ff));
    }

    int size() const { return U.size(); }
    const item_t &operator()(int i) const { return U[i]; }
    const item_t &getItem(int i) const { return U[i]; }
    item_t &getItem(int i) { return U[i]; }
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
template <typeMesh M> struct TestFunction {

    using mesh_t     = M;
    using fespace_t  = GFESpace<mesh_t>;
    using Rd         = typename mesh_t::Rd;
    using item_t     = ItemTestFunction<mesh_t>;
    using itemList_t = ItemList<mesh_t>;

    static const int D = mesh_t::D;

    std::vector<std::vector<itemList_t>> A; ///< Matrix of ItemLists

    TestFunction() {} ///< Default constructor

    /**
     * @brief Construct a new Test Function object
     *
     * @param Vh Finite Element space
     * @param d Number of components of test function
     */
    TestFunction(const fespace_t &Vh, int d) {
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
    TestFunction(const fespace_t &Vh, int d, int comp0) {
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
    TestFunction(const fespace_t &Vh, int d, int comp0, int domm) {
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

template <typeMesh mesh_t> void simplify(TestFunction<mesh_t> &T) {

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

    auto is_nul = [](typename TestFunction<mesh_t>::item_t item) -> bool {
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

template <typeMesh mesh_t> TestFunction<mesh_t> transpose(TestFunction<mesh_t> &T) {
    TestFunction<mesh_t> Ut;
    for (int i = 0; i < T.nbRow(); ++i) {
        for (int j = 0; j < T.nbCol(); ++j) {
            auto list = T.getList(i, j);
            Ut.push({j, i}, list);
        }
    }
    return Ut;
}

template <typeMesh mesh_t>
TestFunction<mesh_t> operator+(const TestFunction<mesh_t> &f1, const TestFunction<mesh_t> &f2) {
    assert(f1.size() == f2.size());
    TestFunction<mesh_t> sum_fi;
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

template <typeMesh mesh_t>
TestFunction<mesh_t> operator-(const TestFunction<mesh_t> &f1, const TestFunction<mesh_t> &f2) {
    assert(f1.size() == f2.size());
    TestFunction<mesh_t> sub_fi;
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

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, double cc) {
    TestFunction<mesh_t> mulT;
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

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(double cc, const TestFunction<mesh_t> &T) { return T * cc; }

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const BaseVector &vec) {
    auto [N, M] = T.size();
    int D       = mesh_t::D;
    TestFunction<mesh_t> Un;

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

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const Conormal &vec) {
    auto [N, M] = T.size();

    TestFunction<mesh_t> Un;
    int D = mesh_t::D;
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

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const Normal_Component &cc) {
    TestFunction<mesh_t> multU;
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

template <typeMesh mesh_t, typename Expr>
TestFunction<mesh_t> operator*(const std::list<std::shared_ptr<Expr>> &fh, const TestFunction<mesh_t> &T) {
    auto [N, M] = T.size();
    assert(M == 1);
    int D = mesh_t::D;
    TestFunction<mesh_t> Un;
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

template <typeMesh mesh_t>
TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const std::list<std::shared_ptr<ExpressionVirtual>> &fh) {
    return fh * T;
}

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(ExpressionFunFEM<mesh_t> &fh, const TestFunction<mesh_t> &T) {
    auto fh_p                                               = std::make_shared<ExpressionFunFEM<mesh_t>>(fh);
    std::list<std::shared_ptr<ExpressionFunFEM<mesh_t>>> ff = {fh_p};
    return ff * T;
}

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, ExpressionFunFEM<mesh_t> &fh) {
    return fh * T;
}

template <typeMesh mesh_t>
TestFunction<mesh_t> operator*(const std::shared_ptr<ExpressionVirtual> &fh, const TestFunction<mesh_t> &T) {
    std::list<std::shared_ptr<ExpressionVirtual>> ff = {fh};
    return ff * T;
}

template <typeMesh mesh_t>
TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const std::shared_ptr<ExpressionVirtual> &fh) {
    std::list<std::shared_ptr<ExpressionVirtual>> ff = {fh};
    return ff * T;
}

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const VirtualParameter &cc, const TestFunction<mesh_t> &T) {
    TestFunction<mesh_t> multU;
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

template <typeMesh mesh_t> TestFunction<mesh_t> operator*(const TestFunction<mesh_t> &T, const VirtualParameter &cc) {
    return cc * T;
}

template <typeMesh mesh_t> TestFunction<mesh_t> dt(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);

    TestFunction<mesh_t> dt_U;
    for (auto comp_i : T.A) {
        for (auto &item_list : comp_i) {
            for (auto &item : item_list.U)
                item.dtu = op_dx;
        }
        dt_U.push(comp_i);
    }
    return dt_U;
}

template <typeMesh mesh_t> TestFunction<mesh_t> dx(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);

    TestFunction<mesh_t> dx_U;
    for (auto comp_i : T.A) {
        for (auto &item_list : comp_i) {
            for (auto &item : item_list.U)
                item.du = op_dx;
        }
        dx_U.push(comp_i);
    }
    return dx_U;
}

template <typeMesh mesh_t> TestFunction<mesh_t> dy(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);

    TestFunction<mesh_t> dx_U;
    for (auto comp_i : T.A) {
        for (auto &item_list : comp_i) {
            for (auto &item : item_list.U)
                item.du = op_dy;
        }
        dx_U.push(comp_i);
    }
    return dx_U;
}

template <typeMesh mesh_t> TestFunction<mesh_t> dz(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);

    TestFunction<mesh_t> dx_U;
    for (auto comp_i : T.A) {
        for (auto &item_list : comp_i) {
            for (auto &item : item_list.U)
                item.du = op_dz;
        }
        dx_U.push(comp_i);
    }
    return dx_U;
}

/**
 * @brief Compute the gradient of a TestFunction
 * 
 * @tparam D Physical dimension
 * @param T TestFunction
 * @return TestFunction<D> 
 */
template <typeMesh mesh_t> TestFunction<mesh_t> grad(const TestFunction<mesh_t> &T) {
    auto [N, M] = T.size(); // N = number of components of T, M = number of columns of T
    int D       = mesh_t::D;
    TestFunction<mesh_t> gradU;

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1); // cannot take gradient of a matrix

        // Iterate over physical dimensions
        for (int d = 0; d < D; ++d) {

            // Get the ItemList corresponding to component i
            auto new_list = T.getList(i, 0);   // e.g. if T is a sum of two test functions, T = T1 + T2, new_list will have two components {T1, T2}

            // Iterate over the test functions in the item list
            for (auto &item : new_list.U) {
                item.du = nextDerivative(d, item.du);   // take first derivative if test function not already differentiated, otherwise take second derivative (or mixed derivative etc.)
            }

            int irow = T.isScalar() ? d : i;            // if T is a scalar, then grad(T) is a vector and the derivatives of T are placed in the rows, 
            int jrow = T.isScalar() ? 0 : d;            // otherwise grad(T) is a matrix and the components of T are placed in the rows, and the derivatives in the columns
            gradU.push({irow, jrow}, new_list);         // push the differentiated item list into the gradU matrix
        }
    }
    return gradU;
}

/**
 * @brief Compute the divergence of a TestFunction
 * 
 * @tparam D Physical dimension
 * @param T TestFunction
 * @return TestFunction<D> with 1 component
 */template <typeMesh mesh_t> TestFunction<mesh_t> div(const TestFunction<mesh_t> &T) {
    auto [N, M] = T.size(); // N = number of components of T
    int D       = mesh_t::D;
    assert(M == 1);         // We do not want to take div of matrices
    assert(N == D);         // The number of components of T must be equal to the physical dimension
    TestFunction<mesh_t> divU;   // Temporary variable to store the divergence of T

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);    // cannot take divergence of a matrix

        // Get the ItemList corresponding to component i
        auto new_list = T.getList(i, 0);

        // Iterate over the test functions in the item list
        for (auto &item : new_list.U) {

            // Take the derivative with the same index as the component of the test function
            item.du = D1(i);
        }

        divU.push({0, 0}, new_list);    // Result is a scalar function, so we push the differentiated item list into the same spot in the matrix
    }

    return divU;
}

/**
 * @brief Compute the three-dimensional curl of a TestFunction
 * 
 * @tparam D Physical dimension (=3)
 * @param T TestFunction
 * @return TestFunction<D> with 3 components
 */
template <typeMesh mesh_t> TestFunction<mesh_t> curl(const TestFunction<mesh_t> &T) {
    // T = [u0, u1, u2]
    int D       = mesh_t::D;
    auto [N, M] = T.size();
    assert(M == 1);             // we do not want to take rot of matrices
    
    // Particularly for 3D
    assert(N == 3);
    assert(D == 3);

    // Temporary test functions
    TestFunction<mesh_t> rotU1;
    TestFunction<mesh_t> rotU2;

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);

        // In component i, we want to return duj/dxk - duk/dxj
        int j = (i + 1) % 3; // (if i = x then j = y, k = z, if i=y then j=z, k=x, if i=z then j=x, k=y)
        int k = (i + 2) % 3;

        auto uj = T.getList(j, 0); // component j of T
        auto uk = T.getList(k, 0); // component k of T

        // if uj = uj1 + uj2 + ... + ukn, we need to differentiate each term separately
        assert(uj.U.size() == uk.U.size());
        for (int l = 0; l < uj.U.size(); ++l) {

            auto &itemj = uj.U[l]; // ujl
            auto &itemk = uk.U[l]; // ukl

            itemj.du = nextDerivative(k, itemj.du); // take derivative dxk of ujl
            itemk.du = nextDerivative(j, itemk.du); // take derivative dxj of ukl
        }

        // result should go into row i
        int irow = i;
        int jrow = 0;

        rotU1.push({irow, jrow}, uj); // push duj/dxk
        rotU2.push({irow, jrow}, uk); // push duk/dxj
    }

    return rotU1 - rotU2; // return duj/dxk - duk/dxj
}

/**
 * @brief Compute two-dimensional curl of a TestFunction
 * @note This operator is different if TestFunction is a vector or a scalar.
 * If u vector: rotgrad(u) = -dy(ux) + dx(uy)
 * If u scalar: rotgrad(u) = [-dy(u), dx(u)]
 * @tparam mesh_t Mesh
 * @param T TestFunction
 * @return TestFunction<D> with 2 components if T is a scalar, and 1 component if T is a vector
 */
template <typeMesh mesh_t> TestFunction<mesh_t> rotgrad(const TestFunction<mesh_t> &T) {
    auto [N, M] = T.size();     // N = number of components of T
    int D       = mesh_t::D;
    assert(M == 1);             // We do not want to take rotgrad of matrices
    TestFunction<mesh_t> gradU;      // Temporary variable to store rotgrad of T

    assert(T.isScalar());       // only implemented when T is scalar (for now)

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);

        // Iterate over physical dimensions
        for (int d = 0; d < D; ++d) {

            // Get the ItemList corresponding to component i
            auto new_list = T.getList(i, 0);

            // Iterate over the test functions in the item list
            for (auto &item : new_list.U) {

                // Take dy if d=0 and dx if d=1
                item.du = rotgradD(d);

                // If d=0, multiply by -1
                if (d == 0)
                    item.c *= -1;
            }

            // If T is a scalar, we want to return a vector, so we push the differentiated item list into the row corresponding to the physical dimension
            // If T is a vector, we want to return a scalar, so we push the differentiated item list into the same spot in the matrix
            //! Erik, is the vector case correct here? Shouldn't it be
            /*
            int irow = T.isScalar() ? d : 0;    
            int jrow = T.isScalar() ? 0 : 0;    
            */
            int irow = T.isScalar() ? d : 0;    
            int jrow = T.isScalar() ? 0 : 0;    
            gradU.push({irow, jrow}, new_list);
        }
    }
    return gradU;
}

/**
 * @brief Compute the symmetric part of the strain tensor of a TestFunction //?
 * 
 * @tparam mesh_t Mesh
 * @param T Test function
 * @return TestFunction<mesh_t> 
 */
template <typeMesh mesh_t> TestFunction<mesh_t> Eps(const TestFunction<mesh_t> &T) {
    auto [N, M] = T.size();
    int D       = mesh_t::D;
    assert(N == D && M == 1);
    TestFunction<mesh_t> gradU = grad(T);
    TestFunction<mesh_t> epsU  = 0.5 * gradU + 0.5 * transpose(gradU);
    return epsU;
}

/**
 * @brief Compute average of a TestFunction
 * 
 * @tparam mesh_t Mesh
 * @param T TestFunction
 * @param v1 Weight of function on first side of the face
 * @param v2 Weight of function on second side of the face
 * @return TestFunction<mesh_t> with the same number of components as T
 */
template <typeMesh mesh_t>
TestFunction<mesh_t> average(const TestFunction<mesh_t> &T, double v1 = 0.5, double v2 = 0.5) {
    double vv[2] = {v1, v2};        // Array of weights
    int N        = T.nbRow();       // Number of components of T
    TestFunction<mesh_t> jumpU;     // Temporary variable to store average of T

    // Iterate over number of components of T
    for (int row = 0; row < N; ++row) {
        int M = T.nbCol();                // Number of columns of T

        // Iterate over columns of T
        for (int col = 0; col < M; ++col) {

            // Iterate over the two sides of the face
            for (int i = 0; i <= 1; ++i) {

                // Get the ItemList corresponding to component (row, col) in the matrix of T
                auto new_list = T.getList(row, col);

                // Iterate over the test functions in the item list
                for (auto &item : new_list.U) {

                    // Multiply by the weight
                    item.c *= vv[i];        

                    // Set the face_side_ to i
                    item.face_side_ = i;    //? Are we not able to take average across interface instead of face?
                }

                jumpU.push({row, col}, new_list);
            }
        }
    }

    return jumpU;
}

/**
 * @brief Compute jump of a TestFunction
 * @note This is the same as average, but with weights 1 and -1
 * @tparam D Physical dimension
 * @param T Test function
 * @return TestFunction<D> with the same number of components as T
 */
template <typeMesh mesh_t> TestFunction<mesh_t> jump(const TestFunction<mesh_t> &T) { return average(T, 1, -1); }

template <typeMesh mesh_t>
TestFunction<mesh_t> operator*(const CutFEM_Rd<mesh_t::D> &cc, const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);
    int N       = T.nbRow();
    int D       = mesh_t::D;
    bool scalar = (N == 1);
    int r       = (scalar) ? D : 1;
    TestFunction<mesh_t> resU;
    if (scalar) {
        for (int j = 0; j < D; ++j) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.addParameter(cc.get_parameter(j));
            }
            resU.push({j, 0}, new_list);
        }
    } else {
        assert(N == D); // column
        for (int j = 0; j < D; ++j) {
            auto new_list = T.getList(j, 0);
            for (auto &item : new_list.U) {
                item.addParameter(cc.get_parameter(j));
            }
            resU.push({0, 0}, new_list);
        }
    }
    return resU;
}

template <typeMesh mesh_t>
TestFunction<mesh_t> operator*(const typename mesh_t::Rd &cc, const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);
    int N       = T.nbRow();
    int D       = mesh_t::D;
    bool scalar = (N == 1);
    int r       = (scalar) ? D : 1; // dim of the result
    TestFunction<mesh_t> resU;
    if (scalar) {
        for (int j = 0; j < D; ++j) {
            auto new_list = T.getList(0, 0);
            for (auto &item : new_list.U) {
                item.c *= cc[j];
            }
            resU.push({j, 0}, new_list);
        }
    } else {
        assert(N == D); // column
        for (int j = 0; j < D; ++j) {
            auto new_list = T.getList(j, 0);
            for (auto &item : new_list.U) {
                item.c *= cc[j];
            }
            resU.push({0, 0}, new_list);
        }
    }
    return resU;
}

/**
 * @brief Compute the trangential gradient of a TestFunction
 * @note gradS = (I - nn^T) grad
 * @tparam mesh_t Mesh
 * @param T TestFunction
 * @return TestFunction<D> with D components
 */
template <typeMesh mesh_t> TestFunction<mesh_t> gradS(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1);             // T must be a vector or scalar
    int D       = mesh_t::D;
    auto [N, M] = T.size();             // Get number of components of T
    TestFunction<mesh_t> gradSU;             // Temporary variable to store gradS of T
    TestFunction<mesh_t> gradU = grad(T);    // Compute grad of T

    // If T is a scalar, gradS is a vector
    if (T.nbRow() == 1) {

        // Iterate over the physical dimensions
        for (int i = 0; i < D; ++i) {

            // Iterate over the physical dimensions
            for (int j = 0; j < D; ++j) {

                // Get the ItemList corresponding to component j in the gradU vector
                auto new_list = gradU.getList(j, 0);

                // Compute the I * grad part
                if (i == j) {
                    gradSU.push({i, 0}, new_list);
                }

                // Compute the (- nn^T) * grad part
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

template <typeMesh mesh_t> TestFunction<mesh_t> dxS(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1 && T.nbRow() == 1);
    auto [N, M] = T.size();
    int D       = mesh_t::D;
    TestFunction<mesh_t> gradSU;
    TestFunction<mesh_t> gradU = grad(T);

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

template <typeMesh mesh_t> TestFunction<mesh_t> dyS(const TestFunction<mesh_t> &T) {
    assert(T.nbCol() == 1 && T.nbRow() == 1);
    auto [N, M] = T.size();
    int D       = mesh_t::D;

    TestFunction<mesh_t> gradSU;
    TestFunction<mesh_t> gradU = grad(T);

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

/**
 * @brief Compute the cross product of a normal and a TestFunction
 * 
 * @tparam mesh_t Mesh
 * @param n Normal
 * @param T Test function
 * @return TestFunction<mesh_t> of dimension 3 
 */
template <typeMesh mesh_t> TestFunction<mesh_t> cross(const Normal& n, const TestFunction<mesh_t> &T) {
    // T = [u0, u1, u2]
    int D       = mesh_t::D;
    auto [N, M] = T.size();
    assert(M == 1);             // we do not want to take rot of matrices
    
    // Particularly for 3D
    assert(N == 3);
    assert(D == 3);

    // Temporary test functions
    TestFunction<mesh_t> rotU1;
    TestFunction<mesh_t> rotU2;

    // Iterate over number of components of T
    for (int i = 0; i < N; ++i) {
        assert(T.nbCol(i) == 1);

        // In component i, we want to return uj*nk - uk*nj
        int j = (i + 1) % 3; // (if i = x then j = y, k = z, if i=y then j=z, k=x, if i=z then j=x, k=y)
        int k = (i + 2) % 3;

        auto uj = T.getList(j, 0); // component j of T
        auto uk = T.getList(k, 0); // component k of T

        // if uj = uj1 + uj2 + ... + ukn, we need to differentiate each term separately
        assert(uj.U.size() == uk.U.size());
        for (int l = 0; l < uj.U.size(); ++l) {

            auto &itemj = uj.U[l]; // ujl
            auto &itemk = uk.U[l]; // ukl

            itemj.addNormal(k); // add nk to ujl
            itemk.addNormal(j); // add nj to ukl
        }

        // result should go into row i
        int irow = i;
        int jrow = 0;

        rotU1.push({irow, jrow}, uj); // push duj/dxk
        rotU2.push({irow, jrow}, uk); // push duk/dxj
    }

    return rotU1 - rotU2; // return duj/dxk - duk/dxj

}

template <typeMesh mesh_t> TestFunction<mesh_t> cross(const TestFunction<mesh_t> &T, const Normal& n) {return -1*cross(n, T);}

template <typeMesh mesh_t>
TestFunction<mesh_t> average(const TestFunction<mesh_t> &U, const TestFunction<mesh_t> &V, int c1, int c2) {
    assert(U.nbCol() == V.nbCol());
    assert(U.nbRow() == V.nbRow());
    auto [N, M] = U.size();
    TestFunction<mesh_t> jumpU;

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

template <typeMesh mesh_t> TestFunction<mesh_t> jump(const TestFunction<mesh_t> &U, const TestFunction<mesh_t> &V) {
    return average(U, V, 1, -1);
}

template <typeMesh mesh_t>
TestFunction<mesh_t> average(const TestFunction<mesh_t> &T, const VirtualParameter &para1,
                             const VirtualParameter &para2) {

    int N = T.nbRow();
    TestFunction<mesh_t> jumpU;
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

template <typeMesh mesh_t> TestFunction<mesh_t> average(const TestFunction<mesh_t> &T, const VirtualParameter &para1) {
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

// typedef TestFunction<2> TestFunction2;
// typedef TestFunction<3> TestFunction3;

#endif
