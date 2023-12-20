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
#ifndef _ITEMVF_HPP
#define _ITEMVF_HPP
#include "testFunction.hpp"
#include <list>

template <typeMesh M> struct ItemVF {

    using mesh_t              = M;
    using fespace_t           = GFESpace<mesh_t>;
    using item_testfunction_t = ItemTestFunction<mesh_t>;
    using Rd                  = typename mesh_t::Rd;
    static const int D        = mesh_t::D;

    double c;
    int cu, du, cv, dv;
    std::vector<int> ar_nu, ar_nv;
    std::vector<int> conormalU_, conormalV_;
    int face_sideU_, face_sideV_;
    int domainU_id_, domainV_id_;

    std::vector<const VirtualParameter *> coefu, coefv;
    int dtu, dtv;
    std::shared_ptr<ExpressionVirtual> expru = nullptr;
    std::shared_ptr<ExpressionVirtual> exprv = nullptr;

    fespace_t const *fespaceU = nullptr;
    fespace_t const *fespaceV = nullptr;

    void (*pfunU)(RNMK_ &, int, int) = f_id;
    void (*pfunV)(RNMK_ &, int, int) = f_id;

    ItemVF()
        : c(0.), cu(-1), du(-1), cv(-1), dv(-1), face_sideU_(-1), face_sideV_(-1), domainU_id_(-1), domainV_id_(-1),
          dtu(-1), dtv(-1) {}
    ItemVF(double cc, int i, int j, int k, int l)
        : c(cc), cu(i), du(j), cv(k), dv(l), face_sideU_(-1), face_sideV_(-1), domainU_id_(-1), domainV_id_(-1),
          dtu(-1), dtv(-1) {}
    ItemVF(double cc, int i, int j, int k, int l, const std::vector<int> &nn, const std::vector<int> &mm)
        : c(cc), cu(i), du(j), cv(k), dv(l), ar_nu(nn), ar_nv(mm), face_sideU_(-1), face_sideV_(-1), domainU_id_(-1),
          domainV_id_(-1), dtu(-1), dtv(-1) {}
    ItemVF(const ItemVF &U) : ItemVF(U.c, U.cu, U.du, U.cv, U.dv, U.ar_nu, U.ar_nv) {
        conormalU_ = U.conormalU_;
        conormalV_ = U.conormalV_;

        face_sideU_ = U.face_sideU_;
        face_sideV_ = U.face_sideV_;

        domainU_id_ = U.domainU_id_;
        domainV_id_ = U.domainV_id_;

        for (int i = 0; i < U.coefu.size(); ++i)
            coefu.push_back(U.coefu[i]);
        for (int i = 0; i < U.coefv.size(); ++i)
            coefv.push_back(U.coefv[i]);

        dtu      = U.dtu;
        dtv      = U.dtv;
        expru    = U.expru;
        exprv    = U.exprv;
        fespaceU = U.fespaceU;
        fespaceV = U.fespaceV;
        pfunU    = U.pfunU;
        pfunV    = U.pfunV;
    }

    ItemVF(const item_testfunction_t &U, const item_testfunction_t &V)
        : ItemVF(U.c * V.c, U.cu, U.du, V.cu, V.du, U.ar_nu, V.ar_nu) {
        conormalU_ = U.conormal;
        conormalV_ = V.conormal;

        face_sideU_ = U.face_side_;
        face_sideV_ = V.face_side_;

        domainU_id_ = U.domain_id_;
        domainV_id_ = V.domain_id_,

        coefu = U.coefu;
        coefv = V.coefu;
        dtu   = U.dtu;
        dtv   = V.dtu;

        expru = U.expru;
        exprv = V.expru;

        fespaceU = U.fespace;
        fespaceV = V.fespace;

        pfunU = U.pfun;
        pfunV = V.pfun;
    }

    bool operator==(const ItemVF &F) {

        if (cu == F.cu && cv == F.cv && du == F.du && dv == F.dv && F.face_sideU_ == face_sideU_ &&
            face_sideV_ == F.face_sideV_ && dtu == F.dtu && dtv == F.dtv && domainU_id_ == F.domainU_id_ &&
            domainV_id_ == F.domainV_id_ && fespaceU == F.fespaceU && fespaceV == F.fespaceV &&
            expru.get() == F.expru.get() && exprv.get() == F.exprv.get() && ar_nu == F.ar_nu && ar_nv == F.ar_nv &&
            conormalU_ == F.conormalU_ && conormalV_ == F.conormalV_) {
        } else
            return false;
        return true;
    }

  public:
    void applyFunNL(RNMK_ &bfu, RNMK_ &bfv) const {
        pfunU(bfu, cu, du);
        pfunV(bfv, cv, dv);
    }

    // FOR NEW VERSION
    double computeCoefElement(double h, double meas, double measK, double measCut, int domain) const {
        R val = 1;
        for (int l = 0; l < 2; ++l) {
            const std::vector<const VirtualParameter *> &listCoef = (l == 0) ? coefu : coefv;
            for (int i = 0; i < listCoef.size(); ++i) {
                val *= listCoef[i]->evaluate(domain, h, meas, measK, measCut);
            }
        }
        return val;
    }
    double computeCoefInterface(double h, double meas, double measK, double measCut, int domi, int domj) const {
        R val = 1;
        for (int l = 0; l < 2; ++l) {
            const std::vector<const VirtualParameter *> &listCoef = (l == 0) ? coefu : coefv;
            int dom                                               = (l == 0) ? domi : domj;
            for (int i = 0; i < listCoef.size(); ++i) {
                val *= listCoef[i]->evaluate(dom, h, meas, measK, measCut);
            }
        }
        return val;
    }
    double computeCoefInterface(double h, double meas, double measK, std::array<double, 2> measCut,
                                std::pair<int, int> domi) const {
        R val = 1;
        for (int l = 0; l < 2; ++l) {
            const std::vector<const VirtualParameter *> &listCoef = (l == 0) ? coefu : coefv;
            int dom                                               = (l == 0) ? domi.first : domi.second;
            for (int i = 0; i < listCoef.size(); ++i) {
                val *= listCoef[i]->evaluate(dom, h, meas, measK, measCut[i]);
            }
        }
        return val;
    }
    double evaluateFunctionOnBackgroundMesh(int k, int dom, Rd mip, double tid, const R *normal = nullptr) const {
        return ((expru) ? expru->GevalOnBackMesh(k, dom, mip, tid, normal) : 1) *
               ((exprv) ? exprv->GevalOnBackMesh(k, dom, mip, tid, normal) : 1);
    }
    double evaluateFunctionOnBackgroundMesh(const std::pair<int, int> &k, const std::pair<int, int> &dom, Rd mip,
                                            double tid, const R *normal = nullptr) const {
        return ((expru) ? expru->GevalOnBackMesh(k.first, dom.first, mip, tid, normal) : 1) *
               ((exprv) ? exprv->GevalOnBackMesh(k.second, dom.second, mip, tid, normal) : 1);
    }
    double evaluateFunctionOnBackgroundMesh(const std::pair<int, int> &k, const std::pair<int, int> &dom, Rd mip,
                                            double tid, const std::pair<Rd, Rd> normal) const {
        return ((expru) ? expru->GevalOnBackMesh(k.first, dom.first, mip, tid, normal.first) : 1) *
               ((exprv) ? exprv->GevalOnBackMesh(k.second, dom.second, mip, tid, normal.second) : 1);
    }

    double computeCoefFromNormal(const Rd normal) const {
        return computeCoeffFromArray(ar_nu, normal) * computeCoeffFromArray(ar_nv, normal);
        ;
    }
    double computeCoefFromConormal(const Rd conormal) const {
        return computeCoeffFromArray(conormalU_, conormal) * computeCoeffFromArray(conormalV_, conormal);
    }
    double computeCoefFromNormal(const Rd normal0, const Rd normal1) const {
        return computeCoeffFromArray(ar_nu, normal0) * computeCoeffFromArray(ar_nv, normal1);
        ;
    }
    double computeCoefFromConormal(const Rd conormalu, const Rd conormalv) const {
        return computeCoeffFromArray(conormalU_, conormalu) * computeCoeffFromArray(conormalV_, conormalv);
    }
    int onWhatElementIsTrialFunction(std::vector<int> k) const {
        if (k.size() == 1) {
            return k[0];
        } else {
            assert(face_sideU_ != -1);
            return (face_sideU_ == 0) ? k[0] : k[1];
        }
    }
    int onWhatElementIsTestFunction(std::vector<int> k) const {
        if (k.size() == 1) {
            return k[0];
        } else {
            assert(face_sideV_ != -1);
            return (face_sideV_ == 0) ? k[0] : k[1];
        }
    }
    template <typename T> T onWhatElementIsTrialFunction(T ki, T kj) const {
        assert(face_sideU_ != -1);
        return (face_sideU_ == 0) ? ki : kj;
    }
    template <typename T> T onWhatElementIsTestFunction(T ki, T kj) const {
        assert(face_sideV_ != -1);
        return (face_sideV_ == 0) ? ki : kj;
    }
    int get_domain_test_function() const { return domainV_id_; }
    int get_domain_trial_function() const { return domainU_id_; }

    bool on(int d) const { return ((domainU_id_ == domainV_id_) && (domainU_id_ == -1 || domainU_id_ == d)); }

  private:
    double computeCoeffFromArray(const std::vector<int> &array_idx, const double *v) const {
        double val = 1;
        for (const auto &i : array_idx)
            val *= v[i];
        return val;
    }

  public:
    // ------------------------------------------

    friend std::ostream &operator<<(std::ostream &f, const ItemVF &u) {
        std::string n[3] = {"nx", "ny", "nz"};
        // f << " FESpaces => " << u.fespaceU << " and " << u.fespaceV << "\t";
        f << u.c << "\t" << whichOperator(u.dtu) << whichOperator(u.du, u.cu);
        for (int i = 0; i < u.ar_nu.size(); ++i)
            f << " * " << n[u.ar_nu[i]];
        // for(int i=0;i<u.coefu.size();++i) f << " * " << u.coefu[i];
        f << " * " << whichOperator(u.dtv) << whichOperatorV(u.dv, u.cv);
        for (int i = 0; i < u.ar_nv.size(); ++i)
            f << " * " << n[u.ar_nv[i]];
        // for(int i=0;i<u.coefv.size();++i) f << " * " << u.coefv[i];
        if (u.face_sideU_ == u.face_sideV_ && u.face_sideU_ != -1)
            f << "\t in Omega_" << u.face_sideU_ + 1;

        f << std::endl;
        return f;
    }
};

template <typeMesh mesh_t> ItemVF<mesh_t> &operator*=(R cc, ItemVF<mesh_t> &F) {
    F.c *= cc;
    return F;
}
template <typeMesh mesh_t> ItemVF<mesh_t> &operator*=(ItemVF<mesh_t> &F, R cc) {
    F.c *= cc;
    return F;
}
template <typeMesh mesh_t> ItemVF<mesh_t> operator-(const ItemVF<mesh_t> &X) {
    ItemVF<mesh_t> F(X);
    F.c *= -1;
    return F;
}

template <typeMesh M> class ListItemVF {

    using mesh_t    = M;
    using fespace_t = GFESpace<mesh_t>;
    using Rd        = typename mesh_t::Rd;
    using item_t    = ItemVF<mesh_t>;
    using this_t    = ListItemVF<mesh_t>;

    static const int D = mesh_t::D;

  public:
    std::vector<item_t> VF;

    bool isRHS_ = true;
    ListItemVF(int l) : VF(l){};
    const item_t &operator()(int i) const { return VF[i]; }
    const item_t &operator[](int i) const { return VF[i]; }
    item_t &operator()(int i) { return VF[i]; }
    item_t &operator[](int i) { return VF[i]; }
    int size() const { return VF.size(); }

    ListItemVF &operator*(const double cc) {
        for (int i = 0; i < VF.size(); ++i)
            (*this)(i).c *= cc;
        return *this;
    }
    ListItemVF &operator+(const ListItemVF &L) {
        int n0 = VF.size();
        int n  = L.size() + this->size();
        this->VF.resize(n);
        for (int i = n0, j = 0; i < n; ++i, ++j)
            (*this)(i) = L(j);
        return *this;
    }
    ListItemVF &operator-(const ListItemVF &L) {
        int n0 = VF.size();
        int n  = L.size() + this->size();
        this->VF.resize(n);
        for (int i = n0, j = 0; i < n; ++i, ++j)
            (*this)(i) = -L(j);

        return *this;
    }
    ListItemVF &operator-() {
        for (int i = 0; i < this->size(); ++i)
            (*this)(i) = -(*this)(i);
        return *this;
    }
    ListItemVF &operator+() { return *this; }

    void reduce() {
        // get size new list
        int l = VF.size();
        KN<int> s(VF.size(), -1);
        KN<int> s2k(VF.size(), -1);

        for (int i = 0; i < VF.size(); ++i) {
            if (VF[i].c == 0) {
                l -= 1;
                continue;
            }
            if (s(i) != -1)
                continue;
            for (int j = i + 1; j < VF.size(); ++j) {
                if (VF.at(i).operator==(VF.at(j))) {
                    s(j) = i;
                    l -= 1;
                }
            }
        }
        if (l == VF.size())
            return;
        std::vector<item_t> u(l);
        int k = 0;

        for (int i = 0; i < VF.size(); ++i) {
            if (VF[i].c == 0) {
                continue;
            }
            if (s(i) == -1) {
                s2k[i] = k;
                u[k++] = VF[i];
            } else {
                u[s2k(s(i))].c += VF[i].c;
            }
        }
    }

    int get_lastOp() const {
        int n = 0;
        for (int i = 0; i < this->size(); ++i)
            n = Max(n, getLastop(VF[i].du, VF[i].dv));
        return n;
    }

    bool isRHS() const { return isRHS_; }
    const fespace_t &get_spaceU(int i) const {
        assert(VF[i].fespaceU || isRHS_);
        return (VF[i].fespaceU) ? *VF[i].fespaceU : *VF[i].fespaceV;
    }
    const fespace_t &get_spaceV(int i) const {
        assert(VF[i].fespaceV);
        return *VF[i].fespaceV;
    }

    friend std::ostream &operator<<(std::ostream &f, const this_t &u) {
        f << u.VF;
        return f;
    }
};

template <typeMesh mesh_t> ListItemVF<mesh_t> jump(const ListItemVF<mesh_t> &L) {
    int n0 = L.VF.size();
    int n  = 2 * n0;

    ListItemVF<mesh_t> item(n);
    item.isRHS_ = L.isRHS_;
    for (int i = 0; i < n0; ++i) {
        item(i)             = L(i);
        item(i).face_sideU_ = 0;
        item(i).face_sideV_ = 0;
    }
    for (int i = n0, j = 0; i < n; ++i, ++j) {
        item(i) = L(j);
        item(i).c *= -1;
        item(i).face_sideU_ = 1;
        item(i).face_sideV_ = 1;
    }
    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> average(const ListItemVF<mesh_t> &L, int v1, int v2) {
    int n0 = L.VF.size();
    int n  = 2 * n0;

    ListItemVF<mesh_t> item(n);
    item.isRHS_ = L.isRHS_;
    for (int i = 0; i < n0; ++i) {
        item(i) = L(i);
        item(i).c *= v1;
        item(i).face_sideU_ = 0;
        item(i).face_sideV_ = 0;
    }
    for (int i = n0, j = 0; i < n; ++i, ++j) {
        item(i) = L(j);
        item(i).c *= v2;
        item(i).face_sideU_ = 1;
        item(i).face_sideV_ = 1;
    }
    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> operator*(R cc, const ListItemVF<mesh_t> &F) {
    ListItemVF<mesh_t> multc(F.size());
    multc.isRHS_ = F.isRHS_;
    int i        = 0;
    for (const auto &item : F.VF) {
        multc.VF[i] = item;
        multc.VF[i++].c *= cc;
    }
    return multc;
}

template <typeMesh mesh_t>
ListItemVF<mesh_t> operator,(const TestFunction<mesh_t> &uu, const TestFunction<mesh_t> &vv) {
    assert(uu.nbRow() == vv.nbRow());
    assert(uu.nbCol() == vv.nbCol()); // && nbCol()==1);
    int D = mesh_t::D;
    int l = 0;
    for (int i = 0; i < uu.nbRow(); ++i) {
        for (int j = 0; j < uu.nbCol(); ++j) {
            l += uu(i, j).size() * vv(i, j).size();
        }
    }

    ListItemVF<mesh_t> item(l);
    item.isRHS_ = false;
    int k = 0, kloc = 0;
    for (int i = 0; i < uu.nbRow(); ++i) {
        for (int j = 0; j < uu.nbCol(); ++j) {
            for (int ui = 0; ui < uu(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &u(uu(i, j).getItem(ui));
                for (int vi = 0; vi < vv(i, j).size(); ++vi) {
                    const ItemTestFunction<mesh_t> &v(vv(i, j).getItem(vi));
                    item(k) = ItemVF<mesh_t>(u, v);
                    k++;
                }
            }
        }
    }
    item.reduce();
    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> operator,(const R c, const TestFunction<mesh_t> &F) {
    int l = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            l += F(i, j).size();
        }
    }

    ListItemVF<mesh_t> item(l);
    int k = 0, kloc = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
                item(k)             = ItemVF<mesh_t>(v.c * c, v.cu, 0, v.cu, v.du, std::vector<int>(), v.ar_nu);
                item(k).face_sideU_ = v.face_side_;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_;
                item(k).coefv       = v.coefu;
                item(k).dtu         = -1;
                item(k).dtv         = v.dtu;
                item(k).expru       = v.expru;
                item(k).fespaceV    = v.fespace;
                k++;
            }
        }
    }
    item.reduce();

    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> operator,(const Rnm &c, const TestFunction<mesh_t> &F) {
    assert(c.N() == F.nbRow() && c.M() == F.nbCol());
    int l = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            l += F(i, j).size();
        }
    }

    ListItemVF<mesh_t> item(l);
    int k = 0, kloc = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
                item(k)             = ItemVF<mesh_t>(v.c * c(i, j), v.cu, 0, v.cu, v.du, std::vector<int>(), v.ar_nu);
                item(k).face_sideU_ = v.face_side_;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_, item(k).coefv = v.coefu;
                item(k).dtu      = -1;
                item(k).dtv      = v.dtu;
                item(k).exprv    = v.expru;
                item(k).fespaceV = v.fespace;

                k++;
            }
        }
    }
    item.reduce();

    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> operator,(const Projection &c, const TestFunction<mesh_t> &F) {
    int l = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            l += F(i, j).size() * (1 + (i == j));
        }
    }

    ListItemVF<mesh_t> item(l);
    int k = 0, kloc = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));

                if (i == j) {
                    item(k)             = ItemVF<mesh_t>(v.c, v.cu, 0, v.cu, v.du, 0, v.ar_nu);
                    item(k).face_sideU_ = v.face_side_;
                    item(k).face_sideV_ = v.face_side_;
                    item(k).domainU_id_ = v.domain_id_;
                    item(k).domainV_id_ = v.domain_id_, item(k).coefv = v.coefu;
                    item(k).dtu      = -1;
                    item(k).dtv      = v.dtu;
                    item(k).exprv    = v.expru;
                    item(k).fespaceV = v.fespace;
                    k++;
                }
                item(k)             = ItemVF<mesh_t>(v.c * (-1), v.cu, 0, v.cu, v.du, c(i, j), v.ar_nu);
                item(k)             = ItemVF<mesh_t>(v.c * (-1), v.cu, 0, v.cu, v.du, c(i, j), v.ar_nu);
                item(k).face_sideU_ = v.face_side_;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_, item(k).coefv = v.coefu;
                item(k).dtu = -1;
                item(k).dtv = v.dtu;

                item(k).exprv    = v.expru;
                item(k).fespaceV = v.fespace;

                k++;
            }
        }
    }
    item.reduce();

    return item;
}

template <typeMesh mesh_t> ListItemVF<mesh_t> operator,(const ExpressionAverage &fh, const TestFunction<mesh_t> &F) {
    int l = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            l += F(i, j).size();
        }
    }
    l *= 2;

    ListItemVF<mesh_t> item(l);
    int k = 0, kloc = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
                item(k)             = ItemVF<mesh_t>(v.c * fh.k1, 0, -1, v.cu, v.du, std::vector<int>(), v.ar_nu);
                item(k).face_sideU_ = 0;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_;
                item(k).coefv       = v.coefu;
                item(k).dtu         = 0;
                item(k).dtv         = v.dtu;
                item(k).expru       = fh.fun1;
                item(k).exprv       = v.expru;
                item(k).fespaceV    = v.fespace;

                k++;
            }
        }
    }

    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
                item(k)             = ItemVF<mesh_t>(v.c * fh.k2, 0, -1, v.cu, v.du, std::vector<int>(), v.ar_nu);
                item(k).face_sideU_ = 1;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_;
                item(k).coefv       = v.coefu;
                item(k).dtu         = 0;
                item(k).dtv         = v.dtu;
                item(k).expru       = fh.fun1;
                item(k).exprv       = v.expru;
                item(k).fespaceV    = v.fespace;

                k++;
            }
        }
    }

    item.reduce();

    return item;
}

template <typeMesh mesh_t, typename Expr>
ListItemVF<mesh_t> operator,(const std::vector<std::shared_ptr<Expr>> &fh, const TestFunction<mesh_t> &F) {
    if (F.nbRow() != fh.size()) {
        std::cout << "size expression \t" << fh.size() << std::endl;
        std::cout << "size test function \t" << F.nbRow() << std::endl;
    }
    assert(F.nbRow() == fh.size());
    int l = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            l += F(i, j).size();
        }
    }

    ListItemVF<mesh_t> item(l);
    int k = 0, kloc = 0;
    for (int i = 0; i < F.nbRow(); ++i) {
        for (int j = 0; j < F.nbCol(); ++j) {
            for (int ui = 0; ui < F(i, j).size(); ++ui) {
                const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
                item(k)             = ItemVF<mesh_t>(v.c, i, -1, v.cu, v.du, std::vector<int>(), v.ar_nu);
                item(k).face_sideU_ = v.face_side_;
                item(k).face_sideV_ = v.face_side_;
                item(k).domainU_id_ = v.domain_id_;
                item(k).domainV_id_ = v.domain_id_;
                item(k).coefv       = v.coefu;
                item(k).dtu         = 0;
                item(k).dtv         = v.dtu;
                item(k).expru       = fh[i];
                item(k).exprv       = v.expru;
                item(k).fespaceV    = v.fespace;

                k++;
            }
        }
    }

    item.reduce();

    return item;
}

template <typeMesh mesh_t>
ListItemVF<mesh_t> innerProduct(const std::shared_ptr<ExpressionVirtual> &fh, const TestFunction<mesh_t> &F) {

    std::vector<std::shared_ptr<ExpressionVirtual>> l;
    l.push_back(fh);
    return (l, F);
}

template <typeMesh mesh_t> ListItemVF<mesh_t> innerProduct(const ExpressionAverage &fh, const TestFunction<mesh_t> &F) {
    return (fh, F);
}

template <typeMesh mesh_t> ListItemVF<mesh_t> innerProduct(double c, const TestFunction<mesh_t> &F) {
    return operator,(c, F);
}

template <typeMesh mesh_t, typename Expr>
ListItemVF<mesh_t> innerProduct(const std::vector<std::shared_ptr<Expr>> &fh, const TestFunction<mesh_t> &F) {
    return operator,(fh, F);
}

template <typeMesh mesh_t>
ListItemVF<mesh_t> innerProduct(const TestFunction<mesh_t> &F1, const TestFunction<mesh_t> &F2) {
    return (F1, F2);
}

template <typeMesh mesh_t>
ListItemVF<mesh_t> contractProduct(const TestFunction<mesh_t> &F1, const TestFunction<mesh_t> &F2) {
    return (F1, F2);
}

template <typeMesh mesh_t> ListItemVF<mesh_t> contractProduct(const Rnm &F1, const TestFunction<mesh_t> &F2) {
    return (F1, F2);
}

template <typeMesh mesh_t> ListItemVF<mesh_t> contractProduct(const Projection &F1, const TestFunction<mesh_t> &F2) {
    return (F1, F2);
}

#endif

// ListItemVF<mesh_t>
// operator,(std::list<ExpressionFunFEM<typename MeshType<mesh_t>::mesh_t> *> fh,
//           const TestFunction<mesh_t> &F) {
//    if (F.nbRow() != fh.size()) {
//       std::cout << "size expression \t" << fh.size() << std::endl;
//       std::cout << "size test function \t" << F.nbRow() << std::endl;
//    }
//    assert(F.nbRow() == fh.size());
//    int l = 0;
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          l += F(i, j).size();
//       }
//    }

//    ListItemVF<mesh_t> item(l);
//    int k = 0, kloc = 0;
//    auto it = fh.begin();
//    for (int i = 0; i < F.nbRow(); ++i, ++it) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < F(i, j).size(); ++ui) {
//             const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
//             item(k) =
//                 ItemVF<mesh_t>(v.c, i, -1, v.cu, v.du, std::vector<int>(),
//                 v.ar_nu);
//             item(k).face_sideU_ = v.face_side_;
//             item(k).face_sideV_ = v.face_side_;
//             item(k).domainU_id_ = v.domain_id_;
//             item(k).domainV_id_ = v.domain_id_, item(k).coefv = v.coefu;
//             item(k).dtu      = 0;
//             item(k).dtv      = v.dtu;
//             item(k).expru    = *it;
//             item(k).exprv    = v.expru;
//             item(k).fespaceV = v.fespace;

//             k++;
//          }
//       }
//    }

//    item.reduce();

//    return item;
// }

// template <typeMesh mesh_t>
// ListItemVF<mesh_t> innerProduct(const ExpressionVirtual &fh,
//                            const TestFunction<mesh_t> &F) {
//    return operator,(fh, F);
// }

// template <typeMesh mesh_t>
// ListItemVF<mesh_t> innerProduct(const ExpressionAverage& fh, const
// TestFunction<mesh_t>& F) {
//  return (fh,F);
// }

// template <typeMesh mesh_t>
// ListItemVF<mesh_t> operator,(const ExpressionVirtual &fh, const TestFunction<mesh_t>
// &F) {
//    int l = 0;
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          l += F(i, j).size();
//       }
//    }

//    ListItemVF<mesh_t> item(l);
//    int k = 0, kloc = 0;
//    for (int i = 0; i < F.nbRow(); ++i) {
//       for (int j = 0; j < F.nbCol(); ++j) {
//          for (int ui = 0; ui < F(i, j).size(); ++ui) {
//             const ItemTestFunction<mesh_t> &v(F(i, j).getItem(ui));
//             item(k) =
//                 ItemVF<mesh_t>(v.c, 0, -1, v.cu, v.du, std::vector<int>(),
//                 v.ar_nu);
//             item(k).face_sideU_ = v.face_side_;
//             item(k).face_sideV_ = v.face_side_;
//             item(k).domainU_id_ = v.domain_id_;
//             item(k).domainV_id_ = v.domain_id_, item(k).coefv = v.coefu;
//             item(k).dtu      = 0;
//             item(k).dtv      = v.dtu;
//             item(k).expru    = &fh;
//             item(k).exprv    = v.expru;
//             item(k).fespaceV = v.fespace;

//             k++;
//          }
//       }
//    }

//    item.reduce();

//    return item;
// }
