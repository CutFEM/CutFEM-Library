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
#ifndef _EXPRESSION_HPP
#define _EXPRESSION_HPP
#include "../parallel/cfomp.hpp"
#include "../common/memoryPool.hpp"
#include "FESpace.hpp"
#include <cmath>
#include <memory>
#include <list>
#include "../algoim/uvector.hpp"
#include "../algoim/interval.hpp"



struct Normal_Component {
    virtual int component() const = 0;
};
struct Normal_Component_X : public Normal_Component {
    virtual int component() const { return 0; }
};
struct Normal_Component_Y : public Normal_Component {
    virtual int component() const { return 1; }
};
struct Normal_Component_Z : public Normal_Component {
    virtual int component() const { return 2; }
};

struct BaseVector {
    virtual int operator[](int i) const = 0;
    virtual int cst(int i) const        = 0;
};
struct Normal : public BaseVector {
    Normal_Component_X x;
    Normal_Component_Y y;
    Normal_Component_Z z;
    int operator[](int i) const { return i; }
    int cst(int i) const { return 1; }
};
struct Conormal : public BaseVector {
    Normal_Component_X x;
    Normal_Component_Y y;
    Normal_Component_Z z;
    int operator[](int i) const { return i; }
    int cst(int i) const { return 1; }
};
struct Tangent : public BaseVector {
    int operator[](int i) const { return (i == 0); }
    int cst(int i) const { return (i == 0) ? -1 : 1; }
};
struct Projection {
    Normal normal;
    KN<int> operator()(int i, int j) const {
        KN<int> ar(2);
        ar(0) = i;
        ar(1) = j;
        return ar;
    }
};

class ExpressionVirtual;
template <typename M> class ExpressionFunFEM;
class ParameterCutFEM;

/**
 * @brief Base class for all expressions
 *
 */
class FunFEMVirtual {
  public:
    // data used if function allocate memory
    std::vector<double> data_;

    // view on data_, either passed from outside or allocated by the function
    std::span<double> v;

    // memory pool for basis functions
    std::unique_ptr<MemoryPool> pool_databf;

    FunFEMVirtual() : v(data_) {}
    FunFEMVirtual(int df) : data_(df, 0.), v(data_) {}
    FunFEMVirtual(std::span<double> u) : v(u) {}

    virtual double eval(const int k, const R *x, int cu = 0, int op = 0) const {
        assert(0);
        return 0.;
    };
    virtual double eval(const int k, const R *x, const R t, int cu, int op, int opt) const {
        assert(0);
        return 0.;
    };
    virtual double evalOnBackMesh(const int k, int dom, const R *x, int cu = 0, int op = 0) const {
        assert(0);
        return 0.;
    };
    virtual double evalOnBackMesh(const int k, int dom, const R *x, const R t, int cu, int op, int opt) const {
        assert(0);
        return 0.;
    };
    virtual int idxElementFromBackMesh(int, int = 0) const {
        assert(0);
        return 0.;
    };

    std::span<double> array() const { return v; }

    // FunFEMVirtual& operator+=(const FunFEMVirtual &other) {
    //     std::vector<double> u_dof_n(data_uh.begin() + n * Wh.get_nb_dof(),
    //     data_uh.begin() + (n + 1) * Wh.get_nb_dof());
    //     std::transform(sol.begin(), sol.end(), u_dof_n.begin(), sol.begin(),
    //     std::plus<double>()); // sum up all dofs

    // }
};

/**
 * @brief FunFEM represents functions defined on a finite element space
 * @note
 * @tparam M Mesh
 */
template <typename M> class FunFEM : public FunFEMVirtual {
  public:
    typedef M Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;

  public:
    FESpace const *Vh  = nullptr;
    TimeSlab const *In = nullptr;

  public:
    FunFEM() : FunFEMVirtual() {}

    // START OF CHATGPT ADDITION
    // // Deep copy
    // FunFEM(const FunFEM& other)
    //     : FunFEMVirtual(), Vh(other.Vh), In(other.In)
    // {
    //     data_.assign(other.v.begin(), other.v.end());                 // allocate & copy coefficients
    //     v = std::span<double>(data_.data(), data_.size());            // bind to *this* storage

    //     if (Vh) {
    //         std::size_t databf_size = (*Vh)[0].NbDoF() * Vh->N * 10;
    //         std::size_t n_chunk     = cutfem_get_max_threads();
    //         pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    //     }
    // }

    // // Deep-copy assignment: same idea
    // FunFEM& operator=(const FunFEM& other) {
    //     if (this == &other) return *this;
    //     Vh = other.Vh;
    //     In = other.In;

    //     data_.assign(other.v.begin(), other.v.end());                 // copy from the view, not data_
    //     v = std::span<double>(data_.data(), data_.size());            // rebind to our storage

    //     if (Vh) {
    //         std::size_t databf_size = (*Vh)[0].NbDoF() * Vh->N * 10;
    //         std::size_t n_chunk     = cutfem_get_max_threads();
    //         pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    //     } else {
    //         pool_databf.reset();
    //     }
    //     return *this;
    // }

    // // Moves are fine to default (unique_ptr is movable)
    // FunFEM(FunFEM&&) noexcept = default;
    // FunFEM& operator=(FunFEM&&) noexcept = default;
    FunFEM(FunFEM&&) noexcept = delete;
    FunFEM& operator=(FunFEM&&) noexcept = delete;

    ~FunFEM() = default;

    // END OF CHATGPT ADDITION

    explicit FunFEM(const FESpace &vh) : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {

        // ndof_per_componenent * nb_component * number_of_possible_operator (op_id, op_dx etc)
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    FunFEM(const FESpace &vh, const TimeSlab &in) : FunFEMVirtual(vh.NbDoF() * in.NbDoF()), Vh(&vh), In(&in) {
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    template <Vector vector_t>
    FunFEM(const FESpace &vh, const TimeSlab &in, vector_t &u) : FunFEMVirtual(u), Vh(&vh), In(&in) {
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    template <Vector vector_t> FunFEM(const FESpace &vh, vector_t &u) : FunFEMVirtual(u), Vh(&vh) {
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    template <typename fct_t>
        requires FunctionLevelSet<fct_t> || FunctionDomain<fct_t> || FunctionScalar<fct_t>
    FunFEM(const FESpace &vh, fct_t f) : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {

        // allocate memory for basis functions
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

        interpolate(*Vh, this->v, f);
    }
    FunFEM(const FESpace &vh, double f) : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {
        // allocate memory for basis functions
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

        std::ranges::fill(this->v, f);
    }

    template <typename fun_t> FunFEM(const FESpace &vh, fun_t f, R tid) : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {

        // allocate memory for basis functions
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

        interpolate(*Vh, this->v, f, tid);
    }

    template <typename fct_t>
        requires FunctionLevelSetTime<fct_t> || FunctionDomainTime<fct_t> || FunctionScalar<fct_t>
    FunFEM(const FESpace &vh, const TimeSlab &in, fct_t f) : FunFEMVirtual(vh.NbDoF() * in.NbDoF()), Vh(&vh), In(&in) {

        // allocate memory for basis functions
        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

        interpolate(*Vh, *In, this->v, f);
    }
    FunFEM(const FESpace &vh, const ExpressionVirtual &fh);
    FunFEM(const FESpace &vh, const ExpressionVirtual &fh1, const ExpressionVirtual &fh2);

    void init(const FESpace &vh) {
        Vh = &vh;
        data_.resize(Vh->NbDoF(), 0.);
        v = data_;

        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }
    void init(const KN_<R> &a) {
        assert(v.size() == a.size());
        for (int i = 0; i < v.size(); ++i)
            v[i] = a(i);

        assert(pool_databf);
    }

    template <typename fct_t>
        requires FunctionLevelSet<fct_t> || FunctionLevelSetTime<fct_t> || FunctionDomain<fct_t> ||
                 FunctionScalar<fct_t>
    void init(const FESpace &vh, fct_t f) {
        Vh = &vh;
        data_.resize(Vh->NbDoF(), 0.);
        v = data_;
        interpolate(*Vh, v, f);

        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    void init(const FESpace &vh, R (*f)(double *, int, R), R tid) {
        Vh = &vh;
        data_.resize(Vh->NbDoF(), 0.);
        v = data_;
        interpolate(*Vh, v, f, tid);

        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    void init(const FESpace &vh, R (*f)(R2, int, R), R tid) {

        Vh = &vh;
        data_.resize(Vh->NbDoF(), 0.);
        v = data_;
        interpolate(*Vh, v, f, tid);

        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    template <typename fct_t>
        requires FunctionLevelSet<fct_t> || FunctionDomain<fct_t> || FunctionScalar<fct_t>
    void init(const FESpace &vh, const TimeSlab &in, fct_t f) {

        Vh = &vh;
        In = &in;
        data_.resize(Vh->NbDoF() * In->NbDoF(), 0.);
        v = data_;
        interpolate(*Vh, *In, v, f);

        std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
        std::size_t n_chunk     = cutfem_get_max_threads();
        pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);
    }

    // friend void swap(FunFEM &f, FunFEM &g) {
    //     assert(g.v.size() == f.v.size());
    //     double *temp = g.data_;
    //     g.data_      = f.data_;
    //     f.data_      = temp;
    //     g.v.set(g.data_, g.Vh->NbDoF());
    //     f.v.set(f.data_, f.Vh->NbDoF());
    // }

    // friend void swap(FunFEM& f, FunFEM& g) noexcept {
    //     // Optional consistency checks (keep if these invariants are required)
    //     assert(f.Vh && g.Vh);
    //     assert(f.data_.size() == static_cast<size_t>(f.Vh->NbDoF()));
    //     assert(g.data_.size() == static_cast<size_t>(g.Vh->NbDoF()));

    //     std::swap(f.data_, g.data_);   // swap the owning buffers

    //     // Rebind spans to the new underlying storage
    //     f.v = std::span<double>(f.data_.data(), f.data_.size());
    //     g.v = std::span<double>(g.data_.data(), g.data_.size());
    // }

    double &operator()(int i) { return v[i]; }
    double operator()(int i) const { return v[i]; }
    operator Rn() const { return Rn(v); }

    double eval(const int k, const R *x, int cu = 0, int op = 0) const;
    double eval(const int k, const R *x, const R t, int cu = 0, int op = 0, int opt = 0) const;
    // double operator()(const int k, const R *x) { return eval(k, x, 0, op_id); }
    void eval(R *u, const int k) const;

    double evalOnBackMesh(const int kb, int dom, const R *x, int cu, int op) const;
    double evalOnBackMesh(const int kb, int dom, const R *x, const R t, int cu, int op, int opt) const;

    int size() const { return Vh->NbDoF() * ((In) ? In->NbDoF() : 1); }
    int size(int k) const { return (*Vh)[k].NbDoF(); }
    void print() const;
    int idxElementFromBackMesh(int kb, int dd = 0) const { return Vh->idxElementFromBackMesh(kb, dd); }
    std::vector<int> getAllDomainId(int k) const { return Vh->getAllDomainId(k); };

    const FESpace &getSpace() const { return *Vh; }
    BasisFctType getBasisFctType() const { return Vh->basisFctType; }
    int getPolynomialOrder() const { return Vh->polynomialOrder; }
    std::shared_ptr<ExpressionFunFEM<M>> expr(int i0 = 0) const;
    std::vector<std::shared_ptr<ExpressionFunFEM<M>>> exprList(int n = -1) const;
    std::vector<std::shared_ptr<ExpressionFunFEM<M>>> exprList(int n, int i0) const;

    FunFEM<M> &operator+=(const FunFEM<M> &other);
    FunFEM<M> &operator-=(const FunFEM<M> &other);
    FunFEM<M> &operator*=(const double c);
    FunFEM<M> &operator/=(const double c);


    // Algoim additions
    using Itv = algoim::Interval<Rd::d>;

    mutable int algoim_k_ = -1;
    mutable double t_algoim_ = -1.0;
    
    void setTime(double t) const { t_algoim_ = t; }
    void setElementFromBackMesh(int kb, int dom = 0) const {
        algoim_k_ = Vh->idxElementFromBackMesh(kb, dom);
    }

    // Evaluation in Rd points
    // double operator()(const algoim::uvector<double, Rd::d>& x) const {
    //     return eval(algoim_k_, x.data(), 0, op_id);
    // }
    // algoim::uvector<double, Rd::d> grad(const algoim::uvector<double, Rd::d>& x) const {
    //     algoim::uvector<double, Rd::d> g;
    //     for (int i = 0; i < Rd::d; ++i)
    //         g(i) = eval(algoim_k_, x.data(), 0, op_dx + i);  // op_dx=1, op_dy=2, op_dz=3
    //     return g;
    // }
    double operator()(const Rd& x) const {
        if (In && t_algoim_ >= 0.0)
            return eval(algoim_k_, x, t_algoim_, 0, op_id, op_id);
        return eval(algoim_k_, x, 0, op_id);
    }
    double operator()(const algoim::uvector<double, Rd::d>& x) const {
        if (In && t_algoim_ >= 0.0)
            return eval(algoim_k_, x.data(), t_algoim_, 0, op_id, op_id);
        return eval(algoim_k_, x.data(), 0, op_id);
    }
    algoim::uvector<double, Rd::d> grad(const algoim::uvector<double, Rd::d>& x) const {
        algoim::uvector<double, Rd::d> g;
        for (int i = 0; i < Rd::d; ++i)
            g(i) = (In && t_algoim_ >= 0.0)
                ? eval(algoim_k_, x.data(), t_algoim_, 0, op_dx + i, op_id)
                : eval(algoim_k_, x.data(), 0, op_dx + i);
        return g;
    }

    
    // Evaluation in interval arithmetic
    Itv operator()(const algoim::uvector<Itv, Rd::d>& x) const;   // implemented in .tpp
    algoim::uvector<Itv, Rd::d> grad(const algoim::uvector<Itv, Rd::d>& x) const;   // .tpp

    Rd normal(Rd x) const {
        algoim::uvector<double, Rd::d> xv;
        for (int i = 0; i < Rd::d; ++i) xv(i) = x[i];
        auto g = grad(xv);
        double nrm = 0; for (int i = 0; i < Rd::d; ++i) nrm += g(i)*g(i);
        nrm = std::sqrt(nrm);
        Rd n; for (int i = 0; i < Rd::d; ++i) n[i] = g(i)/nrm;
        return n;
    }

    // ~FunFEM() {}

};

class ExpressionVirtual {
  public:
    int cu, op, opt;
    KN<int> ar_normal;
    int domain = -1;

    ExpressionVirtual() : cu(0), op(0), opt(0) {}
    ExpressionVirtual(int cc, int opp) : cu(cc), op(opp), opt(0) {}
    ExpressionVirtual(int cc, int opp, int oppt) : cu(cc), op(opp), opt(oppt) {}
    ExpressionVirtual(int cc, int opp, int oppt, int dom) : cu(cc), op(opp), opt(oppt), domain(dom) {}
    virtual R operator()(long i) const                                                  = 0;
    virtual R eval(const int k, const R *x, const R *normal = nullptr) const            = 0;
    virtual R eval(const int k, const R *x, const R t, const R *normal = nullptr) const = 0;

  public:
    virtual R evalOnBackMesh(const int k, int dom, const R *x, const R *normal = nullptr) const            = 0;
    virtual R evalOnBackMesh(const int k, int dom, const R *x, const R t, const R *normal = nullptr) const = 0;
    virtual int idxElementFromBackMesh(int kb, int dd = 0) const                                           = 0;
    virtual std::vector<int> getAllDomainId(int k) const { assert(0); return std::vector<int>();};

  public:
    R GevalOnBackMesh(const int k, int dom, const R *x, const R *normal) const {
        int theDomain = (domain == -1) ? dom : domain;
        return evalOnBackMesh(k, theDomain, x, normal);
    }
    R GevalOnBackMesh(const int k, int dom, const R *x, const R t, const R *normal) const {
        int theDomain = (domain == -1) ? dom : domain;
        return evalOnBackMesh(k, theDomain, x, t, normal);
    }

    virtual int size() const { assert(0); return 0;};

    ExpressionVirtual &operator=(const ExpressionVirtual &L) {
        cu  = L.cu;
        op  = L.op;
        opt = L.opt;
        ar_normal.init(L.ar_normal);
        return *this;
    }

    double computeNormal(const R *normal) const {
        if (normal == nullptr)
            return 1.;
        R val = 1;
        for (int i = 0; i < ar_normal.size(); ++i)
            val *= normal[ar_normal(i)];
        return val;
    }
    void addNormal(int i) {
        int l = ar_normal.size();
        ar_normal.resize(l + 1);
        ar_normal(l) = i;
    }

    virtual void whoAmI() const { std::cout << " I am virtual class Expression" << std::endl; }
    virtual ~ExpressionVirtual() {}
};

class DefaultExpression : public ExpressionVirtual {

    virtual R operator()(long i) const override {
        assert(0);
        return 0;
    };
    virtual R eval(const int k, const R *x, const R *normal) const override {
        assert(0);
        return 0;
    };
    virtual R eval(const int k, const R *x, const R t, const R *normal) const override {
        assert(0);
        return 0;
    };

  public:
    virtual R evalOnBackMesh(const int k, int dom, const R *x, const R *normal) const override {
        assert(0);
        return 0;
    };

    virtual R evalOnBackMesh(const int k, int dom, const R *x, const R t, const R *normal) const override {
        assert(0);
        return 0;
    };

    virtual int idxElementFromBackMesh(int kb, int dd = 0) const override {
        assert(0);
        return 0;
    };
};

template <typename M> class ExpressionFunFEM : public ExpressionVirtual {

  protected:
  public:
    const FunFEM<M> &fun;

  public:
    ExpressionFunFEM(const FunFEM<M> &fh, int cc, int opp) : ExpressionVirtual(cc, opp), fun(fh) {}
    ExpressionFunFEM(const FunFEM<M> &fh, int cc, int opp, int oppt) : ExpressionVirtual(cc, opp, oppt), fun(fh) {}
    ExpressionFunFEM(const FunFEM<M> &fh, int cc, int opp, int oppt, int dom)
        : ExpressionVirtual(cc, opp, oppt, dom), fun(fh) {}

    R operator()(long i) const { return fun(i); }
    virtual int size() const { return fun.size(); }

    void whoAmI() const { std::cout << " I am class ExpressionFunFEM" << std::endl; }

    R eval(const int k, const R *x, const R *normal) const { return fun.eval(k, x, cu, op) * computeNormal(normal); }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        return fun.eval(k, x, t, cu, op, opt) * computeNormal(normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        return fun.evalOnBackMesh(k, dom, x, cu, op) * computeNormal(normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        return fun.evalOnBackMesh(k, dom, x, t, cu, op, opt) * computeNormal(normal);
    }

  public:
    const GFESpace<M> &getSpace() const { return *fun.Vh; }

    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }

    ExpressionFunFEM(const ExpressionFunFEM &L) : ExpressionVirtual(L.cu, L.op, L.opt, L.domain), fun(L.fun) {
        ar_normal.init(L.ar_normal);
    }

    ExpressionFunFEM &operator=(const ExpressionFunFEM &L) {
        cu  = L.cu;
        op  = L.op;
        opt = L.opt;
        ar_normal.init(L.ar_normal);
        fun(L.fun);
        domain = L.domain;
        return *this;
    }

    ExpressionFunFEM operator*(const Normal &n) {
        ExpressionFunFEM ff(*this);
        ff.addNormal(cu);
        return ff;
    }

    // friend ExpressionFunFEM dx<M>(const ExpressionFunFEM<M> &u);
    // friend ExpressionFunFEM dy<M>(const ExpressionFunFEM<M> &u);
    // friend ExpressionFunFEM dz<M>(const ExpressionFunFEM<M> &u);
    // friend ExpressionFunFEM dt<M>(const ExpressionFunFEM<M> &u);
};

template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dx(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dx, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dy(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dy, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dz(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dz, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dt(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, u->op, op_dx, u->domain);
}

/* Sebastian's implementation, is this correct?
It generates an assertion error */
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dxx(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dxx, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dyy(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dyy, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dzz(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dzz, u->opt, u->domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dxy(const std::shared_ptr<ExpressionFunFEM<M>> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u->fun, u->cu, op_dxy, u->opt, u->domain);
}

template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dx(const ExpressionFunFEM<M> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u.fun, u.cu, op_dx, u.opt, u.domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dy(const ExpressionFunFEM<M> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u.fun, u.cu, op_dy, u.opt, u.domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dz(const ExpressionFunFEM<M> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u.fun, u.cu, op_dz, u.opt, u.domain);
}
template <typename M> std::shared_ptr<ExpressionFunFEM<M>> dt(const ExpressionFunFEM<M> &u) {
    return std::make_shared<ExpressionFunFEM<M>>(u.fun, u.cu, u.op, op_dx, u.domain);
}

template<typename M>
std::shared_ptr<ExpressionVirtual> dot(
    const std::vector<std::shared_ptr<ExpressionFunFEM<M>>>& u,
    const std::vector<std::shared_ptr<ExpressionFunFEM<M>>>& v)
{
    assert(u.size() == v.size() && "Vectors must have the same size for scalar product.");

    std::shared_ptr<ExpressionVirtual> result = u[0] * v[0];
    for (std::size_t i = 1; i < u.size(); ++i) {
        result = result + (u[i] * v[i]);
    }
    return result;
}

// Forward declarations for generic differential operators
class ExpressionSum;
class ExpressionProduct;
std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionSum> &u);
std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionSum> &u);
std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionSum> &u);
std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionSum> &u);
std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionProduct> &u);
std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionProduct> &u);
std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionProduct> &u);
std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionProduct> &u);

// Generic differential operators that dispatch based on runtime type
inline std::shared_ptr<ExpressionVirtual> dx(const std::shared_ptr<ExpressionVirtual> &u) {
    // Try casting to ExpressionFunFEM (need to try different mesh types)
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh3>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshQuad2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshHexa>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    // Try casting to ExpressionSum
    if (auto ptr = std::dynamic_pointer_cast<ExpressionSum>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    // Try casting to ExpressionProduct
    if (auto ptr = std::dynamic_pointer_cast<ExpressionProduct>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dx(ptr));
    }
    // If none match, error
    std::cerr << "Error: dx not implemented for this expression type" << std::endl;
    assert(false && "dx not implemented for this expression type");
    return nullptr;
}

inline std::shared_ptr<ExpressionVirtual> dy(const std::shared_ptr<ExpressionVirtual> &u) {
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh3>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshQuad2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshHexa>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionSum>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionProduct>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dy(ptr));
    }
    std::cerr << "Error: dy not implemented for this expression type" << std::endl;
    assert(false && "dy not implemented for this expression type");
    return nullptr;
}

inline std::shared_ptr<ExpressionVirtual> dz(const std::shared_ptr<ExpressionVirtual> &u) {
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh3>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshQuad2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshHexa>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionSum>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionProduct>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dz(ptr));
    }
    std::cerr << "Error: dz not implemented for this expression type" << std::endl;
    assert(false && "dz not implemented for this expression type");
    return nullptr;
}

inline std::shared_ptr<ExpressionVirtual> dt(const std::shared_ptr<ExpressionVirtual> &u) {
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<Mesh3>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshQuad2>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionFunFEM<MeshHexa>>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionSum>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    if (auto ptr = std::dynamic_pointer_cast<ExpressionProduct>(u)) {
        return std::static_pointer_cast<ExpressionVirtual>(dt(ptr));
    }
    std::cerr << "Error: dt not implemented for this expression type" << std::endl;
    assert(false && "dt not implemented for this expression type");
    return nullptr;
}


class ExpressionMultConst : public ExpressionVirtual {
    const std::shared_ptr<ExpressionVirtual> fun1;
    const double c;
    const bool nx, ny, nz;
    const R2 p;

  public:
    ExpressionMultConst(const std::shared_ptr<ExpressionVirtual> &fh1, const double &cc)
        : fun1(fh1), c(cc), nx(false), ny(false), nz(false), p(R2(1, 1)) {}
    ExpressionMultConst(const std::shared_ptr<ExpressionVirtual> &fh1, const Normal_Component_X &nnx)
        : fun1(fh1), c(1.), nx(true), ny(false), nz(false), p(R2(1, 1)) {}
    ExpressionMultConst(const std::shared_ptr<ExpressionVirtual> &fh1, const Normal_Component_Y &nny)
        : fun1(fh1), c(1.), nx(false), ny(true), nz(false), p(R2(1, 1)) {}
    ExpressionMultConst(const std::shared_ptr<ExpressionVirtual> &fh1, const Normal_Component_Z &nnz)
        : fun1(fh1), c(1.), nx(false), ny(true), nz(true), p(R2(1, 1)) {}
    R operator()(long i) const { return c * ((*fun1)(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        double compN = ((nx) ? normal[0] : 1) * ((ny) ? normal[1] : 1) * ((nz) ? normal[2] : 1);
        return fun1->eval(k, x, normal) * c * compN * p[0];
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        double compN = ((nx) ? normal[0] : 1) * ((ny) ? normal[1] : 1) * ((nz) ? normal[2] : 1);
        return fun1->eval(k, x, t, normal) * c * compN * p[0];
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double compN = ((nx) ? normal[0] : 1) * ((ny) ? normal[1] : 1) * ((nz) ? normal[2] : 1);
        return fun1->evalOnBackMesh(k, dom, x, normal) * c * compN * p[dom == 1];
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double compN = ((nx) ? normal[0] : 1) * ((ny) ? normal[1] : 1) * ((nz) ? normal[2] : 1);
        return fun1->evalOnBackMesh(k, dom, x, t, normal) * c * compN * p[dom == 1];
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
    ~ExpressionMultConst() {}
};
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1, double cc);
std::shared_ptr<ExpressionMultConst> operator*(double cc, const std::shared_ptr<ExpressionVirtual> &f1);
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                               const Normal_Component_X &cc);
std::shared_ptr<ExpressionMultConst> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                               const Normal_Component_Y &cc);
std::shared_ptr<ExpressionMultConst> operator*(const ParameterCutFEM &v, const std::shared_ptr<ExpressionVirtual> &f1);

class ExpressionAbs : public ExpressionVirtual {
    const std::shared_ptr<ExpressionVirtual> fun1;

  public:
    ExpressionAbs(const std::shared_ptr<ExpressionVirtual> &fh1) : fun1(fh1) {}

    R operator()(long i) const { return fabs((*fun1)(i)); }

    R eval(const int k, const R *x, const R *normal) const { return fabs(fun1->eval(k, x, normal)); }
    R eval(const int k, const R *x, const R t, const R *normal) const { return fabs(fun1->eval(k, x, t, normal)); }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        return fabs(fun1->evalOnBackMesh(k, dom, x, normal));
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        return fabs(fun1->evalOnBackMesh(k, dom, x, t, normal));
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
    ~ExpressionAbs() {}
};
std::shared_ptr<ExpressionAbs> fabs(const std::shared_ptr<ExpressionVirtual> &f1);

class ExpressionProduct : public ExpressionVirtual {
    const std::shared_ptr<ExpressionVirtual> fun1;
    const std::shared_ptr<ExpressionVirtual> fun2;

    // Friend declarations for differential operators
    friend std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionProduct> &u);
    friend std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionProduct> &u);
    friend std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionProduct> &u);
    friend std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionProduct> &u);

  public:
    // Template constructor that accepts any types derived from ExpressionVirtual
    template<typename T1, typename T2>
    ExpressionProduct(std::shared_ptr<T1> fh1, std::shared_ptr<T2> fh2)
        : fun1(std::static_pointer_cast<ExpressionVirtual>(fh1)), 
          fun2(std::static_pointer_cast<ExpressionVirtual>(fh2)) {}

    R operator()(long i) const { return (*fun1)(i) * (*fun2)(i); }

    R eval(const int k, const R *x, const R *normal) const {
        return fun1->eval(k, x, normal) * fun2->eval(k, x, normal);
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        return fun1->eval(k, x, t, normal) * fun2->eval(k, x, t, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        return fun1->evalOnBackMesh(k, dom, x, normal) * fun2->evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        return fun1->evalOnBackMesh(k, dom, x, t, normal) * fun2->evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
    ~ExpressionProduct() {}
};
std::shared_ptr<ExpressionVirtual> operator*(const std::shared_ptr<ExpressionVirtual> &f1,
                                             const std::shared_ptr<ExpressionVirtual> &f2);

// Template overload for operator* that accepts derived types
template<typename T1, typename T2,
         typename = std::enable_if_t<std::is_base_of_v<ExpressionVirtual, T1> && 
                                      std::is_base_of_v<ExpressionVirtual, T2>>>
inline std::shared_ptr<ExpressionVirtual> operator*(const std::shared_ptr<T1> &f1,
                                                    const std::shared_ptr<T2> &f2) {
    return std::static_pointer_cast<ExpressionVirtual>(std::make_shared<ExpressionProduct>(f1, f2));
}

// Differential operators for ExpressionProduct (applying product rule: d/dx(f*g) = df/dx * g + f * dg/dx)
// Note: These return ExpressionSum because the product rule produces a sum
inline std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionProduct> &u) {
    auto fun1_dx = dx(u->fun1);
    auto fun2_dx = dx(u->fun2);
    // d/dx(f*g) = df/dx * g + f * dg/dx
    auto term1 = std::make_shared<ExpressionProduct>(fun1_dx, u->fun2);
    auto term2 = std::make_shared<ExpressionProduct>(u->fun1, fun2_dx);
    return std::make_shared<ExpressionSum>(term1, term2);
}

inline std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionProduct> &u) {
    auto fun1_dy = dy(u->fun1);
    auto fun2_dy = dy(u->fun2);
    // d/dy(f*g) = df/dy * g + f * dg/dy
    auto term1 = std::make_shared<ExpressionProduct>(fun1_dy, u->fun2);
    auto term2 = std::make_shared<ExpressionProduct>(u->fun1, fun2_dy);
    return std::make_shared<ExpressionSum>(term1, term2);
}

inline std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionProduct> &u) {
    auto fun1_dz = dz(u->fun1);
    auto fun2_dz = dz(u->fun2);
    // d/dz(f*g) = df/dz * g + f * dg/dz
    auto term1 = std::make_shared<ExpressionProduct>(fun1_dz, u->fun2);
    auto term2 = std::make_shared<ExpressionProduct>(u->fun1, fun2_dz);
    return std::make_shared<ExpressionSum>(term1, term2);
}

inline std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionProduct> &u) {
    auto fun1_dt = dt(u->fun1);
    auto fun2_dt = dt(u->fun2);
    // d/dt(f*g) = df/dt * g + f * dg/dt
    auto term1 = std::make_shared<ExpressionProduct>(fun1_dt, u->fun2);
    auto term2 = std::make_shared<ExpressionProduct>(u->fun1, fun2_dt);
    return std::make_shared<ExpressionSum>(term1, term2);
}

/// @brief Class that compute an function or expression to a power of n
/// @tparam D The dimension
class ExpressionPow : public ExpressionVirtual {

  public:
    typedef std::shared_ptr<ExpressionVirtual> ptr_expr_t;

  private:
    const ptr_expr_t fun1;
    const double n;

  public:
    ExpressionPow(const ptr_expr_t &f1h, const double nn) : fun1(f1h), n(nn) {}

    R operator()(long i) const { return pow((*fun1)(i), n); }

    R eval(const int k, const R *x, const R *normal) const {
        const double val = fun1->eval(k, x, normal);
        return pow(val, n);
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        const double val = fun1->eval(k, x, t, normal);
        return pow(val, n);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        const double val = fun1->evalOnBackMesh(k, dom, x, normal);
        return pow(val, n);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        const double val = fun1->evalOnBackMesh(k, dom, x, t, normal);
        return pow(val, n);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
};

std::shared_ptr<ExpressionPow> pow(const std::shared_ptr<ExpressionVirtual> &f1, const double nn);
std::shared_ptr<ExpressionPow> operator^(const std::shared_ptr<ExpressionVirtual> &f1, const double nn);

std::shared_ptr<ExpressionPow> sqrt(const std::shared_ptr<ExpressionVirtual> &f1);

class ExpressionDivision : public ExpressionVirtual {
    const std::shared_ptr<ExpressionVirtual> fun1;
    const std::shared_ptr<ExpressionVirtual> fun2;

  public:
    ExpressionDivision(const std::shared_ptr<ExpressionVirtual> &fh1, const std::shared_ptr<ExpressionVirtual> &fh2)
        : fun1(fh1), fun2(fh2) {}

    R operator()(long i) const { return (*fun1)(i) / ((*fun2)(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        double v = fun2->eval(k, x, normal);
        assert(fabs(v) > 1e-15);
        return fun1->eval(k, x, normal) / v;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        double v = fun2->eval(k, x, t, normal);
        assert(fabs(v) > 1e-15);
        return fun1->eval(k, x, t, normal) / v;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double v = fun2->evalOnBackMesh(k, dom, x, normal);

        assert(fabs(v) > 1e-15);
        return fun1->evalOnBackMesh(k, dom, x, normal) / v;
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double v = fun2->evalOnBackMesh(k, dom, x, t, normal);
        assert(fabs(v) > 1e-15);
        return fun1->evalOnBackMesh(k, dom, x, t, normal) / v;
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
    ~ExpressionDivision() {}
};
std::shared_ptr<ExpressionDivision> operator/(const std::shared_ptr<ExpressionVirtual> &f1,
                                              const std::shared_ptr<ExpressionVirtual> &f2);

class ExpressionSum : public ExpressionVirtual {
    const std::shared_ptr<ExpressionVirtual> fun1;
    const std::shared_ptr<ExpressionVirtual> fun2;

    // Friend declarations for differential operators
    friend std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionSum> &u);
    friend std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionSum> &u);
    friend std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionSum> &u);
    friend std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionSum> &u);

  public:
    // Template constructor that accepts any types derived from ExpressionVirtual
    template<typename T1, typename T2>
    ExpressionSum(std::shared_ptr<T1> fh1, std::shared_ptr<T2> fh2)
        : fun1(std::static_pointer_cast<ExpressionVirtual>(fh1)), 
          fun2(std::static_pointer_cast<ExpressionVirtual>(fh2)) {}
    R operator()(long i) const { return (*fun1)(i) + (*fun2)(i); }

    R eval(const int k, const R *x, const R *normal) const {
        return fun1->eval(k, x, normal) + fun2->eval(k, x, normal);
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        return fun1->eval(k, x, t, normal) + fun2->eval(k, x, t, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        return fun1->evalOnBackMesh(k, dom, x, normal) + fun2->evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        return fun1->evalOnBackMesh(k, dom, x, t, normal) + fun2->evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
    ~ExpressionSum() {}
};

std::shared_ptr<ExpressionVirtual> operator+(const std::shared_ptr<ExpressionVirtual> &f1,
                                         const std::shared_ptr<ExpressionVirtual> &f2);

std::shared_ptr<ExpressionVirtual> operator-(const std::shared_ptr<ExpressionVirtual> &f1,
                                         const std::shared_ptr<ExpressionVirtual> &f2);

// Template overloads for operator+ and operator- that accept derived types
template<typename T1, typename T2, 
         typename = std::enable_if_t<std::is_base_of_v<ExpressionVirtual, T1> && 
                                      std::is_base_of_v<ExpressionVirtual, T2>>>
inline std::shared_ptr<ExpressionVirtual> operator+(const std::shared_ptr<T1> &f1,
                                                const std::shared_ptr<T2> &f2) {
    return std::static_pointer_cast<ExpressionVirtual>(std::make_shared<ExpressionSum>(f1, f2));
}

template<typename T1, typename T2,
         typename = std::enable_if_t<std::is_base_of_v<ExpressionVirtual, T1> && 
                                      std::is_base_of_v<ExpressionVirtual, T2>>>
inline std::shared_ptr<ExpressionVirtual> operator-(const std::shared_ptr<T1> &f1,
                                                const std::shared_ptr<T2> &f2) {
    auto neg_f2 = std::make_shared<ExpressionMultConst>(f2, -1.0);
    return std::static_pointer_cast<ExpressionVirtual>(std::make_shared<ExpressionSum>(f1, neg_f2));
}

// Differential operators for ExpressionSum (applying chain rule: d/dx(f+g) = df/dx + dg/dx)
inline std::shared_ptr<ExpressionSum> dx(const std::shared_ptr<ExpressionSum> &u) {
    auto fun1_dx = dx(u->fun1);
    auto fun2_dx = dx(u->fun2);
    return std::make_shared<ExpressionSum>(fun1_dx, fun2_dx);
}

inline std::shared_ptr<ExpressionSum> dy(const std::shared_ptr<ExpressionSum> &u) {
    auto fun1_dy = dy(u->fun1);
    auto fun2_dy = dy(u->fun2);
    return std::make_shared<ExpressionSum>(fun1_dy, fun2_dy);
}

inline std::shared_ptr<ExpressionSum> dz(const std::shared_ptr<ExpressionSum> &u) {
    auto fun1_dz = dz(u->fun1);
    auto fun2_dz = dz(u->fun2);
    return std::make_shared<ExpressionSum>(fun1_dz, fun2_dz);
}

inline std::shared_ptr<ExpressionSum> dt(const std::shared_ptr<ExpressionSum> &u) {
    auto fun1_dt = dt(u->fun1);
    auto fun2_dt = dt(u->fun2);
    return std::make_shared<ExpressionSum>(fun1_dt, fun2_dt);
}

class ExpressionNormal2 : public ExpressionVirtual {
    typedef Mesh2 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uxnx, uyny;
    double c0 = 1;

  public:
    ExpressionNormal2(const FunFEM<M> &fh1, const Normal n)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(0);
        uyny.addNormal(1);
    }
    ExpressionNormal2(const FunFEM<M> &fh1, const Tangent t)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(1);
        uyny.addNormal(0);
        c0 = -1;
    }
    ExpressionNormal2(const FunFEM<M> &fh1, const Conormal n)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(0);
        uyny.addNormal(1);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return c0 * uxnx.evalOnBackMesh(k, dom, x, normal) + uyny.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return c0 * uxnx.evalOnBackMesh(k, dom, x, t, normal) + uyny.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionNormal2() {}
};
std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Normal &n);
std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Tangent &n);
std::shared_ptr<ExpressionNormal2> operator*(const FunFEM<Mesh2> &f1, const Conormal &n);

class ExpressionNormal2Q : public ExpressionVirtual {
    typedef MeshQuad2 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uxnx, uyny;
    double c0 = 1;

  public:
    ExpressionNormal2Q(const FunFEM<M> &fh1, const Normal n)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(0);
        uyny.addNormal(1);
    }
    ExpressionNormal2Q(const FunFEM<M> &fh1, const Tangent t)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(1);
        uyny.addNormal(0);
        c0 = -1;
    }
    ExpressionNormal2Q(const FunFEM<M> &fh1, const Conormal n)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(0);
        uyny.addNormal(1);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return c0 * uxnx.evalOnBackMesh(k, dom, x, normal) + uyny.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return c0 * uxnx.evalOnBackMesh(k, dom, x, t, normal) + uyny.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionNormal2Q() {}
};
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Normal &n);
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Tangent &n);
std::shared_ptr<ExpressionNormal2Q> operator*(const FunFEM<MeshQuad2> &f1, const Conormal &n);

class ExpressionNormal3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uxnx, uyny, uznz;

  public:
    ExpressionNormal3(const FunFEM<M> &fh1)
        : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0), uznz(fh1, 2, op_id, 0, 0) {
        assert(fh1.Vh->N != 1);
        uxnx.addNormal(0);
        uyny.addNormal(1);
        uznz.addNormal(2);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return uxnx.evalOnBackMesh(k, dom, x, normal) + uyny.evalOnBackMesh(k, dom, x, normal) +
               uznz.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return uxnx.evalOnBackMesh(k, dom, x, t, normal) + uyny.evalOnBackMesh(k, dom, x, t, normal) +
               uznz.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionNormal3() {}
};
std::shared_ptr<ExpressionNormal3> operator*(const FunFEM<Mesh3> &f1, const Normal &n);

/**
 * @brief n x u, where u is an expression of a 3D vector field
 *
 */
class ExpressionNormalCrossX3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uzny, uynz;

  public:
    ExpressionNormalCrossX3(const FunFEM<M> &fh1) : fun(fh1), uynz(fh1, 1, op_id, 0, 0), uzny(fh1, 2, op_id, 0, 0) {
        assert(fh1.Vh->N == 3); // assert that the function is a 3D vector field

        uzny.addNormal(1);
        uynz.addNormal(2);
    }

    virtual R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const { return uzny.eval(k, x, normal) - uynz.eval(k, x, normal); }

    R eval(const int k, const R *x, const R t, const R *normal) const {
        return uzny.eval(k, x, t, normal) - uynz.eval(k, x, t, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return uzny.evalOnBackMesh(k, dom, x, normal) - uynz.evalOnBackMesh(k, dom, x, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return uzny.evalOnBackMesh(k, dom, x, t, normal) - uynz.evalOnBackMesh(k, dom, x, t, normal);
    }

    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
};

class ExpressionNormalCrossY3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uxnz, uznx;

  public:
    ExpressionNormalCrossY3(const FunFEM<M> &fh1) : fun(fh1), uxnz(fh1, 0, op_id, 0, 0), uznx(fh1, 2, op_id, 0, 0) {
        assert(fh1.Vh->N == 3); // assert that the function is a 3D vector field

        uznx.addNormal(0);
        uxnz.addNormal(2);
    }

    virtual R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const { return uxnz.eval(k, x, normal) - uznx.eval(k, x, normal); }

    R eval(const int k, const R *x, const R t, const R *normal) const {
        return uxnz.eval(k, x, t, normal) - uznx.eval(k, x, t, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return uxnz.evalOnBackMesh(k, dom, x, normal) - uznx.evalOnBackMesh(k, dom, x, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return uxnz.evalOnBackMesh(k, dom, x, t, normal) - uznx.evalOnBackMesh(k, dom, x, t, normal);
    }

    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
};

class ExpressionNormalCrossZ3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> uynx, uxny;

  public:
    ExpressionNormalCrossZ3(const FunFEM<M> &fh1) : fun(fh1), uynx(fh1, 1, op_id, 0, 0), uxny(fh1, 0, op_id, 0, 0) {
        assert(fh1.Vh->N == 3); // assert that the function is a 3D vector field

        uynx.addNormal(0);
        uxny.addNormal(1);
    }

    virtual R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const { return uynx.eval(k, x, normal) - uxny.eval(k, x, normal); }

    R eval(const int k, const R *x, const R t, const R *normal) const {
        return uynx.eval(k, x, t, normal) - uxny.eval(k, x, t, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return uynx.evalOnBackMesh(k, dom, x, normal) - uxny.evalOnBackMesh(k, dom, x, normal);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return uynx.evalOnBackMesh(k, dom, x, t, normal) - uxny.evalOnBackMesh(k, dom, x, t, normal);
    }

    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
};

// class ExpressionNormalCross3 : public ExpressionVirtual {
//     typedef Mesh3 M;
//     const FunFEM<M> &fun;

//     ExpressionFunFEM<M> uxny, uxnz, uynx, uynz, uznx, uzny;

//   public:
//     // ExpressionNormal3(const FunFEM<M> &fh1)
//     //     : fun(fh1), uxnx(fh1, 0, op_id, 0, 0), uyny(fh1, 1, op_id, 0, 0), uznz(fh1, 2, op_id, 0, 0) {
//     //     assert(fh1.Vh->N != 1);
//     //     uxnx.addNormal(0);
//     //     uyny.addNormal(1);
//     //     uznz.addNormal(2);
//     // }

//     ExpressionNormalCross3(const FunFEM<M> &fh1)
//         : fun(fh1), uxny(fh1, 0, op_id, 0, 0), uxnz(fh1, 0, op_id, 0, 0), uynx(fh1, 1, op_id, 0, 0),
//           uynz(fh1, 1, op_id, 0, 0), uznx(fh1, 2, op_id, 0, 0), uzny(fh1, 2, op_id, 0, 0) {
//         assert(fh1.Vh->N == 3); // assert that the function is a 3D vector field

//         uxny.addNormal(1);
//         uxnz.addNormal(2);
//         uynx.addNormal(0);
//         uynz.addNormal(2);
//         uznx.addNormal(0);
//         uzny.addNormal(1);

//         fun_vector.at(0) = std::make_shared<ExpressionFunFEM<M>>(std::make_shared<ExpressionFunFEM<M>>(uynz) -
//                                                                  std::make_shared<ExpressionFunFEM<M>>(uzny));
//         fun_vector.at(1) = std::make_shared<ExpressionFunFEM<M>>(std::make_shared<ExpressionFunFEM<M>>(uznx) -
//                                                                  std::make_shared<ExpressionFunFEM<M>>(uxnz));
//         fun_vector.at(2) = std::make_shared<ExpressionFunFEM<M>>(std::make_shared<ExpressionFunFEM<M>>(uxny) -
//                                                                  std::make_shared<ExpressionFunFEM<M>>(uynx));
//     }

//     R operator()(long i) const {
//         assert(0);
//         return 0;
//     };

//     // // Do we evaluate vectors component-wise?
//     // R eval(const int k, const R *x, const R *normal) const {
//     //     std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
//     //     assert(0);
//     //     return 0;
//     // }
//     // R eval(const int k, const R *x, const R t, const R *normal) const {
//     //     std::cout << " evaluating f*n expression withoutr giving the normal as input " << std::endl;
//     //     assert(0);
//     //     return 0;
//     // }

//     // R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
//     //     assert(normal);
//     //     return uxnx.evalOnBackMesh(k, dom, x, normal) + uyny.evalOnBackMesh(k, dom, x, normal) +
//     //            uznz.evalOnBackMesh(k, dom, x, normal);
//     // }
//     // R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
//     //     assert(normal);
//     //     return uxnx.evalOnBackMesh(k, dom, x, t, normal) + uyny.evalOnBackMesh(k, dom, x, t, normal) +
//     //            uznz.evalOnBackMesh(k, dom, x, t, normal);
//     // }
//     int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
//     ~ExpressionNormalCross3() {}
// };

// std::vector<std::shared_ptr<ExpressionVirtual>> cross(const FunFEM<Mesh3> &f1, const Normal &n);
std::vector<std::shared_ptr<ExpressionVirtual>> cross(const Normal &n, const FunFEM<Mesh3> &f1);

/* 3D Curl of a FunFEM */
class ExpressionCurl3D {
  public:
    typedef Mesh3 M;
    const FunFEM<M> &fun;

    ExpressionCurl3D(const FunFEM<M> &fh1) : fun(fh1) {}

    std::vector<std::shared_ptr<ExpressionVirtual>> operator()() const {
        return {
            dy(fun.expr(2)) - dz(fun.expr(1)), // d/dy(u_z) - d/dz(u_y)
            dz(fun.expr(0)) - dx(fun.expr(2)), // d/dz(u_x) - d/dx(u_z)
            dx(fun.expr(1)) - dy(fun.expr(0))  // d/dx(u_y) - d/dy(u_x)
        };
    }
};

// Function to create and return the curl components directly
std::vector<std::shared_ptr<ExpressionVirtual>> curl(const FunFEM<Mesh3> &uh);

// divS for 2d
template <typeMesh M> class ExpressionDSx2 : public ExpressionVirtual {
    // typedef Mesh2 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> dxu1, dxu1nxnx, dyu1nxny;

  public:
    ExpressionDSx2(const FunFEM<M> &fh1, int ci = 0)
        : fun(fh1), dxu1(fh1, ci, op_dx, 0, 0), dxu1nxnx(fh1, ci, op_dx, 0, 0),
          dyu1nxny(fh1, ci, op_dy, 0, 0) {
        dxu1nxnx.addNormal(0);
        dxu1nxnx.addNormal(0);
        dyu1nxny.addNormal(0);
        dyu1nxny.addNormal(1);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dxu1.evalOnBackMesh(k, dom, x, normal) - dxu1nxnx.evalOnBackMesh(k, dom, x, normal) -
               dyu1nxny.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dxu1.evalOnBackMesh(k, dom, x, t, normal) - dxu1nxnx.evalOnBackMesh(k, dom, x, t, normal) -
               dyu1nxny.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDSx2() {}
};

template <typeMesh M> std::shared_ptr<ExpressionDSx2<M>> dxS(const FunFEM<M> &f1);
std::shared_ptr<ExpressionDSx2<Mesh2>> dxS(const FunFEM<Mesh2> &f1, int ci);

template <typeMesh M> class ExpressionDSy2 : public ExpressionVirtual {
    // typedef Mesh2 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> dxu2, dxu2nxny, dyu2nyny;

  public:
    ExpressionDSy2(const FunFEM<M> &fh1, int ci = 1)
        : fun(fh1), dxu2(fh1, ci, op_dy, 0, 0), dxu2nxny(fh1, ci, op_dx, 0, 0),
          dyu2nyny(fh1, ci, op_dy, 0, 0) {
        dxu2nxny.addNormal(0);
        dxu2nxny.addNormal(1);
        dyu2nyny.addNormal(1);
        dyu2nyny.addNormal(1);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dxu2.evalOnBackMesh(k, dom, x, normal) - dxu2nxny.evalOnBackMesh(k, dom, x, normal) -
               dyu2nyny.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dxu2.evalOnBackMesh(k, dom, x, t, normal) - dxu2nxny.evalOnBackMesh(k, dom, x, t, normal) -
               dyu2nyny.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDSy2() {}
};

template <typeMesh M> std::shared_ptr<ExpressionDSy2<M>> dyS(const FunFEM<M> &f1);
std::shared_ptr<ExpressionDSy2<Mesh2>> dyS(const FunFEM<Mesh2> &f1, int ci);

template <typeMesh M> class ExpressionDivS2 : public ExpressionVirtual {
    // typedef Mesh2 M;
    const FunFEM<M> &fun;
    const std::shared_ptr<ExpressionDSx2<M>> dx;
    const std::shared_ptr<ExpressionDSy2<M>> dy;
    // const ExpressionDSz<M> dz;

  public:
    ExpressionDivS2(const FunFEM<M> &fh1) : fun(fh1), dx(dxS(fh1)), dy(dyS(fh1)) {}

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dx->evalOnBackMesh(k, dom, x, normal) + dy->evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dx->evalOnBackMesh(k, dom, x, t, normal) + dy->evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDivS2() {}
};

template <typeMesh M> std::shared_ptr<ExpressionDivS2<M>> divS(const FunFEM<M> &f1);

// divs for 3D
class ExpressionDSx3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> dxu1, dxu1nxnx, dyu1nxny, dzu1nxnz;

  public:
    ExpressionDSx3(const FunFEM<M> &fh1,int ci)
        : fun(fh1), dxu1(fh1, ci, op_dx, 0, 0), dxu1nxnx(fh1, ci, op_dx, 0, 0), dyu1nxny(fh1, ci, op_dy, 0, 0),
          dzu1nxnz(fh1, ci, op_dz, 0, 0) {
        dxu1nxnx.addNormal(0);
        dxu1nxnx.addNormal(0);
        dyu1nxny.addNormal(0);
        dyu1nxny.addNormal(1);
        dzu1nxnz.addNormal(0);
        dzu1nxnz.addNormal(2);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dxu1.evalOnBackMesh(k, dom, x, normal) - dxu1nxnx.evalOnBackMesh(k, dom, x, normal) -
               dyu1nxny.evalOnBackMesh(k, dom, x, normal) - dzu1nxnz.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dxu1.evalOnBackMesh(k, dom, x, t, normal) - dxu1nxnx.evalOnBackMesh(k, dom, x, t, normal) -
               dyu1nxny.evalOnBackMesh(k, dom, x, t, normal) - dzu1nxnz.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDSx3() {}
};
std::shared_ptr<ExpressionDSx3> dxS(const FunFEM<Mesh3> &f1, int ci);

class ExpressionDSy3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> dxu2, dxu2nxny, dyu2nyny, dzu2nynz;

  public:
    ExpressionDSy3(const FunFEM<M> &fh1, int ci)
        : fun(fh1), dxu2(fh1, ci, op_dy, 0, 0), dxu2nxny(fh1, ci, op_dx, 0, 0), dyu2nyny(fh1, ci, op_dy, 0, 0),
          dzu2nynz(fh1, ci, op_dz, 0, 0) {


        dxu2nxny.addNormal(0);
        dxu2nxny.addNormal(1);
        dyu2nyny.addNormal(1);
        dyu2nyny.addNormal(1);
        dzu2nynz.addNormal(1);
        dzu2nynz.addNormal(2);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dxu2.evalOnBackMesh(k, dom, x, normal) - dxu2nxny.evalOnBackMesh(k, dom, x, normal) -
               dyu2nyny.evalOnBackMesh(k, dom, x, normal) - dzu2nynz.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dxu2.evalOnBackMesh(k, dom, x, t, normal) - dxu2nxny.evalOnBackMesh(k, dom, x, t, normal) -
               dyu2nyny.evalOnBackMesh(k, dom, x, t, normal) - dzu2nynz.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDSy3() {}
};
std::shared_ptr<ExpressionDSy3> dyS(const FunFEM<Mesh3> &f1, int ci);

class ExpressionDSz3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    ExpressionFunFEM<M> dxu3, dxu3nxnz, dyu3nynz, dzu3nznz;

  public:
    ExpressionDSz3(const FunFEM<M> &fh1, int ci)
        : fun(fh1), dxu3(fh1, ci, op_dz, 0, 0), dxu3nxnz(fh1, ci, op_dx, 0, 0), dyu3nynz(fh1, ci, op_dy, 0, 0),
          dzu3nznz(fh1, ci, op_dz, 0, 0) {
        dxu3nxnz.addNormal(0);
        dxu3nxnz.addNormal(2);
        dyu3nynz.addNormal(1);
        dyu3nynz.addNormal(2);
        dzu3nznz.addNormal(2);
        dzu3nznz.addNormal(2);
    }

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dxu3.evalOnBackMesh(k, dom, x, normal) - dxu3nxnz.evalOnBackMesh(k, dom, x, normal) -
               dyu3nynz.evalOnBackMesh(k, dom, x, normal) - dzu3nznz.evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dxu3.evalOnBackMesh(k, dom, x, t, normal) - dxu3nxnz.evalOnBackMesh(k, dom, x, t, normal) -
               dyu3nynz.evalOnBackMesh(k, dom, x, t, normal) - dzu3nznz.evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionDSz3() {}
};
std::shared_ptr<ExpressionDSz3> dzS(const FunFEM<Mesh3> &f1, int ci);

class ExpressionDivS3 : public ExpressionVirtual {
    typedef Mesh3 M;
    const FunFEM<M> &fun;
    const std::shared_ptr<ExpressionDSx3> dx;
    const std::shared_ptr<ExpressionDSy3> dy;
    const std::shared_ptr<ExpressionDSz3> dz;

  public:
    ExpressionDivS3(const FunFEM<M> &fh1) : fun(fh1), dx(dxS(fh1,0)), dy(dyS(fh1,1)), dz(dzS(fh1,2)) {}

    R operator()(long i) const {
        assert(0);
        return 0;
    };

    R eval(const int k, const R *x, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        std::cout << " evaluating DSx expression withoutr giving the normal as input " << std::endl;
        assert(0);
        return 0;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        return dx->evalOnBackMesh(k, dom, x, normal) + dy->evalOnBackMesh(k, dom, x, normal) +
               dz->evalOnBackMesh(k, dom, x, normal);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        return dx->evalOnBackMesh(k, dom, x, t, normal) + dy->evalOnBackMesh(k, dom, x, t, normal) +
               dz->evalOnBackMesh(k, dom, x, t, normal);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }

    ~ExpressionDivS3() {}
};
std::shared_ptr<ExpressionDivS3> divS(const FunFEM<Mesh3> &f1);

// template<typename M>
// class ExpressionAverage { //}: public ExpressionVirtual{
//   public:
//     // const ExpressionFunFEM<M> fun1;
//     // const ExpressionVirtual &fun1;
//     std::shared_ptr<ExpressionVirtual> fun1;
//     const R k1, k2;

//     // ExpressionAverage(const ExpressionFunFEM<M>& fh, double kk1, double kk2)
//     // : fun1(fh.fun,fh.cu, fh.op, fh.opt, -1), k1(kk1), k2(kk2){
//     // }
//     ExpressionAverage(const std::shared_ptr<ExpressionVirtual> &fh1, double kk1, double kk2)
//         : fun1(fh1), k1(kk1), k2(kk2) {}

//     R operator()(long i) const {
//         assert(0 && " cannot use this ");
//         return k1 * (*fun1)(i);
//     }

//     R eval(const int k, const R *x, const R *normal) const {
//         assert(0 && " need to be evaluated on backmesh");
//         return fun1->eval(k, x, normal) + fun1->eval(k, x, normal);
//     }
//     R eval(const int k, const R *x, const R t, const R *normal) const {
//         assert(0 && " need to be evaluated on backmesh");
//         return fun1->eval(k, x, t, normal) + fun1->eval(k, x, t, normal);
//     }

//     R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
//         return k1 * fun1->evalOnBackMesh(k, 0, x, normal) + k2 * fun1->evalOnBackMesh(k, 1, x, normal);
//     }
//     R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
//         return k1 * fun1->evalOnBackMesh(k, 0, x, t, normal) + k2 * fun1->evalOnBackMesh(k, 1, x, t, normal);
//     }
//     int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }

//     ~ExpressionAverage() {}
// };
class ExpressionAverage : public ExpressionVirtual {
    public:
      // const ExpressionFunFEM<M> fun1;
      // const ExpressionVirtual &fun1;
      std::shared_ptr<ExpressionVirtual>  fun1;
      const R k1, k2;
  
      // ExpressionAverage(const ExpressionFunFEM<M>& fh, double kk1, double kk2)
      // : fun1(fh.fun,fh.cu, fh.op, fh.opt, -1), k1(kk1), k2(kk2){
      // }
      ExpressionAverage(const std::shared_ptr<ExpressionVirtual> &fh1, double kk1, double kk2)
          : fun1(fh1), k1(kk1), k2(kk2) {}
  
      R operator()(long i) const {
          assert(0 && " cannot use this ");
          return k1 * (*fun1)(i);
      }
  
      R eval(const int k, const R *x, const R *normal) const {
          assert(0 && " need to be evaluated on backmesh");
          return fun1->eval(k, x, normal) + fun1->eval(k, x, normal);
      }
      R eval(const int k, const R *x, const R t, const R *normal) const {
          assert(0 && " need to be evaluated on backmesh");
          return fun1->eval(k, x, t, normal) + fun1->eval(k, x, t, normal);
      }
  
      R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(fun1.get());
            
          return k1 * fun1->evalOnBackMesh(k, 0, x, normal) + k2 * fun1->evalOnBackMesh(k, 1, x, normal);
      }
      R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
          return k1 * fun1->evalOnBackMesh(k, 0, x, t, normal) + k2 * fun1->evalOnBackMesh(k, 1, x, t, normal);
      }
      int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1->idxElementFromBackMesh(kb, dd); }
  
      ~ExpressionAverage() {}
  };
std::shared_ptr<ExpressionAverage> average(const std::shared_ptr<ExpressionVirtual> &fh1, const double kk1 = 0.5,
                                           const double kk2 = 0.5);
std::shared_ptr<ExpressionAverage> jump(const std::shared_ptr<ExpressionVirtual> &fh1, const double kk1 = 1,
                                        const double kk2 = -1);
std::shared_ptr<ExpressionAverage> operator*(double c, const ExpressionAverage &fh);
std::shared_ptr<ExpressionAverage> operator*(const ExpressionAverage &fh, double c);

// template<typename M>
// ExpressionAverage<M> average(const ExpressionVirtual & fh1, const double
// kk1=0.5, const double kk2=0.5) {
//   return ExpressionAverage<M>(fh1,kk1,kk2);
// }
// template<typename M>
// ExpressionAverage<M> jump(const ExpressionVirtual & fh1, const double kk1=1,
// const double kk2=-1){
//   return ExpressionAverage<M>(fh1,1,-1);
// }
// template<typename M>
// ExpressionAverage<M> operator* (double c, const ExpressionAverage<M>& fh){
//   return ExpressionAverage<M>(fh.fun1,c*fh.k1, c*fh.k2);
// }
// template<typename M>
// ExpressionAverage<M> operator* ( const ExpressionAverage<M>& fh, double c){
//   return ExpressionAverage<M>(fh.fun1,c*fh.k1, c*fh.k2);
// }

class ExpressionBurgerFlux : public ExpressionVirtual {
    const ExpressionVirtual &fun1;

  public:
    ExpressionBurgerFlux(const ExpressionVirtual &fh1) : fun1(fh1) {}

    R operator()(long i) const { return fabs(fun1(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        double val = fun1.eval(k, x, normal);
        return 0.5 * val * val;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        double val = fun1.eval(k, x, t, normal);
        return 0.5 * val * val;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double val = fun1.evalOnBackMesh(k, dom, x, normal);
        return 0.5 * val * val;
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double val = fun1.evalOnBackMesh(k, dom, x, t, normal);
        return 0.5 * val * val;
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1.idxElementFromBackMesh(kb, dd); }
    ~ExpressionBurgerFlux() {}
};
class ExpressionNormalBurgerFlux : public ExpressionVirtual {
    const ExpressionVirtual &fun1;

  public:
    ExpressionNormalBurgerFlux(const ExpressionVirtual &fh1) : fun1(fh1) {}

    R operator()(long i) const { return fabs(fun1(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        assert(normal);
        double val = fun1.eval(k, x, normal);
        return 0.5 * val * val * (normal[0] + normal[1]);
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        assert(normal);
        double val = fun1.eval(k, x, t, normal);
        return 0.5 * val * val * (normal[0] + normal[1]);
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        assert(normal);
        double val = fun1.evalOnBackMesh(k, dom, x, normal);
        return 0.5 * val * val * (normal[0] + normal[1]);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        assert(normal);
        double val = fun1.evalOnBackMesh(k, dom, x, t, normal);
        return 0.5 * val * val * (normal[0] + normal[1]);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun1.idxElementFromBackMesh(kb, dd); }
    ~ExpressionNormalBurgerFlux() {}
};
ExpressionBurgerFlux burgerFlux(const ExpressionVirtual &f1);
ExpressionNormalBurgerFlux burgerFlux(const ExpressionVirtual &f1, const Normal &n);

template <typename M> class ExpressionLinearSurfaceTension : public ExpressionVirtual {
    const FunFEM<M> &fun;
    const double sigma0;
    const double beta;
    const double tid;

  public:
    ExpressionLinearSurfaceTension(const FunFEM<M> &fh, double ssigma0, double bbeta, double ttid)
        : fun(fh), sigma0(ssigma0), beta(bbeta), tid(ttid) {}

    R operator()(long i) const { return fabs(fun(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        assert(0);
        return 0.;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        assert(0);
        return 0.;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double val = fun.evalOnBackMesh(k, dom, x, tid, 0, op_id, op_id);
        return sigma0 * (1 - beta * val);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double val = fun.evalOnBackMesh(k, dom, x, tid, 0, op_id, op_id);
        return sigma0 * (1 - beta * val);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionLinearSurfaceTension() {}
};
template <typename M> class ExpressionNonLinearSurfaceTension : public ExpressionVirtual {
    const FunFEM<M> &fun;
    const double sigma0;
    const double beta;
    const double tid;

  public:
    ExpressionNonLinearSurfaceTension(const FunFEM<M> &fh, double ssigma0, double bbeta, double ttid)
        : fun(fh), sigma0(ssigma0), beta(bbeta), tid(ttid) {}

    R operator()(long i) const { return fabs(fun(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        assert(0);
        return 0.;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        assert(0);
        return 0.;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double val = fun.evalOnBackMesh(k, dom, x, tid, 0, op_id, op_id);
        return sigma0 * (1 + beta * std::log(1 - val));
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double val = fun.evalOnBackMesh(k, dom, x, tid, 0, op_id, op_id);
        return sigma0 * (1 + beta * std::log(1 - val));
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionNonLinearSurfaceTension() {}
};

template <typename M> class ExpressionF : public ExpressionVirtual {
    const FunFEM<M> &fun;
    const FunFEM<M> &fun0;
    const double time;

  public:
    ExpressionF(const FunFEM<M> &ch, const FunFEM<M> &c0h, double time_)
        : fun(ch), fun0(c0h), time(time_) {}

    R operator()(long i) const { return fabs(fun(i)); }

    R eval(const int k, const R *x, const R *normal) const {
        assert(0);
        return 0.;
    }
    R eval(const int k, const R *x, const R t, const R *normal) const {
        assert(0);
        return 0.;
    }

    R evalOnBackMesh(const int k, const int dom, const R *x, const R *normal) const {
        double c  = fun.evalOnBackMesh(k, dom, x, time, 0, op_id, op_id);
        double c0 = fun0.evalOnBackMesh(k, dom, x, time, 0, op_id, op_id);
        assert(fabs(c0) > 1e-15);
        return 2.*c*c/(c0*c0 + c*c);
    }
    R evalOnBackMesh(const int k, const int dom, const R *x, const R t, const R *normal) const {
        double c = fun.evalOnBackMesh(k, dom, x, time, 0, op_id, op_id);
        double c0 = fun0.evalOnBackMesh(k, dom, x, time, 0, op_id, op_id);
        assert(fabs(c0) > 1e-15);

        return 2.*c*c/(c0*c0 + c*c);
    }
    int idxElementFromBackMesh(int kb, int dd = 0) const { return fun.idxElementFromBackMesh(kb, dd); }
    ~ExpressionF() {}
};

// Build a heap-stored level set to hold by unique_ptr inside the algoim cut
// classes (which need a stable, non-aliasing handle but cannot store a
// non-copyable FunFEM by value). For a FunFEM we build an independent non-owning
// *view* over the caller's dof data (shares the values, distinct object); for
// any other (analytic) level-set type we copy. The viewed/copied level set must
// outlive the holder. Shared by AlgoimInterface and AlgoimMacro.
template <typename M, typename L> std::unique_ptr<L> make_level_set_view(const L &src) {
    if constexpr (std::is_same_v<std::remove_cvref_t<L>, FunFEM<M>>) {
        std::span<double> data = src.array();
        return std::make_unique<L>(*src.Vh, data);
    } else {
        return std::make_unique<L>(src);
    }
}

#include "expression.tpp"

template <typeMesh M> std::shared_ptr<ExpressionDSx2<M>> dxS(const FunFEM<M> &f1) {
    return std::make_shared<ExpressionDSx2<M>>(f1);
}
template <typeMesh M> std::shared_ptr<ExpressionDSy2<M>> dyS(const FunFEM<M> &f1) {
    return std::make_shared<ExpressionDSy2<M>>(f1);
}
template <typeMesh M> std::shared_ptr<ExpressionDivS2<M>> divS(const FunFEM<M> &f1) {
    return std::make_shared<ExpressionDivS2<M>>(f1);
}

typedef FunFEM<Mesh2> Fun2_h;
typedef ExpressionFunFEM<Mesh2> Expression2;
typedef FunFEM<Mesh3> Fun3_h;
typedef ExpressionFunFEM<Mesh3> Expression3;

#endif
