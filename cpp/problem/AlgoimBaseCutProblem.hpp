

#ifndef BASE_CUTPROBLEM_SAYE_HPP
#define BASE_CUTPROBLEM_SAYE_HPP

#include "../algoim/quadrature_general.hpp"

template <typename M, typename L> class AlgoimBaseCutFEM : public BaseCutFEM<M> {

    using mesh_t        = M;
    using fespace_t     = GFESpace<mesh_t>;
    using itemVFlist_t  = ListItemVF<mesh_t>;
    using FElement      = typename fespace_t::FElement;
    using Rd            = typename FElement::Rd;
    using QF            = typename FElement::QF;
    using QFB           = typename FElement::QFB;
    using Element       = typename mesh_t::Element;
    using BorderElement = typename mesh_t::BorderElement;

    int quadrature_order = 5;

    L phi;

  public:

    void addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                double cst_time) override;

    void addInterfaceContribution(const itemVFlist_t &VF, const Interface<mesh_t> &interface, int ifac, double tid,
                                  const TimeSlab *In, double cst_time, int itq) override;

    void addLagrangeContribution(const itemVFlist_t &VF, const Interface<mesh_t> &interface, const int iface) override;


    // template <typename Fct>
    // void addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<M> &interface, int ifac,
    //                               double tid, const TimeSlab *In, double cst_time, int itq) override;

    AlgoimBaseCutFEM(const QuadratureFormular1d &qt, L &phi_, const ProblemOption &option, int np)
        : BaseCutFEM<mesh_t>(qt, option, np), phi(phi_) {}

    AlgoimBaseCutFEM(const fespace_t &vh, L &phi_, const ProblemOption &option, int np)
        : BaseCutFEM<mesh_t>(vh, option, np), phi(phi_) {}
};

template <meshQuadrilateral M, typename L> class AlgoimCutFEM : public AlgoimBaseCutFEM<M, L>, public Solver {
    typedef GFESpace<M> fespace_t;
    typedef std::map<std::pair<int, int>, R> Matrix;

  public:
    // AlgoimCutFEM(const ProblemOption &option = defaultProblemOption) : BaseCutFEM<mesh_t>(option, 1), Solver(option)
    // {}

    AlgoimCutFEM(const QuadratureFormular1d &qt, L &phi, const ProblemOption &option = defaultProblemOption)
        : AlgoimBaseCutFEM<M, L>(qt, phi, option, 1), Solver(option) {}

    AlgoimCutFEM(const fespace_t &vh, L &phi, const ProblemOption &option = defaultProblemOption)
        : AlgoimBaseCutFEM<M, L>(vh, phi, option, 1), Solver(option) {}

    // AlgoimCutFEM(const QuadratureFormular1d &qt, int np, const ProblemOption &option = defaultProblemOption)
    //     : BaseCutFEM<mesh_t>(qt, option, np), Solver(option) {}

    // AlgoimCutFEM(const fespace_t &vh, int np, const ProblemOption &option = defaultProblemOption)
    //     : BaseCutFEM<mesh_t>(vh, option, np), Solver(option) {}

    void solve() {
        gather(this->mat_);
        Solver::solve(this->mat_[0], this->rhs_);
    }
    void solve(std::string solverName) {
        this->solver_name_ = solverName;
        gather(this->mat_);
        Solver::solve(this->mat_[0], this->rhs_);
    }
    // void solve(std::map<std::pair<int, int>, R> &A, Rn &b) { Solver::solve(A, b); }
    // void solve(std::vector<Matrix> &A, Rn &b, std::string solverName) {
    //     gather(A);
    //     Solver::solve(A[0], b);
    // }
};

#include "AlgoimBaseCutProblem.tpp"

#endif