

#ifndef BASE_CUTPROBLEM_SAYE_HPP
#define BASE_CUTPROBLEM_SAYE_HPP

#include "../algoim/quadrature_general.hpp"

template <typename Mesh, typename L> class AlgoimBaseCutFEM : public BaseCutFEM<Mesh> {

    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename FElement::QFB QFB;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;

    L phi;

  public:

    void addElementContribution(const ListItemVF<Rd::d> &VF, const int k, const TimeSlab *In, int itq, double cst_time);

	void addInterfaceContribution(const ListItemVF<Rd::d> &VF, const Interface<Mesh> &interface, int ifac,
									double tid, const TimeSlab *In, double cst_time, int itq);

    AlgoimBaseCutFEM(const QuadratureFormular1d &qt, L &phi_, const ProblemOption &option, int np) : BaseCutFEM<Mesh>(qt, option, np), phi(phi_) {}

    AlgoimBaseCutFEM(const FESpace &vh, L &phi_, const ProblemOption &option, int np) : BaseCutFEM<Mesh>(vh, option, np), phi(phi_) {}

};

template <meshQuadrilateral M, typename L>
class AlgoimCutFEM : public AlgoimBaseCutFEM<M, L>, public Solver {
    typedef GFESpace<M> FESpace;
    typedef std::map<std::pair<int, int>, R> Matrix;

  public:
    // AlgoimCutFEM(const ProblemOption &option = defaultProblemOption) : BaseCutFEM<Mesh>(option, 1), Solver(option) {}

    AlgoimCutFEM(const QuadratureFormular1d &qt, L &phi, const ProblemOption &option = defaultProblemOption)
        : AlgoimBaseCutFEM<M, L>(qt, phi, option, 1), Solver(option) {}

    AlgoimCutFEM(const FESpace &vh, L &phi, const ProblemOption &option = defaultProblemOption)
        : AlgoimBaseCutFEM<M, L>(vh, phi, option, 1), Solver(option) {}

    // AlgoimCutFEM(const QuadratureFormular1d &qt, int np, const ProblemOption &option = defaultProblemOption)
    //     : BaseCutFEM<Mesh>(qt, option, np), Solver(option) {}

    // AlgoimCutFEM(const FESpace &vh, int np, const ProblemOption &option = defaultProblemOption)
    //     : BaseCutFEM<Mesh>(vh, option, np), Solver(option) {}

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