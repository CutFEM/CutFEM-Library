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
#ifndef BASE_CUTPROBLEM_HPP
#define BASE_CUTPROBLEM_HPP

// #include "baseProblem.hpp"

/**
 * @brief The actual solution object in CutFEM problems
 *
 * @tparam Mesh
 */
template <typename Mesh> class BaseCutFEM : public BaseFEM<Mesh> {

    typedef ActiveMesh<Mesh> CutMesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename FElement::QFB QFB;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;

    using mesh_t       = Mesh;
    using itemVFlist_t = ListItemVF<mesh_t>;

    int number_of_stabilized_edges;

  public:
    BaseCutFEM(const ProblemOption &option, int np) : BaseFEM<Mesh>(option, np) {}
    BaseCutFEM(const QuadratureFormular1d &qt, const ProblemOption &option, int np) : BaseFEM<Mesh>(qt, option, np) {}
    BaseCutFEM(const FESpace &vh, const ProblemOption &option, int np) : BaseFEM<Mesh>(vh, option, np) {}

    // Integral on K
    void addBilinear(const itemVFlist_t &, const CutMesh &);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, int itq, const TimeSlab &In);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, int itq);
    void addLinear(const itemVFlist_t &VF, const CutMesh &);
    void addLinear(const itemVFlist_t &VF, const CutMesh &, int itq, const TimeSlab &In);
    void addLinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In);
    void addLinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, int itq);
    virtual void addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                        double cst_time);

    void addBilinear(const itemVFlist_t &, const CutMesh &, const CExtension &, const int);
    void addLinear(const itemVFlist_t &, const CutMesh &, const CExtension &, const int);
    void addElementContributionOtherSide(const itemVFlist_t &, const int, const TimeSlab *, int, double);

    // integral on innerFace
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b, const TimeSlab &In);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b, const TimeSlab &In, int itq);
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b);
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b, const TimeSlab &In);
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CFacet &b, const TimeSlab &In, int itq);
    virtual void addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1,
                                     const std::pair<int, int> &e2, const TimeSlab *In, int itq, double cst_time);

    // integral on boundary
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, const TimeSlab &In,
                     std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, const TimeSlab &In, int itq,
                     std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, const TimeSlab &In,
                   std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const CutMesh &, const CBorder &b, const TimeSlab &In, int itq,
                   std::list<int> label = {});
    virtual void addBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int ifac,
                                       const TimeSlab *In, int itq, double cst_time);

    void setDirichlet(const FunFEM<Mesh> &gh, const CutMesh &Th, std::list<int> label = {});

    // integral on interface
    using BaseFEM<Mesh>::addBilinear;
    using BaseFEM<Mesh>::addLinear;
    using BaseFEM<Mesh>::addLagrangeMultiplier;

    // integral on inner Ridge / intersection with interface
    void addBilinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const CRidge &innerRidge,
                     std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const CRidge &innerRidge,
                     const TimeSlab &In, std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const CRidge &innerRidge,
                     const TimeSlab &In, int itq, std::list<int> label = {});

    void addLinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const CRidge &innerRidge,
                   std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const CRidge &innerRidge,
                   const TimeSlab &In, std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const CRidge &innerRidge,
                   const TimeSlab &In, int itq, std::list<int> label = {});
    void addInterfaceRidgeContribution(const itemVFlist_t &VF, const Interface<Mesh> &interface, int ifac,
                                       const TimeSlab *In, int itq, double cst_time);

    // Face stabilization
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In, int itq);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, int itq, const TimeSlab &In);
    void addFaceStabilizationSpecial(const itemVFlist_t &VF, const CutMesh &, int itq, const TimeSlab &In);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const MacroElement<Mesh> &);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const TimeMacroElement<Mesh> &);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const TimeMacroElement2<Mesh> &);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const MacroElementPartition<Mesh> &);
    template <typename L>
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const AlgoimMacro<Mesh, L> &);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const TimeMacroElementSurface<Mesh> &);
    void addFaceStabilizationRHS(const itemVFlist_t &VF, const CutMesh &Th, const MacroElement<Mesh> &macro);

    template <typename L>
    void addPatchStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In,
                              const AlgoimMacro<Mesh, L> &);
    void addPatchStabilization(const itemVFlist_t &VF, const CutMesh &);
    void addPatchStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In);

    // Lagrange multiplier
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &, const int k);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const TimeSlab &In);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const TimeSlab &In, int itq,
                               bool init = true);
    void addLagrangeContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq, double cst_time);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &, const CBorder &b,
                               std::list<int> label = {});
    void addLagrangeBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int ifac,
                                       const TimeSlab *In, int itq, double cst_time);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const CExtension &ext,
                               const int epsE);
    void addLagrangeContributionOtherSide(const itemVFlist_t &VF, const int k, const int epsE);

    void addLagrangeVecToRowAndCol(const std::span<double> vecRow, const std::span<double> vecCol, const R val_rhs);

    void removeDofForHansbo(const FESpace &Vh);

    // For time problem
    template <typename V>
        requires NonAllocVector<V> || std::is_same_v<V, KN<typename V::element_type>>
    void saveSolution(const V);

    template <typename V>
        requires NonAllocVector<V> || std::is_same_v<V, KN<typename V::element_type>>
    void initialSolution(V &);

    int get_number_of_stabilized_edges() { return number_of_stabilized_edges; }
};

/**
 * @brief Create main object for a CutFEM problem
 *
 * @tparam Mesh
 */
template <typename Mesh> class CutFEM : public BaseCutFEM<Mesh>, public Solver {
    typedef GFESpace<Mesh> FESpace;
    typedef std::map<std::pair<int, int>, R> Matrix;

  public:
    CutFEM(const ProblemOption &option = defaultProblemOption) : BaseCutFEM<Mesh>(option, 1), Solver(option) {}

    CutFEM(const QuadratureFormular1d &qt, const ProblemOption &option = defaultProblemOption)
        : BaseCutFEM<Mesh>(qt, option, 1), Solver(option) {}

    CutFEM(const FESpace &vh, const ProblemOption &option = defaultProblemOption)
        : BaseCutFEM<Mesh>(vh, option, 1), Solver(option) {}

    CutFEM(const QuadratureFormular1d &qt, int np, const ProblemOption &option = defaultProblemOption)
        : BaseCutFEM<Mesh>(qt, option, np), Solver(option) {}

    CutFEM(const FESpace &vh, int np, const ProblemOption &option = defaultProblemOption)
        : BaseCutFEM<Mesh>(vh, option, np), Solver(option) {}

    void solve() {
        gather(this->mat_);
        Solver::solve(this->mat_[0], this->rhs_);
    }
    void solve(std::string solverName) {
        this->solver_name_ = solverName;
        gather(this->mat_);
        Solver::solve(this->mat_[0], this->rhs_);
    }
    void solve(std::map<std::pair<int, int>, R> &A, Rn &b) { Solver::solve(A, b); }
    void solve(std::vector<Matrix> &A, Rn &b, std::string solverName) {
        gather(A);
        Solver::solve(A[0], b);
    }
};

#include "baseCutProblem_Function.hpp"
#endif
