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
#ifndef BASE_PROBLEM_HPP
#define BASE_PROBLEM_HPP

#include "problem.hpp"

template <typename Mesh> class BaseFEM : public ShapeOfProblem<Mesh>, public QuadratureOfProblem<Mesh> {

    typedef ActiveMesh<Mesh> CutMesh;
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename FElement::QFB QFB;
    typedef typename Mesh::BorderElement BorderElement;

    using mesh_t       = Mesh;
    using itemVF_t     = ItemVF<mesh_t>;
    using itemVFlist_t = ListItemVF<mesh_t>;

  protected:
    int N_component_max_ = 0;
    int df_loc_max_      = 0;
    double *databf_      = nullptr;
    double *databf_time_ = nullptr;
    long offset_bf_      = 0;
    long offset_bf_time  = 0;

  public:
    ~BaseFEM() {
        delete[] databf_;
        delete[] databf_time_;
    }
    BaseFEM(const BaseFEM&)            = delete;              // or implement deep copy
    BaseFEM& operator=(const BaseFEM&) = delete;              // or implement deep copy
    BaseFEM(BaseFEM&& other) noexcept
        : N_component_max_(other.N_component_max_), df_loc_max_(other.df_loc_max_),
        databf_(other.databf_), databf_time_(other.databf_time_),
        offset_bf_(other.offset_bf_), offset_bf_time(other.offset_bf_time) {
        other.databf_ = nullptr;
        other.databf_time_ = nullptr;
    }
    BaseFEM& operator=(BaseFEM&& other) noexcept {
        if (this != &other) {
            delete[] databf_;
            delete[] databf_time_;
            N_component_max_ = other.N_component_max_;
            df_loc_max_      = other.df_loc_max_;
            databf_          = other.databf_;
            databf_time_     = other.databf_time_;
            offset_bf_       = other.offset_bf_;
            offset_bf_time   = other.offset_bf_time;
            other.databf_ = nullptr;
            other.databf_time_ = nullptr;
        }
        return *this;
    }

    BaseFEM(const ProblemOption &option) : ShapeOfProblem<Mesh>(), QuadratureOfProblem<Mesh>(option) {}
    BaseFEM(const QuadratureFormular1d &qt, const ProblemOption &option)
        : ShapeOfProblem<Mesh>(), QuadratureOfProblem<Mesh>(qt, option) {}
    BaseFEM(const FESpace &vh, const ProblemOption &option) : BaseFEM<Mesh>(option) {
        this->mapIdx0_.clear();
        this->mapIdx0_[&vh] = 0;
        int ndf             = vh.NbDoF();
        this->init(ndf);

        N_component_max_  = vh.N;
        df_loc_max_       = vh[0].NbDoF();
        offset_bf_        = 5 * df_loc_max_ * N_component_max_ * op_DDall;
        long size_data_bf = this->thread_count_max_ * offset_bf_;
        databf_           = new double[size_data_bf];
    }

    void initSpace(const FESpace &Vh, const TimeSlab &In) {

        this->mapIdx0_.clear();
        this->mapIdx0_[&Vh] = 0;

        this->init(Vh.NbDoF() * In.NbDoF(), In.NbDoF());

        N_component_max_ = Vh.N;
        df_loc_max_      = Vh[0].NbDoF();

        // free old space buffers (array delete)    
        delete[] databf_;

        offset_bf_        = 5 * df_loc_max_ * N_component_max_ * op_DDall;
        long size_data_bf = this->thread_count_max_ * offset_bf_;
        databf_     = new double[size_data_bf];

        // free old time buffers (array delete)
        delete[] databf_time_;

        offset_bf_time     = In.NbDoF() * op_DDall;
        size_data_bf       = this->thread_count_max_ * offset_bf_time;
        databf_time_ = new double[size_data_bf];
    }

    void initSpace(const FESpace &Vh) {
        this->mapIdx0_.clear();
        this->mapIdx0_[&Vh] = 0;

        this->init(Vh.NbDoF());

        N_component_max_ = Vh.N;
        df_loc_max_      = Vh[0].NbDoF();

        delete[] databf_;

        offset_bf_        = 5 * df_loc_max_ * N_component_max_ * op_DDall;
        long size_data_bf = this->thread_count_max_ * offset_bf_;
        databf_     = new double[size_data_bf];
    }

    void add(const FESpace &Qh) {
        this->mapIdx0_[&Qh] = this->get_nb_dof();
        int ndf             = this->get_nb_dof() + Qh.NbDoF();
        this->init(ndf);

        N_component_max_ = std::max(N_component_max_, Qh.N);
        df_loc_max_      = std::max(df_loc_max_, Qh[0].NbDoF());

        // was: delete databf_;  -> must be array delete
        delete[] databf_;

        offset_bf_        = 5 * df_loc_max_ * N_component_max_ * op_DDall;
        long size_data_bf = this->thread_count_max_ * offset_bf_;
        databf_           = new double[size_data_bf];
    }
    void add(const FESpace &Vh, const TimeSlab &In) {
        assert(this->nb_dof_time_ == In.NbDoF());

        this->mapIdx0_[&Vh] = this->get_nb_dof();
        int ndf             = this->get_nb_dof() + (Vh.NbDoF() * In.NbDoF());
        this->init(ndf);

        N_component_max_ = std::max(N_component_max_, Vh.N);
        df_loc_max_      = std::max(df_loc_max_, Vh[0].NbDoF());

        // was: delete databf_;  -> must be array delete
        delete[] databf_;

        offset_bf_        = 5 * df_loc_max_ * N_component_max_ * op_DDall;
        long size_data_bf = this->thread_count_max_ * offset_bf_;
        databf_           = new double[size_data_bf];
    }

    // add local to Matrix or rhs
    void addToMatrix(const itemVF_t &VF, const FElement &FKu, const FElement &FKv, const RNMK_ &fu, const RNMK_ &fv,
                     double Cint);
    // void addToMatrix(const itemVF_t& VFi, const FElement& FKu, const
    // FElement& FKv, const RNMK_& fu, const RNMK_& fv, double Cint, int
    // id_thread) ;
    void addToMatrix_Opt(const itemVF_t &VF, const FElement &FK, const RNMK_ &fv, double Cint);

    void addToMatrix(const itemVF_t &VF, const TimeSlab &In, const FElement &FKu, const FElement &FKv, const RNMK_ &fu,
                     const RNMK_ &fv, double Cint);

    void addToMatrixSpecial(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKu, const FElement &FKv,
                            const RNMK_ &fu, const RNMK_ &fv, double Cint, int from_time_dof_u, int from_time_dof_v);

    void addToMatrix(const itemVF_t &VFi, const FElement &FKv, const RNMK_ &fv, double Cint);
    void addToMatrix(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKv, const RNMK_ &fv, double Cint);

    void addToRHS(const itemVF_t &VF, const FElement &FKv, const RNMK_ &fv, double Cint);
    void addToRHS(const itemVF_t &VF, const TimeSlab &In, const FElement &FKv, const RNMK_ &fv, double Cint);

    // add diagonal terms
    void addDiagonal(const FESpace &Qh, double val);
    void setDiagonal(const FESpace &Qh, double val);

    // Integral on K
    void addBilinear(const itemVFlist_t &VF, const Mesh &);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &);
    void addBilinear(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &);
    void addLinear(const itemVFlist_t &VF, const CutMesh &);
    void addLinear(const itemVFlist_t &VF, const Mesh &);
    template <typename Fct>
    void addLinear(const Fct &f, const itemVFlist_t &VF, const Mesh &);
    void addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq, double cst_time);
    template <typename Fct>
    void addElementContribution(const Fct &f, const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                double cst_time);
    void addElementInterfaceContribution(const itemVFlist_t &VF, const Interface<Mesh> &interface, int ifac, int k,
                                          const TimeSlab *In, double cst_time, int itq);

    // integral on innerFace
    void addBilinear(const itemVFlist_t &VF, const Mesh &, const CFacet &b);
    void addLinear(const itemVFlist_t &VF, const Mesh &, const CFacet &b);
    void addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1, const std::pair<int, int> &e2,
                             const TimeSlab *In, int itq, double cst_time);
    void addFaceContributionMixed(const itemVFlist_t &VF, const std::pair<int, int> &e1, const std::pair<int, int> &e2,
                                  const TimeSlab *In, int itq, double cst_time);
    void addFaceContributionSpecial(const itemVFlist_t &VF, const std::pair<int, int> &e1,
                                    const std::pair<int, int> &e2, const TimeSlab *In, int itq, double cst_time);

    // integral on patches
    void addPatchContribution(const itemVFlist_t &VF, const int k, const int kn, const TimeSlab *In, int itq,
                              double cst_time);
    void addPatchContributionMixed(const itemVFlist_t &VF, const int kb, const int kbn, const TimeSlab *In, int itq,
                                   double cst_time);

    // integral on boundary
    void addBilinear(const itemVFlist_t &VF, const Mesh &, const CBorder &b, std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const Mesh &, const CBorder &b, std::list<int> label = {});
    void addBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int ifac,
                               const TimeSlab *In, int itq, double cst_time);

    // integral on active mesh boundary
    void addInnerBorderContribution(const itemVFlist_t &VF, const int k, const int ifac,
                                    const TimeSlab *In, int itq, double cst_time);
    void addInnerBorderContributionMixed(const itemVFlist_t &VF, const int kb, const int ifac,
                                    const TimeSlab *In, int itq, double cst_time);                                    
    void addOuterBorderContribution(const itemVFlist_t &VF, const int k, const int ifac,
                                    const TimeSlab *In, int itq, double cst_time);
    void addOuterBorderContributionMixed(const itemVFlist_t &VF, const int kb, const int ifac,
                                    const TimeSlab *In, int itq, double cst_time);
    
                                    

    void setDirichlet(const FunFEM<Mesh> &, const Mesh &, std::list<int> label = {});
    void setDirichletHone(const FunFEM<Mesh> &gh, const Mesh &Th, std::list<int> label = {});

    // integral on interface
    void addBilinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const TimeSlab &In, int itq,
                     std::list<int> label = {});

    void addBilinearAbsTest(const itemVFlist_t &VF, const Interface<Mesh> &gamma, std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In,
                     std::list<int> label = {});
    void addBilinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In, int itq,
                     std::list<int> label = {});
    template <typename Fct>
    void addLinear(const Fct &f, const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In,
                   std::list<int> label = {});
    template <typename Fct>
    void addLinear(const Fct &f, const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In, int itq,
                   std::list<int> label = {});
    template <typename Fct>
    void addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<Mesh> &gamma, int ifac,
                                  double tid, const TimeSlab *In, double cst_time, int itq);

    void addLinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, std::list<int> label = {});
    template <typename Fct>
    void addLinear(const Fct &f, const itemVFlist_t &VF, const Interface<Mesh> &gamma, std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const TimeSlab &In, int itq,
                   std::list<int> label = {});
    template <typename Fct>
    void addLinear(const Fct &f, const itemVFlist_t &VF, const Interface<Mesh> &gamma, const TimeSlab &In, int itq,
                   std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In,
                   std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const TimeInterface<Mesh> &gamma, const TimeSlab &In, int itq,
                   std::list<int> label = {});
    virtual void addInterfaceContribution(const itemVFlist_t &VF, const Interface<Mesh> &gamma, int ifac, double tid,
                                          const TimeSlab *In, double cst_time, int itq);

    void addBilinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const Mapping<Mesh> &,
                     std::list<int> label = {});
    void addLinear(const itemVFlist_t &VF, const Interface<Mesh> &gamma, const Mapping<Mesh> &,
                   std::list<int> label = {});
    void addInterfaceContribution(const itemVFlist_t &VF, const Interface<Mesh> &gamma, int ifac, double tid,
                                  const TimeSlab *In, double cst_time, int itq, const Mapping<Mesh> &);

    // integral for Lagrange multiplier
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const Mesh &Th);
    void addLagrangeContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq, double cst_time);
    void addLagrangeMultiplier(const itemVFlist_t &VF, double val, const Interface<Mesh> &gamma);
    virtual void addLagrangeContribution(const itemVFlist_t &VF, const Interface<mesh_t> &interface, const int iface);

    void addLagrangeBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int ifac,
                                       const TimeSlab *In, int itq, double cst_time);

    void addLagrangeVecToRowAndCol(const std::span<double> vecRow, const std::span<double> vecCol, const R val_rhs);

    // integral for global face stabilization
    void addPatchStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In);
    void addFaceStabilization(const itemVFlist_t &VF, const CutMesh &, const TimeSlab &In);
};

template <typename Mesh> class FEM : public BaseFEM<Mesh>, public Solver {
    typedef GFESpace<Mesh> FESpace;
    // typedef GenericMapping<Mesh> Mapping;

  public:
    FEM(const QuadratureFormular1d &qt, const ProblemOption &option = defaultProblemOption)
        : BaseFEM<Mesh>(qt, option), Solver(option) {}
    FEM(const FESpace &vh, const ProblemOption &option = defaultProblemOption)
        : BaseFEM<Mesh>(vh, option), Solver(option) {}

    void solve() { Solver::solve(this->mat_, this->rhs_); }
    void solve(std::string solverName) {
        this->solver_name_ = solverName;
        Solver::solve(this->mat_, this->rhs_);
    }
    void solve(std::map<std::pair<int, int>, R> &A, Rn &b) { Solver::solve(A, b); }
};

#include "baseProblem_Function.hpp"
#include "baseCutProblem.hpp"
#include "AlgoimBaseCutProblem.hpp"
#include "AlgoimCutFEM.hpp"

#endif
