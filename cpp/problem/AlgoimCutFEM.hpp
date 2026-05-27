#ifndef ALGOIM_CUTFEM_HPP
#define ALGOIM_CUTFEM_HPP


// #include "../cutfem.hpp"
#include "../algoim/quadrature_general.hpp"
#include "../algoim/quadrature_multipoly.hpp"
#include "../algoim/algoim_quad_rule.hpp"
#include <iostream>


// Specializations for different mesh types - implemented in AlgoimCutFEM.tpp
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenVol(const Mesh2::Element& K, Phi& phi, const ProblemOption& option);
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenSurf(const Mesh2::Element& K, Phi& phi, const ProblemOption& option);
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenFace(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, int ifac);

template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenVol(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option);
template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenSurf(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option);

// Gustaf
template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenVol(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option);
template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenSurf(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option);
template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenFace(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option, int ifac);

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenVol(const Mesh3::Element& K, Phi& phi, const ProblemOption& option);

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenSurf(const Mesh3::Element& K, Phi& phi, const ProblemOption& option);

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenFace(const Mesh3::Element& K, Phi& phi, const ProblemOption& option, int ifac);



template <typeMesh Mesh, typename Phi>
class AlgoimCutFEMUnified : public BaseCutFEM<Mesh>, public Solver {

    using mesh_t        = Mesh;
    using fespace_t     = GFESpace<mesh_t>;
    using itemVFlist_t  = ListItemVF<mesh_t>;
    using fe_element_t  = typename fespace_t::FElement;
    using element_t     = typename mesh_t::Element;
    using Rd            = typename fe_element_t::Rd;
    using matrix_t      = std::map<std::pair<int, int>, R>;

    ProblemOption options_;
    Phi &phi_;

public:

    AlgoimCutFEMUnified(const fespace_t& vh, Phi& phi, const ProblemOption& opt = defaultProblemOption)
      : BaseCutFEM<mesh_t>(vh, opt),
        Solver(opt),
        phi_(phi),
        options_(opt)
        {}

    AlgoimCutFEMUnified(const QuadratureFormular1d &qt, Phi &phi, const ProblemOption &option = defaultProblemOption)
        : BaseCutFEM<mesh_t>(qt, option), 
        Solver(option),
        phi_(phi),
        options_(option)
        {}
        
    // Integral overrides
    void addElementContribution(const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq,
                                double cst_time) override;

    // template <typename Fct>
    // void addElementContributionExact(const Fct &f, const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq,
    //                             double cst_time);

    void addInterfaceContribution(const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
                                  const TimeSlab* In, double cst_time, int itq) override;

    // template <typename Fct>
    // void addInterfaceContributionExact(const Fct &f, const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
    //                               const TimeSlab* In, double cst_time, int itq);

    // Gustaf
    void addBilinear(const itemVFlist_t& VF, const Interface<mesh_t>& interface);
    void addLinear(const itemVFlist_t& VF, const Interface<mesh_t>& interface);

    void addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1, 
                            const std::pair<int, int> &e2, const TimeSlab *In, 
                            int itq, double cst_time) override;

    // void addBilinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th);
    // void addLinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th);
    
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th, const TimeSlab &In);
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th, const int itq, const TimeSlab &In, const double scaling_time = 1.);
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th, const TimeSlab &In, const int itq);
    
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const TimeInterface<mesh_t> &gamma, const TimeSlab &In);
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const Interface<mesh_t> &gamma, const TimeSlab &In, const int itq);
    // template <typename Fct>
    // void addLinearExact(const Fct &f, const itemVFlist_t &VF, const TimeInterface<mesh_t> &gamma, const TimeSlab &In, const int itq);
    
    
    // Solver methods
    void solve() { Solver::solve(this->mat_, this->rhs_); }
    void solve(std::string solverName) {
        this->solver_name_ = solverName;
        Solver::solve(this->mat_, this->rhs_);
    }
    void solve(matrix_t &A, std::span<double> b) { Solver::solve(A, b); }
    void solve(matrix_t &A, std::span<double> b, std::string solverName) { 
        this->solver_name_ = solverName;    
        Solver::solve(A, b); 
    }

};





#include "AlgoimCutFEM.tpp"


#endif // ALGOIM_CUTFEM_HPP