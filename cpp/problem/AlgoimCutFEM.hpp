#ifndef ALGOIM_CUTFEM_HPP
#define ALGOIM_CUTFEM_HPP


// #include "../cutfem.hpp"
#include "../algoim/quadrature_general.hpp"
#include "../algoim/quadrature_multipoly.hpp"
#include <iostream>

template <typeMesh Mesh>
struct AlgoimQuadratureRule {

    std::vector<typename Mesh::Rd> points;    // Physical coordinates
    std::vector<double> weights;              // Quadrature weights
    std::vector<typename Mesh::Rd> normals;   // Normals (for surfaces)
    
    bool empty() const {return points.empty();}
    bool size() const {return points.size();}

    void reserve(size_t n) {
        points.reserve(n);
        weights.reserve(n);
        normals.reserve(n);
    }
};

// Specializations for different mesh types - implemented in AlgoimCutFEM.tpp
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenVol(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, double time = 0.);
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenSurf(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, double time = 0.);
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenFace(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, int ifac, double time = 0.);


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

    AlgoimCutFEMUnified(const fespace_t& vh, Phi& phi, const ProblemOption& opt)
      : BaseCutFEM<mesh_t>(vh, opt),
        Solver(opt),
        phi_(phi),
        options_(opt)
        {}

    // Integral overrides
    void addElementContribution(const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq,
                                double cst_time) override;

    void addInterfaceContribution(const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
                                  const TimeSlab* In, double cst_time, int itq) override;

    void addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1, 
                            const std::pair<int, int> &e2, const TimeSlab *In, 
                            int itq, double cst_time) override;


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