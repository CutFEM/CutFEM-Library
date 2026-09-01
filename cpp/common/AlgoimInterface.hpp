
#ifndef ALGOIM_INTERFACE_HPP
#define ALGOIM_INTERFACE_HPP

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "interface_levelSet.hpp"
#include "../algoim/quadrature_general.hpp"
#include "../solver/solver.hpp"
#include "../algoim/algoim_quad_rule.hpp"

template <typename M, typename L> std::unique_ptr<L> make_level_set_view(const L &src);

/**
 * @brief Interface class that uses the algoim quadrature generation to find cut elements.
 *
 * @tparam M Mesh type
 * @tparam L Algoim level set type
 */
template <typeMesh M, typename L> class AlgoimInterface : public Interface<M> {
    using mesh_t         = M;
    using Element        = typename mesh_t::Element;
    using Rd             = typename mesh_t::Rd;
    static const int nve = Rd::d;
    using Face           = FaceInterface<nve>;
    using ElementIdx     = SortArray<Ubyte, Element::Rd::d + 1>;

    // save cut elements as a map of background mesh index of cut element, and corresponding interface quadrature rule

    // The level set is held through a unique_ptr (built by make_level_set_view)
    // so that a non-copyable FunFEM can be stored as an independent non-owning
    // view; the reference `phi` preserves the original by-value member interface
    // used throughout this class.  The viewed/copied level set must outlive this
    // interface (true for the per-time-node level sets owned by the geometry in
    // the active-surface driver).
    std::unique_ptr<L> phi_owned_;
    L &phi;
    const int quadrature_order = 5;
    int number_of_cut_elements{0};
    // std::map<int, algoim::QuadratureRule<2>> cut_elements;
    std::map<int, AlgoimQuadratureRule<M>> cut_elements;

  public:
    // Narrow-band cut detection.  make_algoim_patch() otherwise runs the
    // expensive algoim::quadGen on EVERY background element just to decide
    // whether it is cut; only the O(surface) cut cells actually need it.  When
    // narrow_band_diag_factor > 0 (and the level set is a FunFEM), an element
    // whose nodal level-set values are all strictly one sign and whose smallest
    // nodal magnitude exceeds narrow_band_diag_factor * (cell diagonal) is
    // skipped without calling quadGen.
    //
    // Rigorously safe (never skips a genuinely cut element) when the level set
    // is |grad phi| <= narrow_band_diag_factor Lipschitz: for any x in the cell
    // there is a node v within one diagonal, so |phi(x)| >= min_node|phi| -
    // factor*diag > 0 with constant sign.  The signed-distance level sets used
    // by the active-surface drivers satisfy this with factor = 1; the default
    // 1.5 adds margin for imperfect reinitialization / P2 curvature.  Set to 0
    // to recover the original quadGen-on-every-element behavior.
    static inline double narrow_band_diag_factor = 1.5;

    AlgoimInterface(const mesh_t &Mesh, const L &phi_, int label = 0);

    // Options-aware constructor: use the SAME bernstein degree / 1D quadrature
    // degree for the cut-classification rules as the assembly will use.  The
    // default constructor keeps the historical (bernstein 2, quad degree 5)
    // classification, which is only consistent with assemblies run at
    // bernstein 2; with higher-degree level sets the cut topology and the
    // assembly rules otherwise disagree on borderline cells.
    AlgoimInterface(const mesh_t &Mesh, const L &phi_, int bernstein_deg, int quad_deg, int label = 0);

    // Full quadrature configuration, including the MeshQuad2 backend.  Use
    // this when the stored cut-classification rule must be identical to the
    // rule used later by assembly and diagnostics.
    AlgoimInterface(const mesh_t &Mesh, const L &phi_, const ProblemOption &option, int label = 0);

    // std::map<int, algoim::QuadratureRule<2>> get_cut_elements() { return cut_elements; }
    std::map<int, AlgoimQuadratureRule<M>> get_cut_elements() { return cut_elements; }
    const AlgoimQuadratureRule<M> *get_cut_quadrature(int k) const;
    int get_nb_cut_elements() { return cut_elements.size(); }
    SignElement<Element> get_SignElement(int k) const override;
    Partition<Element> get_partition(int k) const override;

    Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k,
                                                         int ifac) const override;
    void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                       std::list<int> &erased_element, int sign_part) const override;

    double measure(const Face &f) const override;

    bool isCutFace(int k, int ifac) const override;

    bool isCut(int k) const override;

    Rd normal(int k, std::span<double> x) const override;

    R measure(int i) const override;

    Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const override;

    size_t size() const override { return cut_elements.size(); }

    double get_t() { return phi.t; }

  private:
    // Rule-generation parameters for cut classification (see options-aware
    // constructors); defaults preserve the historical behavior.
    ProblemOption patch_options_;

    void make_algoim_patch(int label);
};

#include "AlgoimInterface.tpp"

#endif
