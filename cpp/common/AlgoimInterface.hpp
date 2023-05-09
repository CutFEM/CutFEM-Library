
#ifndef ALGOIM_INTERFACE_HPP
#define ALGOIM_INTERFACE_HPP

#include "interface_levelSet.hpp"
#include "../algoim/quadrature_general.hpp"

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

    L phi;
    const int quadrature_order = 5;
    int number_of_cut_elements{0};
    std::map<int, algoim::QuadratureRule<2>> cut_elements;

  public:
    AlgoimInterface(const mesh_t &Mesh, const L &phi_, int label = 0);

    std::map<int, algoim::QuadratureRule<2>> get_cut_elements() { return cut_elements; }
    int get_nb_cut_elements() { return cut_elements.size(); }
    SignElement<Element> get_SignElement(int k) const override;
    Partition<Element> get_partition(int k) const override;

    Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k,
                                                         int ifac) const override;
    void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                       std::list<int> &erased_element, int sign_part) const override;

    double measure(const Face &f) const override;

    bool isCutFace(int k, int ifac) const override {
        assert(0);
        return 0;
    };

    bool isCut(int k) const override;

    Rd normal(int k, std::span<double> x) const override;

    R measure(int i) const override;

    Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const override;

    size_t size() const override { return cut_elements.size(); }

  private:
    void make_algoim_patch(int label);
};

#include "AlgoimInterface.tpp"

#endif