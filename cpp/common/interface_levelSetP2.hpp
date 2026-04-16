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
#ifndef COMMON_LEVELSET_INTERFACE_P2_HPP
#define COMMON_LEVELSET_INTERFACE_P2_HPP

#include <limits>
#include "base_interface.hpp"
#include "dataStruct2D.hpp"
#include "dataStruct3D.hpp"
#include "../FESpace/operationTypeFE.hpp"

// ---------------------------------------------------------------------------
// Sub-element connectivity table for P2 elements.
// For each mesh element type, lists the n_sub sub-elements (each with nv
// vertex indices into the local 0..(nv+ne-1) P2 node array) that partition
// the reference element when edge midpoints are added.
//
// 2D Triangle (nv=3, ne=3, 6 nodes total → 4 sub-triangles):
//   nodes 0-2 = vertices, 3=mid(v1,v2), 4=mid(v2,v0), 5=mid(v0,v1)
//
// 3D Tet (nv=4, ne=6, 10 nodes total → 8 sub-tets):
//   nodes 0-3 = vertices,
//   4=mid(v0,v1), 5=mid(v0,v2), 6=mid(v0,v3),
//   7=mid(v1,v2), 8=mid(v1,v3), 9=mid(v2,v3)
//   Inner 4 tets use diagonal 6–7 (mid(v0,v3)–mid(v1,v2)).
// ---------------------------------------------------------------------------
template <typename Element> struct SubElementTable; // undefined base

template <> struct SubElementTable<Triangle2> {
    static constexpr int n_sub       = 4;
    static constexpr int idx[4][3]   = {{0,5,4},{1,3,5},{2,4,3},{5,3,4}};
};

template <> struct SubElementTable<Tet> {
    static constexpr int n_sub       = 8;
    static constexpr int idx[8][4]   = {
        {0,4,5,6},  // corner tet at v0
        {1,4,7,8},  // corner tet at v1
        {2,5,7,9},  // corner tet at v2
        {3,6,8,9},  // corner tet at v3
        {6,7,4,5},  // inner tet (diagonal 6-7), equatorial pair (4,5)
        {6,7,5,9},  // inner tet (diagonal 6-7), equatorial pair (5,9)
        {6,7,9,8},  // inner tet (diagonal 6-7), equatorial pair (9,8)
        {6,7,8,4},  // inner tet (diagonal 6-7), equatorial pair (8,4)
    };
};

template <typeMesh M> class InterfaceLevelSet_P2 : public Interface<M> {

    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    static const int nve = Rd::d;
    typedef FaceInterface<nve> Face;
    typedef SortArray<Ubyte, Element::Rd::d + 1> ElementIdx;

    std::vector<byte> ls_sign;
    std::vector<double> ls_;

    std::vector<Rd> outward_normal_;
    // Per-element P2 levelset values: ls_p2_[k][i] for i in 0..nv+ne-1.
    // Vertices i<nv use local vertex index; edge midpoints i>=nv use nv+local_edge.
    std::vector<std::vector<double>> ls_p2_;

  public:
    template <typeFunFEM Fct> InterfaceLevelSet_P2(const Mesh &MM, const Fct &lss, int label = 0);

    SignElement<Element> get_SignElement(int k) const override;
    Partition<Element> get_partition(int k) const override;
    Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k,
                                                         int ifac) const override;
    bool isCutFace(int k, int ifac) const override;
    bool isCut(int k) const override;

    R measure(const Face &f) const override;

    Rd normal(int k, std::span<double> x = std::span<double>()) const override { return outward_normal_[k]; }
    void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                       std::list<int> &erased_element, int sign_part) const override;

    Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const override;

    size_t size() const override { return this->faces_.size(); }

    R measure(int i) const override { return measure(this->faces_[i]); };

    double measure() const { return this->measure_; }

  private:
    // make_patch is a template to allow evaluating the P2 FunFEM directly.
    template <typeFunFEM Fct>
    void make_patch(const Fct &lss, int label);

    // Returns the local coarse-element edge index (0..ne-1) on which a
    // sub-element edge (g0,g1) lies, or Element::ne if the sub-edge is
    // interior to the coarse element (connects two edge midpoints).
    static Uint coarse_edge_of(int g0, int g1);
};

#include "interface_levelSetP2.tpp"

// template <typename M> void InterfaceLevelSet_P2<M>::make_patch(int label) {

//     assert(this->backMesh);
//     this->faces_.resize(0); // reinitialize arrays
//     this->vertices_.resize(0);
//     this->element_of_face_.resize(0);
//     this->outward_normal_.resize(0);
//     this->face_of_element_.clear();

//     const Mesh &Th = *(this->backMesh);
//     util::copy_levelset_sign(ls_, ls_sign);

//     const Uint nb_vertex_K = Element::nv;
//     double loc_ls[nb_vertex_K];
//     byte loc_ls_sign[nb_vertex_K];

//     for (int k = 0; k < this->backMesh->nbElmts(); k++) { // loop over elements

//         const typename Mesh::Element &K(Th[k]);

//         for (Uint i = 0; i < K.nv; ++i) {
//             loc_ls_sign[i] = ls_sign[Th(K[i])];
//             loc_ls[i]      = ls_[Th(K[i])];
//         }
//         const RefPatch<Element> &cut = RefPatch<Element>::instance(loc_ls_sign);
//         if (cut.empty())
//             continue;

//         for (typename RefPatch<Element>::const_face_iterator it = cut.face_begin(), end = cut.face_end(); it != end;
//              ++it) {
//             this->face_of_element_[k] = this->element_of_face_.size();
//             this->faces_.push_back(make_face(*it, K, loc_ls, label));
//             this->element_of_face_.push_back(k);
//             this->outward_normal_.push_back(make_normal(K, loc_ls));
//         }
//     }
// }

// template <typename M>
// const typename InterfaceLevelSet_P2<M>::Face
// InterfaceLevelSet_P2<M>::make_face(const typename RefPatch<Element>::FaceIdx &ref_tri, const typename Mesh::Element &K,
//                                    const double lset[Element::nv], int label) {

//     Uint loc_vert_num;
//     Uint triIdx[nve];

//     for (Uint j = 0; j < nve; ++j) {
//         loc_vert_num = ref_tri[j];

//         if (loc_vert_num < K.nv) { // zero vertex
//             // const Uint idx = (*this->backMesh)(K[loc_vert_num]);
//             // Rd Q           = (*this->backMesh)(idx);
//             // this->vertices_.push_back(Q);
//             // triIdx[j] = this->vertices_.size() - 1;
//             std::cout << " Interface cutting through a node " << std::endl;
//             exit(EXIT_FAILURE);
//             // assert(0);
//         } else { // genuine edge vertex

//             const Ubyte i0 = Mesh::Element::nvedge[loc_vert_num - K.nv][0],
//                         i1 = Mesh::Element::nvedge[loc_vert_num - K.nv][1];

//             const double t = lset[i0] / (lset[i0] - lset[i1]);
//             Rd Q           = (1.0 - t) * ((Rd)K[i0]) + t * ((Rd)K[i1]); // linear interpolation
//             this->vertices_.push_back(Q);
//             triIdx[j] = this->vertices_.size() - 1;

//             this->edge_of_node_.push_back(loc_vert_num - K.nv);
//         }
//     }
//     return Face(triIdx, label);
// }

// template <typename M>
// typename InterfaceLevelSet_P2<M>::Rd InterfaceLevelSet_P2<M>::make_normal(const typename Mesh::Element &K,
//                                                                          const double lset[Element::nv]) {

//     Rd grad[Element::nv];
//     K.Gradlambda(grad);
//     Rd normal_ls;
//     for (int i = 0; i < Mesh::Element::nv; ++i) {
//         normal_ls += grad[i] * lset[i];
//     }
//     normal_ls /= normal_ls.norm();
//     return normal_ls;
// }

#endif
