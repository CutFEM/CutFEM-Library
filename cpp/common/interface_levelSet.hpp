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
#ifndef COMMON_LEVELSET_INTERFACE_HPP
#define COMMON_LEVELSET_INTERFACE_HPP

#include "base_interface.hpp"

template <typeMesh M> class InterfaceLevelSet : public Interface<M> {

    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    static const int nve = Rd::d;
    typedef FaceInterface<nve> Face;
    typedef SortArray<Ubyte, Element::Rd::d + 1> ElementIdx;

    KN<byte> ls_sign;
    KN<double> ls_;

    std::vector<Rd> outward_normal_;

  public:
    template <typeFunFEM Fct> InterfaceLevelSet(const Mesh &MM, const Fct &lss, int label = 0);

    SignElement<Element> get_SignElement(int k) const override;
    Partition<Element> get_partition(int k) const override;
    Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k,
                                                         int ifac) const override;
    // return Partition<Element>((*this->backMesh)[k], loc_ls);
    // Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k, int ifac) const;
    // {
    //         double loc_ls[Element::Face::nv];
    //         for (int i = 0; i < Element::Face::nv; ++i) {
    //             int j     = Element::nvhyperFace[ifac][i];
    //             int iglb  = this->backMesh->at(k, j);
    //             loc_ls[i] = ls_[iglb];
    //         }
    //         return Partition<typename Element::Face>(face, loc_ls);
    // }
    bool isCutFace(int k, int ifac) const override;
    bool isCut(int k) const override;

    R measure(const Face &f) const override;

    Rd normal(int k, std::span<double> x = std::span<double>()) const override { return outward_normal_[k]; }
    void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                       std::list<int> &erased_element, int sign_part) const override;

    Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const override;

    size_t size() const override { return this->faces_.size(); }

    R measure(int i) const override { return measure(this->faces_[i]); };

  private:
    void make_patch(int label);

    const Face make_face(const typename RefPatch<Element>::FaceIdx &ref_tri, const typename Mesh::Element &K,
                         const double lset[Element::nv], int label);

    Rd make_normal(const typename Mesh::Element &K, const double lset[Element::nv]);
};

#include "interface_levelSet.tpp"

// template <typename M> void InterfaceLevelSet<M>::make_patch(int label) {

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
// const typename InterfaceLevelSet<M>::Face
// InterfaceLevelSet<M>::make_face(const typename RefPatch<Element>::FaceIdx &ref_tri, const typename Mesh::Element &K,
//                                 const double lset[Element::nv], int label) {

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
// typename InterfaceLevelSet<M>::Rd InterfaceLevelSet<M>::make_normal(const typename Mesh::Element &K,
//                                                                     const double lset[Element::nv]) {

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
