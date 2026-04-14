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
#ifndef COMMON_LEVELSET_INTERFACEP2_TPP
#define COMMON_LEVELSET_INTERFACEP2_TPP

template <typeMesh M>
template <typeFunFEM Fct>
InterfaceLevelSet_P2<M>::InterfaceLevelSet_P2(const M &MM, const Fct &lss, int label)
    : Interface<M>(MM), ls_(lss.array().begin(), lss.array().end()) {

    make_patch(lss, label);
}

template <typeMesh M>
SignElement<typename InterfaceLevelSet_P2<M>::Element> InterfaceLevelSet_P2<M>::get_SignElement(int k) const {
    typedef typename InterfaceLevelSet_P2<M>::Element Element;
    byte loc_ls[Element::nv];
    for (int i = 0; i < Element::nv; ++i) {
        int iglb  = this->backMesh->at(k, i);
        loc_ls[i] = ls_sign[iglb];
    }
    return SignElement<Element>(loc_ls);
}

template <typeMesh M> bool InterfaceLevelSet_P2<M>::isCut(int k) const {
    return (this->face_of_element_.find(k) != this->face_of_element_.end());
}

// Uses coarse-vertex signs from ls_ (same as P1 InterfaceLevelSet).
// This is conservative: may return false when only an edge midpoint crosses
// zero, but is consistent with get_SignElement and get_partition.
template <typeMesh M> bool InterfaceLevelSet_P2<M>::isCutFace(int k, int ifac) const {
    if constexpr (M::D == 3) {
        for (int e = 0; e < Element::Face::ne; ++e) {
            const int id_edge = Element::edgeOfFace[ifac][e];
            const int i1      = Element::nvedge[id_edge][0];
            const int i2      = Element::nvedge[id_edge][1];
            if (ls_[this->backMesh->at(k, i1)] * ls_[this->backMesh->at(k, i2)] < 0)
                return true;
        }
        return false;
    } else {
        const int i1 = Element::nvedge[ifac][0];
        const int i2 = Element::nvedge[ifac][1];
        return ls_[this->backMesh->at(k, i1)] * ls_[this->backMesh->at(k, i2)] < 0;
    }
}

template <typeMesh M>
Partition<typename InterfaceLevelSet_P2<M>::Element> InterfaceLevelSet_P2<M>::get_partition(int k) const {
    typedef typename InterfaceLevelSet_P2<M>::Element Element;

    double loc_ls[Element::nv];
    for (int i = 0; i < Element::nv; ++i) {
        int iglb  = this->backMesh->at(k, i);
        loc_ls[i] = ls_[iglb];
    }

    return Partition<Element>((*this->backMesh)[k], loc_ls);
}

template <typeMesh M>
Partition<typename InterfaceLevelSet_P2<M>::Element::Face>
InterfaceLevelSet_P2<M>::get_partition_face(const typename Element::Face &face, int k, int ifac) const {
    typedef typename InterfaceLevelSet_P2<M>::Element Element;

    double loc_ls[Element::Face::nv];
    for (int i = 0; i < Element::Face::nv; ++i) {
        int j     = Element::nvhyperFace[ifac][i];
        int iglb  = this->backMesh->at(k, j);
        loc_ls[i] = ls_[iglb];
    }
    return Partition<typename Element::Face>(face, loc_ls);
}

template <typeMesh M>
void InterfaceLevelSet_P2<M>::cut_partition(Physical_Partition<typename InterfaceLevelSet_P2<M>::Element> &local_partition,
                                         std::vector<ElementIdx> &new_element_idx, std::list<int> &erased_element,
                                         int sign_part) const {
    std::cout << " An element might be cut multiplue time, and it is not "
                 "suppose to happen"
              << std::endl;
    exit(EXIT_FAILURE);
};

template <typeMesh M> R InterfaceLevelSet_P2<M>::measure(const Face &f) const {
    Rd l[nve];
    for (int i = 0; i < nve; ++i)
        l[i] = this->vertices_[f[i]];
    return geometry::measure_hyper_simplex(l);
};

// Rd get_intersection_node(int k, const Rd A, const Rd B) const {
//   double fA = fun.eval(k, A);
//   double fB = fun.eval(k, B);
//   double t = -fA/(fB-fA);
//   return (1-t) * A + t * B;
// }

template <typeMesh M>
typename InterfaceLevelSet_P2<M>::Rd
InterfaceLevelSet_P2<M>::mapToPhysicalFace(int ifac, const typename InterfaceLevelSet_P2<M>::Element::RdHatBord x) const {
    typename InterfaceLevelSet_P2<M>::Rd N[nve];
    for (int i = 0; i < nve; ++i)
        N[i] = this->vertices_[this->faces_[ifac][i]];
    return geometry::map_point_to_simplex(N, x);
}

// Returns the coarse element edge index (0..ne-1) on which the sub-edge
// (g0, g1) lies, or Element::ne if the sub-edge is interior (both endpoints
// are edge midpoints of different coarse edges).
template <typeMesh M>
Uint InterfaceLevelSet_P2<M>::coarse_edge_of(int g0, int g1) {
    using E = typename M::Element;
    // If one endpoint is an edge midpoint, the sub-edge lies on that coarse edge
    // (provided the other endpoint is one of that coarse edge's vertices).
    if (g0 >= E::nv) return static_cast<Uint>(g0 - E::nv);
    if (g1 >= E::nv) return static_cast<Uint>(g1 - E::nv);
    // Both endpoints are original vertices: find the matching coarse edge.
    for (int e = 0; e < E::ne; ++e) {
        if ((E::nvedge[e][0] == g0 && E::nvedge[e][1] == g1) ||
            (E::nvedge[e][0] == g1 && E::nvedge[e][1] == g0))
            return static_cast<Uint>(e);
    }
    // Both are edge midpoints of different coarse edges: interior sub-edge.
    return static_cast<Uint>(E::ne);
}

template <typeMesh M>
template <typeFunFEM Fct>
void InterfaceLevelSet_P2<M>::make_patch(const Fct &lss, int label) {

    using Table = SubElementTable<Element>;
    constexpr int n_nodes = Element::nv + Element::ne;

    assert(this->backMesh);
    this->faces_.resize(0);
    this->vertices_.resize(0);
    this->element_of_face_.resize(0);
    this->outward_normal_.resize(0);
    this->face_of_element_.clear();
    this->measure_ = 0.;

    const M &Th = *(this->backMesh);
    util::copy_levelset_sign(ls_, ls_sign);

    for (int k = 0; k < Th.nbElmts(); ++k) {
        const Element &K(Th[k]);

        // Compute physical positions and P2 levelset values at all node positions.
        // Indices 0..nv-1 are coarse vertices; nv..nv+ne-1 are edge midpoints.
        Rd     nodes[n_nodes];
        double vals[n_nodes];

        for (int i = 0; i < Element::nv; ++i) {
            nodes[i] = (Rd)K[i];
            vals[i]  = lss.eval(k, nodes[i], 0, op_id);
        }
        for (int e = 0; e < Element::ne; ++e) {
            const int i0           = Element::nvedge[e][0];
            const int i1           = Element::nvedge[e][1];
            nodes[Element::nv + e] = 0.5 * (nodes[i0] + nodes[i1]);
            vals[Element::nv + e]  = lss.eval(k, nodes[Element::nv + e], 0, op_id);
        }

        // Process each sub-element.
        for (int s = 0; s < Table::n_sub; ++s) {
            byte   sub_signs[Element::nv];
            double sub_vals[Element::nv];
            for (int j = 0; j < Element::nv; ++j) {
                sub_vals[j]  = vals[Table::idx[s][j]];
                sub_signs[j] = util::sign(sub_vals[j]);
            }

            const RefPatch<Element> &cut = RefPatch<Element>::instance(sub_signs);
            if (cut.empty())
                continue;

            for (auto it = cut.face_begin(); it != cut.face_end(); ++it) {
                Uint triIdx[nve];
                for (Uint j = 0; j < static_cast<Uint>(nve); ++j) {
                    const Uint loc = (*it)[j];
                    if (loc < static_cast<Uint>(Element::nv)) {
                        std::cout << " Interface cutting through a P2 sub-node " << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    // Identify which sub-element edge is cut (loc = nv + local_edge).
                    const int e_sub = static_cast<int>(loc) - Element::nv;
                    const int lv0   = Element::nvedge[e_sub][0]; // local sub-vertex index
                    const int lv1   = Element::nvedge[e_sub][1];
                    const int g0    = Table::idx[s][lv0];        // global P2 node index
                    const int g1    = Table::idx[s][lv1];
                    const double v0 = vals[g0], v1 = vals[g1];
                    const double t  = v0 / (v0 - v1);
                    const Rd Q      = (1.0 - t) * nodes[g0] + t * nodes[g1];
                    this->vertices_.push_back(Q);
                    triIdx[j] = static_cast<Uint>(this->vertices_.size() - 1);
                    // edge_of_node_: coarse edge index, or Element::ne if interior.
                    this->edge_of_node_.push_back(coarse_edge_of(g0, g1));
                }

                this->face_of_element_[k] = this->element_of_face_.size();
                this->faces_.push_back(Face(triIdx, label));
                this->element_of_face_.push_back(k);

                // Normal: gradient of the P2 levelset at the sub-element centroid.
                Rd centroid{};
                for (int j = 0; j < Element::nv; ++j)
                    centroid += (1.0 / Element::nv) * nodes[Table::idx[s][j]];

                Rd normal_ls;
                if constexpr (M::D == 2) {
                    normal_ls = Rd(lss.eval(k, centroid, 0, op_dx),
                                   lss.eval(k, centroid, 0, op_dy));
                } else {
                    normal_ls = Rd(lss.eval(k, centroid, 0, op_dx),
                                   lss.eval(k, centroid, 0, op_dy),
                                   lss.eval(k, centroid, 0, op_dz));
                }
                normal_ls /= normal_ls.norm();
                this->outward_normal_.push_back(normal_ls);
                this->measure_ += measure(this->faces_.back());
            }
        }
    }
}

#endif
