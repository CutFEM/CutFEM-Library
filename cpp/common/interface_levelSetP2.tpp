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

// ---------------------------------------------------------------------------
// get_reference_P2_node<M>(i)
// Returns the reference-element coordinates of the i-th P2 node for mesh M.
// Node ordering matches the SubElementTable and nvedge conventions:
//   2D Triangle (nv=3, ne=3):
//     0=(0,0)  1=(1,0)  2=(0,1)
//     3=mid(v1,v2)=(0.5,0.5)  4=mid(v2,v0)=(0,0.5)  5=mid(v0,v1)=(0.5,0)
//   3D Tet (nv=4, ne=6):
//     0..3 = vertices, 4..9 = edge midpoints (order matches nvedge)
// ---------------------------------------------------------------------------
template <typeMesh M>
typename M::Element::RdHat
get_reference_P2_node(int i) {
    using RdHat = typename M::Element::RdHat;
    if constexpr (M::D == 3) {
        static const RdHat ref_nodes[10] = {
            RdHat(0,0,0), RdHat(1,0,0), RdHat(0,1,0), RdHat(0,0,1), // 0-3: vertices
            RdHat(0.5,0,0), RdHat(0,0.5,0), RdHat(0,0,0.5),          // 4-6: mid(v0,v1/v2/v3)
            RdHat(0.5,0.5,0), RdHat(0.5,0,0.5), RdHat(0,0.5,0.5)     // 7-9: mid(v1,v2/v3), mid(v2,v3)
        };
        return ref_nodes[i];
    } else {
        // nvedge[0]={1,2}, [1]={2,0}, [2]={0,1}  →  node 3=mid(v1,v2), 4=mid(v2,v0), 5=mid(v0,v1)
        static const RdHat ref_nodes[6] = {
            RdHat(0,0), RdHat(1,0), RdHat(0,1),            // 0-2: vertices
            RdHat(0.5,0.5), RdHat(0,0.5), RdHat(0.5,0)     // 3=mid(v1,v2), 4=mid(v2,v0), 5=mid(v0,v1)
        };
        return ref_nodes[i];
    }
}

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
    // Use the P2 LS values at coarse vertices stored during make_patch.
    // ls_sign[iglb] would be wrong here: for a P2 function, lss.array() is
    // indexed by P2 DOF indices, not mesh vertex indices, so ls_sign[vertex_idx]
    // accesses an arbitrary DOF entry and gives the wrong sign.
    for (int i = 0; i < Element::nv; ++i) {
        loc_ls[i] = util::sign(ls_p2_[k][i]);
    }
    return SignElement<Element>(loc_ls);
}

template <typeMesh M> bool InterfaceLevelSet_P2<M>::isCut(int k) const {
    return (this->face_of_element_.find(k) != this->face_of_element_.end());
}

// Check if face ifac of element k is cut by the P2 interface.
// Uses the stored per-element P2 levelset values (ls_p2_[k]).
template <typeMesh M> bool InterfaceLevelSet_P2<M>::isCutFace(int k, int ifac) const {
    double max_val = std::numeric_limits<double>::lowest();
    double min_val = std::numeric_limits<double>::max();

    // Check coarse vertex signs on the face.
    for (int lv : Element::nvhyperFace[ifac]) {
        const double v = ls_p2_[k][lv];
        if (v > max_val) max_val = v;
        if (v < min_val) min_val = v;
    }
    // Check edge-midpoint signs on the face.
    if constexpr (M::D == 2) {
        // In 2D, face ifac IS edge ifac, midpoint at local P2 index nv + ifac.
        const double v = ls_p2_[k][Element::nv + ifac];
        if (v > max_val) max_val = v;
        if (v < min_val) min_val = v;
    } else {
        // In 3D, edgeOfFace[ifac] lists the coarse-element edge indices on the face.
        for (int le : Element::edgeOfFace[ifac]) {
            const double v = ls_p2_[k][Element::nv + le];
            if (v > max_val) max_val = v;
            if (v < min_val) min_val = v;
        }
    }
    return (max_val > 0.0 && min_val < 0.0);
}

// Bulk partition for element k, using the P2 sub-element refinement.
// The coarse element is split into n_sub sub-elements (SubElementTable), each
// cut with a P1 levelset. Resulting simplices are mapped back to coarse reference
// space and collected into a custom Partition.
template <typeMesh M>
Partition<typename InterfaceLevelSet_P2<M>::Element> InterfaceLevelSet_P2<M>::get_partition(int k) const {
    using Element = typename InterfaceLevelSet_P2<M>::Element;
    using Table   = SubElementTable<Element>;
    using RdHat   = typename Element::RdHat;
    constexpr int n_nodes = Element::nv + Element::ne;

    // Retrieve the stored per-element P2 levelset values.
    const double* vals = ls_p2_[k].data();

    // Build the custom (empty) partition for the coarse element.
    Partition<Element> coarse_partition((*this->backMesh)[k], true);

    for (int s = 0; s < Table::n_sub; ++s) {
        // P2 LS values at the sub-element's nv corners (from global P2 nodes).
        double sub_vals[Element::nv];
        for (int j = 0; j < Element::nv; ++j)
            sub_vals[j] = vals[Table::idx[s][j]];

        // P1 cut of the sub-element using its corner LS values.
        Partition<Element> sub_part((*this->backMesh)[k], sub_vals);

        // Reference coordinates of this sub-element's corners in the COARSE ref space.
        RdHat sub_corners[Element::nv];
        for (int j = 0; j < Element::nv; ++j)
            sub_corners[j] = get_reference_P2_node<M>(Table::idx[s][j]);

        // Iterate over all simplices produced by the P1 cut (both signs).
        for (auto it = sub_part.element_begin(0); it != sub_part.element_end(0); ++it) {
            const int sign = sub_part.whatSign(it);

            // Compute coarse reference coordinates of each simplex vertex directly.
            // sub_corners[j] already holds the coarse ref coords of the j-th sub-vertex.
            RdHat mapped_pts[Element::nv];
            for (int i = 0; i < Element::nv; ++i) {
                const Uint idx = (*it)[i];
                if (idx < static_cast<Uint>(Element::nv)) {
                    // Corner vertex: coarse ref coords are sub_corners[idx].
                    mapped_pts[i] = sub_corners[idx];
                } else {
                    // Cut point on local sub-element edge e_loc at parameter t.
                    const int e_loc = static_cast<int>(idx) - Element::nv;
                    const int v0    = Element::nvedge[e_loc][0];
                    const int v1    = Element::nvedge[e_loc][1];
                    const double t  = -sub_vals[v0] / (sub_vals[v1] - sub_vals[v0]);
                    mapped_pts[i]   = (1.0 - t) * sub_corners[v0] + t * sub_corners[v1];
                }
            }
            coarse_partition.add_simplex(sign, mapped_pts);
        }
    }
    return coarse_partition;
}

// Face partition for a cut face, using P2 sub-triangle refinement (3D) or
// P1 fallback (2D).
template <typeMesh M>
Partition<typename InterfaceLevelSet_P2<M>::Element::Face>
InterfaceLevelSet_P2<M>::get_partition_face(const typename Element::Face &face, int k, int ifac) const {
    using Element  = typename InterfaceLevelSet_P2<M>::Element;
    using FaceType = typename Element::Face;
    using RdHatFace = typename FaceType::RdHat;

    if constexpr (M::D == 2) {
        // 2D: face (edge) has nv=2 vertices; use P1 partition with coarse vertex signs.
        double loc_ls[FaceType::nv];
        for (int i = 0; i < FaceType::nv; ++i)
            loc_ls[i] = ls_p2_[k][Element::nvhyperFace[ifac][i]];
        return Partition<FaceType>(face, loc_ls);
    } else {
        // 3D: face is a triangle with 3 vertices + 3 edge midpoints.
        // Use SubElementTable<Triangle2> (same index structure as any triangle face).
        using FaceTable = SubElementTable<Triangle2>;

        // Build face_vals[6] in FaceTable ordering (nodes 0,1,2 = face vertices;
        // nodes 3,4,5 = midpoints per Triangle2::nvedge ordering).
        double face_vals[6];
        for (int i = 0; i < 3; ++i)
            face_vals[i] = ls_p2_[k][Element::nvhyperFace[ifac][i]];

        // For each face edge e (in Triangle2/FaceType nvedge ordering), find the
        // corresponding coarse element edge and read its midpoint LS value.
        for (int e = 0; e < 3; ++e) {
            const int fv0 = FaceType::nvedge[e][0];
            const int fv1 = FaceType::nvedge[e][1];
            const int cv0 = Element::nvhyperFace[ifac][fv0];
            const int cv1 = Element::nvhyperFace[ifac][fv1];
            int ce = -1;
            for (int j = 0; j < Element::ne; ++j) {
                if ((Element::nvedge[j][0] == cv0 && Element::nvedge[j][1] == cv1) ||
                    (Element::nvedge[j][0] == cv1 && Element::nvedge[j][1] == cv0)) {
                    ce = j; break;
                }
            }
            assert(ce >= 0);
            face_vals[3 + e] = ls_p2_[k][Element::nv + ce];
        }

        // Custom partition for the coarse face.
        Partition<FaceType> coarse_face_partition(face, true);

        for (int s = 0; s < FaceTable::n_sub; ++s) {
            double sub_vals[FaceType::nv];
            for (int j = 0; j < FaceType::nv; ++j)
                sub_vals[j] = face_vals[FaceTable::idx[s][j]];

            Partition<FaceType> sub_part(face, sub_vals);

            // Reference coords of sub-triangle corners in the coarse face ref space.
            RdHatFace sub_corners[FaceType::nv];
            for (int j = 0; j < FaceType::nv; ++j)
                sub_corners[j] = get_reference_P2_node<Mesh2>(FaceTable::idx[s][j]);

            for (auto it = sub_part.element_begin(0); it != sub_part.element_end(0); ++it) {
                const int sign = sub_part.whatSign(it);

                // Compute coarse face reference coordinates directly.
                RdHatFace mapped_pts[FaceType::nv];
                for (int i = 0; i < FaceType::nv; ++i) {
                    const Uint idx = (*it)[i];
                    if (idx < static_cast<Uint>(FaceType::nv)) {
                        // Corner vertex: coarse face ref coords are sub_corners[idx].
                        mapped_pts[i] = sub_corners[idx];
                    } else {
                        // Cut point on local sub-face edge e_loc at parameter t.
                        const int e_loc = static_cast<int>(idx) - FaceType::nv;
                        const int v0    = FaceType::nvedge[e_loc][0];
                        const int v1    = FaceType::nvedge[e_loc][1];
                        const double t  = -sub_vals[v0] / (sub_vals[v1] - sub_vals[v0]);
                        mapped_pts[i]   = (1.0 - t) * sub_corners[v0] + t * sub_corners[v1];
                    }
                }
                coarse_face_partition.add_simplex(sign, mapped_pts);
            }
        }
        return coarse_face_partition;
    }
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
    // Both endpoints are edge midpoints of different coarse edges: interior sub-edge.
    if (g0 >= E::nv && g1 >= E::nv) return static_cast<Uint>(E::ne);
    // Exactly one endpoint is an edge midpoint: the sub-edge lies on that coarse edge.
    if (g0 >= E::nv) return static_cast<Uint>(g0 - E::nv);
    if (g1 >= E::nv) return static_cast<Uint>(g1 - E::nv);
    // Both endpoints are original coarse vertices: find the matching coarse edge.
    for (int e = 0; e < E::ne; ++e) {
        if ((E::nvedge[e][0] == g0 && E::nvedge[e][1] == g1) ||
            (E::nvedge[e][0] == g1 && E::nvedge[e][1] == g0))
            return static_cast<Uint>(e);
    }
    return static_cast<Uint>(E::ne); // should be unreachable
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

    // Allocate per-element P2 levelset storage.
    ls_p2_.assign(Th.nbElmts(), std::vector<double>(n_nodes));

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

        // Store P2 levelset values for this element (used by get_partition etc.).
        ls_p2_[k].assign(vals, vals + n_nodes);

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
                Rd   face_v[nve];
                for (Uint j = 0; j < static_cast<Uint>(nve); ++j) {
                    const Uint loc = (*it)[j];
                    if (loc < static_cast<Uint>(Element::nv)) {
                        std::cout << " Interface cutting through a P2 sub-node " << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    // Identify which sub-element edge is cut (loc = nv + local_edge).
                    const int e_sub = static_cast<int>(loc) - Element::nv;
                    const int lv0   = Element::nvedge[e_sub][0];
                    const int lv1   = Element::nvedge[e_sub][1];
                    const int g0    = Table::idx[s][lv0];
                    const int g1    = Table::idx[s][lv1];
                    const double v0 = vals[g0], v1 = vals[g1];
                    const double t  = v0 / (v0 - v1);
                    const Rd Q      = (1.0 - t) * nodes[g0] + t * nodes[g1];
                    this->vertices_.push_back(Q);
                    triIdx[j] = static_cast<Uint>(this->vertices_.size() - 1);
                    face_v[j] = Q;
                    this->edge_of_node_.push_back(coarse_edge_of(g0, g1));
                }

                this->face_of_element_.emplace(k, this->element_of_face_.size());
                this->faces_.push_back(Face(triIdx, label));
                this->element_of_face_.push_back(k);

                // Compute the normal using the SORTED vertex order stored in faces_.back().
                // Face(triIdx, label) is a SortArray and reorders the indices; the normal
                // must be computed from the same order used by mapToPhysicalFace, otherwise
                // an odd-permutation sort silently flips the triangle handedness.
                Rd sv[nve];
                for (int j = 0; j < nve; ++j)
                    sv[j] = this->vertices_[this->faces_.back()[j]];

                Rd normal_ls;
                if constexpr (M::D == 2) {
                    normal_ls = (sv[1] - sv[0]).perp();
                } else {
                    const Rd e1 = sv[1] - sv[0];
                    const Rd e2 = sv[2] - sv[0];
                    normal_ls   = e1 ^ e2;
                }
                // Orient toward the positive-levelset side of the sub-element.
                for (int j = 0; j < Element::nv; ++j) {
                    if (sub_vals[j] > 0.0) {
                        if ((nodes[Table::idx[s][j]] - sv[0], normal_ls) < 0.0)
                            normal_ls = -normal_ls;
                        break;
                    }
                }
                normal_ls /= normal_ls.norm();
                this->outward_normal_.push_back(normal_ls);
                this->measure_ += measure(this->faces_.back());
            }
        }
    }
}

#endif
