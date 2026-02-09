#ifndef CUTFEM_PROBLEM_BOUNDARY_HPP
#define CUTFEM_PROBLEM_BOUNDARY_HPP

#include <vector>
#include "../common/point.hpp"
#include "../FESpace/FESpace.hpp"


template <int D> class DofData {

    using v_t = typename typeRd<D>::Rd;

  public:
    DofData() = default;

    DofData(int k, int dom, int ci, v_t P) : k(k), domain(dom), ci(ci), P(P) {}

    DofData(DofData &&)            = default;
    DofData &operator=(DofData &&) = default;


    int k;
    
    int domain{0};
    
    int ci;
    
    v_t P;
};

template <typename M> class BoundaryDirichlet {

    using Rd           = M::Rd;
    static const int D = Rd::d;
    using mesh_t       = M;
    using elt_t        = typename mesh_t::Element;
    using space_t      = GFESpace<mesh_t>;
    using cutspace_t   = CutFESpace<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;
    using dof_data_t   = DofData<D>;

  public:
    BoundaryDirichlet(const cutspace_t &Vh, std::vector<int> lab = {});  // set strong BC on outer (fitted) boundary of CutFEM problem
    BoundaryDirichlet(const cutspace_t &Vh, const int domain = 0);      // set strong BC on outer boundary of the active mesh corresponding to domain "domain" of a CutFEM problem
    BoundaryDirichlet(const cutspace_t &Vh, const BarycentricActiveMesh2& active_mesh, const int domain = 0);
    BoundaryDirichlet(const space_t &Vh, std::list<int> lab = {});     // set strong BC on the boundary of a standard fitted FEM mesh

    void apply_inhomogeneous(std::map<std::pair<int, int>, double> &A);
    template <typename Fct>
    void apply(std::span<double> b, const Fct &f);      // exact function f
    void apply(std::span<double> b, const fct_t &f);    // fem function f
    void apply(std::span<double> b, double val);
    // Column elimination + RHS shift for inhomogeneous g:
    // for r!=I: b[r] -= A(r,I)*g, and set A(r,I)=0. Keeps A(I,I)=1.
    void finalize_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                                std::span<double> b,
                                size_t dof_start = 0);

    // void apply_inhomogeneous(std::map<std::pair<int, int>, double> &A, std::span<double> b);
    // void apply_homogeneous(std::map<std::pair<int, int>, double> &A, std::span<double> b);

    std::map<int, dof_data_t> boundary_dofs;
};

// template <typename M>
// BoundaryDirichlet<M>::BoundaryDirichlet(const cutspace_t &Vh, std::vector<int> lab) {
//     const auto &Th    = Vh.Th;
//     const auto &cutTh = Vh.get_mesh();

//     // loop over boundary elements
//     for (int k = Th.first_boundary_element(); k < Th.last_boundary_element(); k += Th.next_boundary_element()) {

//         const auto &BE(Th.be(k));
//         auto it_lab = std::find(lab.begin(), lab.end(), BE.lab);

//         if (it_lab == lab.end())
//             continue;

//         auto [elt_idx, face_idx] = Th.getBoundaryElement(k);
//         std::vector<int> idxK = Vh.idxAllElementFromBackMesh(elt_idx, -1); // index of element elt_idx in the cut space
//         assert(idxK.size() == 1);                                          // asserts that the boundary is not cut

//         const auto &T(Th[elt_idx]);
//         // const auto FK(Vh[elt_idx]);
//         const auto &FK(Vh[idxK[0]]);

//         const int domain = FK.get_domain();

//         std::vector<Rd> dof_point(FK.tfe->NbPtforInterpolation);
//         FK.tfe->global_dofs(T, dof_point);

//         // int ndof_node = FK.getBasisFct()->getDofEntities()[0];
//         // int ndof_edge = FK.getBasisFct().getDofEntities()[1];
//         int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

//         for (int ic = 0; ic < Vh.N; ++ic) {
//             for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

//                 // here need a fct that give the item from the dof for P3 for example
//                 int id_item       = df_loc;
//                 bool is_on_border = false;

//                 // case 1: dof is on node. E.g. (0, 1, 2) if triangle
//                 if (id_item < T.nv) {
//                     for (int i = 0; i < elt_t::nva; ++i) {
//                         int i_e = elt_t::nvedge.at(face_idx).at(i);
//                         if (i_e == id_item) {
//                             is_on_border = true;
//                             break;
//                         }
//                     }
//                 }
//                 // case 2: dof is on an edge (or face)
//                 else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

//                     int id_face = (id_item - T.nv) / ndof_edge_per_component;
//                     // std::cout << "id_item = " << id_item << std::endl;
//                     // std::cout << "id_item - T.nv = " << id_item - T.nv << std::endl;
//                     // std::cout << "id_face = " << id_face << std::endl;
//                     // std::cout << "ndof_edge = " << ndof_edge_per_component << std::endl;
//                     if (id_face == face_idx) {
//                         is_on_border = true;
//                     }
//                 }
//                 if (is_on_border) {
//                     Rd P                   = dof_point.at(df_loc);
//                     size_t df_glob         = FK.loc2glb(df);
//                     boundary_dofs[df_glob] = DofData<D>(elt_idx, domain, ic, P);
//                 }
//             }
//         }
//     }
// }

template <typename M>
BoundaryDirichlet<M>::BoundaryDirichlet(const cutspace_t &Vh, std::vector<int> lab) {
    const auto &Th    = Vh.Th;
    const auto &cutTh = Vh.get_mesh();

    // loop over boundary elements
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element(); idx_be += cutTh.next_boundary_element()) {

        int idx_bdry_face; // Index of the boundary face in the boundary element
        const int kb = cutTh.Th.BoundaryElement(idx_be, idx_bdry_face);  // Get the boundary element index in the original mesh

        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1); // Get all the indices of the active elements corresponding to the boundary element in the cut mesh

        // const auto &BE(Th.be(k));
        // auto it_lab = std::find(lab.begin(), lab.end(), BE.lab);

        // if (it_lab == lab.end())
        //     continue;

        // auto [elt_idx, face_idx] = Th.getBoundaryElement(k);
        // std::vector<int> idxK = Vh.idxAllElementFromBackMesh(elt_idx, -1); // index of element elt_idx in the cut space
        // assert(idxK.size() == 1);                                          // asserts that the boundary is not cut
        
        // assert(idxK.size() == 1);
        if (idxK.size() == 0) continue;

        int k = idxK[0];
        const auto &FK(Vh[k]);  // Get the finite element for the current element in the finite element space
        const int domain = FK.get_domain();
        
        // Get the current element and its boundary element
        const auto &T(cutTh.Th[kb]);
        const auto &BE(cutTh.be(idx_be));

        // Check if the current boundary element label is in the list of labels to process, or if all labels are being processed
        auto it_lab = std::find(lab.begin(), lab.end(), BE.lab);
        if (it_lab == lab.end())
            continue;

        // Skip processing if the element is cut
        if (cutTh.isCut(k, 0)) {
            continue;
        }

        std::vector<Rd> dof_point(FK.tfe->NbPtforInterpolation);
        FK.tfe->global_dofs(T, dof_point);

        // int ndof_node = FK.getBasisFct()->getDofEntities()[0];
        // int ndof_edge = FK.getBasisFct().getDofEntities()[1];
        int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

        for (int ic = 0; ic < Vh.N; ++ic) {
            for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

                // here need a fct that give the item from the dof for P3 for example
                int id_item       = df_loc;
                bool is_on_border = false;

                // std::cout << "I'm here\n";

                // case 1: dof is on node. E.g. (0, 1, 2) if triangle
                if (id_item < T.nv) {
                    for (int i = 0; i < elt_t::nva; ++i) {
                        // std::cout << "nva = " << elt_t::nva << ", idx_bdry_face = " << idx_bdry_face << ", i = " << i << "\n";
                        int i_e = elt_t::nvedge.at(idx_bdry_face).at(i);
                        if (i_e == id_item) {
                            is_on_border = true;
                            break;
                        }
                    }
                }
                // case 2: dof is on an edge (or face)
                else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

                    int id_face = (id_item - T.nv) / ndof_edge_per_component;
                    // std::cout << "id_item = " << id_item << std::endl;
                    // std::cout << "id_item - T.nv = " << id_item - T.nv << std::endl;
                    // std::cout << "id_face = " << id_face << std::endl;
                    // std::cout << "ndof_edge = " << ndof_edge_per_component << std::endl;
                    if (id_face == idx_bdry_face) {
                        is_on_border = true;
                    }
                }
                if (is_on_border) {
                    Rd P                   = dof_point.at(df_loc);
                    size_t df_glob         = FK.loc2glb(df);
                    boundary_dofs[df_glob] = DofData<D>(kb, domain, ic, P);
                }
            }
        }
    }
}


// NOTE: This works only for nodal finite elements
template <typename M>
BoundaryDirichlet<M>::BoundaryDirichlet(const cutspace_t &Vh, const int domain) {
    const auto &Th    = Vh.Th;  // background mesh
    const auto &cutTh = Vh.get_mesh();  // active mesh

    // loop over boundary elements
    //for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element(); idx_be += cutTh.next_boundary_element()) {
    for (int k = cutTh.first_element(); k < cutTh.last_element(); k += cutTh.next_element()) {

        if (!cutTh.isCut(k, 0)) {
            continue;
        } 

        const int kb = cutTh.idxElementInBackMesh(k);
        const auto &FK(Vh[k]);  // Get the finite element for the current element in the finite element space
        const auto &T(Th[kb]);
        const int domain = FK.get_domain();

        for (int ifac = 0; ifac < M::Element::nea; ++ifac) { 

            int jfac = ifac;
            int kn   = cutTh.ElementAdj(k, jfac);

            // a neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn != -1)
                continue;
            
            std::vector<Rd> dof_points(FK.tfe->NbPtforInterpolation);
            FK.tfe->global_dofs(T, dof_points);

            int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;
        
            for (int ic = 0; ic < Vh.N; ++ic) {
                for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

                    // here need a fct that give the item from the dof for P3 for example
                    int id_item       = df_loc;
                    bool is_on_border = false;

                    // case 1: dof is on node. E.g. (0, 1, 2) if triangle
                    if (id_item < T.nv) {
                        for (int i = 0; i < elt_t::nva; ++i) {

                            int i_e = elt_t::nvedge.at(ifac).at(i);
                            if (i_e == id_item) {
                                is_on_border = true;
                                break;
                            }
                        }
                    }
                    // case 2: dof is on an edge (or face)
                    else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

                        int id_face = (id_item - T.nv) / ndof_edge_per_component;

                        if (id_face == ifac) {
                            is_on_border = true;
                        }
                    }
                    if (is_on_border) {

                        Rd P                   = dof_points.at(df_loc);
                        //std::cout << "Boundary point P = " << P << "\n";
                        size_t df_glob         = FK.loc2glb(df);
                        boundary_dofs[df_glob] = DofData<D>(kb, domain, ic, P);
                    }
                }
            }
            
        }

    }

}



template <typename M> BoundaryDirichlet<M>::BoundaryDirichlet(const space_t &Vh, std::list<int> lab) {
    const mesh_t &Th = Vh.Th;

    bool all_label = (lab.size() == 0);

    // loop over boundary elements
    for (int k = Th.first_boundary_element(); k < Th.last_boundary_element(); k += Th.next_boundary_element()) {

        const auto &BE(Th.be(k));

        if (!(util::contain(lab, BE.lab) || all_label)) {
            continue;
        }

        auto [elt_idx, face_idx] = Th.getBoundaryElement(k);
        const auto &T(Th[elt_idx]);
        const auto &FK(Vh[elt_idx]);

        const int domain = FK.get_domain();

        std::vector<Rd> dof_point(FK.tfe->NbPtforInterpolation);
        FK.tfe->global_dofs(T, dof_point);

        // std::cout << "typeid(*FK.tfe).name() = " << typeid(*FK.tfe).name() << std::endl;
        
        // std::vector<Rd> dof_point;
        // for (int p = 0; p < FK.tfe->NbPtforInterpolation; p++) {
        //     Rd P(FK.Pt(p));
        //     dof_point.push_back(P);
        // }
        // if (T.EdgeOrientation(0) < 0) {
        //     std::swap(dof_point[3], dof_point[4]); // 3,4
        // }
        // if (T.EdgeOrientation(1) < 0) {
        //     std::swap(dof_point[5], dof_point[6]); // 5,6
        // }
        // if (T.EdgeOrientation(2) < 0) {
        //     std::swap(dof_point[7], dof_point[8]); // 7,8
        // }


        int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

        for (int ic = 0; ic < Vh.N; ++ic) {
            for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

                // here need a fct that give the item from the dof for P3 for example
                int id_item       = df_loc;
                bool is_on_border = false;

                // case 1: dof is on node. E.g. (0, 1, 2) if triangle
                if (id_item < T.nv) {
                    for (int i = 0; i < elt_t::nva; ++i) {
                        int i_e = elt_t::nvedge.at(face_idx).at(i);
                        if (i_e == id_item) {
                            is_on_border = true;
                            break;
                        }
                    }
                }
                // case 2: dof is on an edge (or face)
                else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

                    int id_face = (id_item - T.nv) / ndof_edge_per_component;
                    if (id_face == face_idx) {
                        is_on_border = true;
                    }
                }
                if (is_on_border) {
                    Rd P                   = dof_point.at(df_loc);
                    size_t df_glob         = FK.loc2glb(df);
                    boundary_dofs[df_glob] = DofData<D>(elt_idx, domain, ic, P);
                }
            }
        }
    }
}

template <typename M> void BoundaryDirichlet<M>::apply_inhomogeneous(std::map<std::pair<int, int>, double> &A_map) {

    using it_t = std::map<std::pair<int, int>, double>::iterator;

    size_t nz_rm  = 0;
    it_t it_start = A_map.begin();
    for (auto &[df, dof_data] : boundary_dofs) {

        std::pair<int, int> min_val = std::make_pair(df, 0);
        std::pair<int, int> max_val = std::make_pair(df + 1, 0);

        auto it_begin = std::find_if(it_start, A_map.end(), [&min_val](auto &a) { return a.first >= min_val; });
        auto it_end   = std::find_if(it_begin, A_map.end(), [&max_val](auto &a) { return a.first >= max_val; });

        nz_rm += std::distance(it_begin, it_end);
        A_map.erase(it_begin, it_end);

        A_map[std::make_pair(df, df)] = 1.0;

        it_start = it_end;
        nz_rm--;
    }
    // std::cout << A_map.size() << " elements left" << std::endl;
    // std::cout << nz_rm << " elements removed" << std::endl;
    // std::cout << A_map.size() + nz_rm << " ?= " << nz << std::endl;
}

template <typename M> void BoundaryDirichlet<M>::apply(std::span<double> b, const fct_t &f) {

    for (auto &[df, dof_data] : boundary_dofs) {
        b[df] = f.evalOnBackMesh(dof_data.k, dof_data.domain, dof_data.P, dof_data.ci, op_id);
    }
}

template <typename M> 
template <typename Fct>
void BoundaryDirichlet<M>::apply(std::span<double> b, const Fct &f) {

    for (auto &[df, dof_data] : boundary_dofs) {
        b[df] = f(dof_data.P, dof_data.ci);
    }
}

template <typename M> void BoundaryDirichlet<M>::apply(std::span<double> b, double val) {

    for (auto &[df, dof_data] : boundary_dofs) {
        b[df] = val;
    }
}

// template <typename M>
// void BoundaryDirichlet<M>::apply_homogeneous(std::map<std::pair<int, int>, double> &A_map, std::span<double> b) {
//     for (const auto &dof : boundary_dofs) {
//         int i        = dof.first;
//         double value = dof.second.val;

//         // Erase the row
//         for (auto it = A_map.begin(); it != A_map.end();) {
//             if (it->first.first == i) {
//                 it = A_map.erase(it);
//             } else {
//                 ++it;
//             }
//         }

//         // Erase the column
//         for (auto it = A_map.begin(); it != A_map.end();) {
//             if (it->first.second == i) {
//                 it = A_map.erase(it);
//             } else {
//                 ++it;
//             }
//         }

//         // Set the diagonal element to 1
//         A_map[{i, i}] = 1.0;

//         // Update the RHS vector
//         b[i] = 0.;
//     }
// }

// template <typename M>
// void BoundaryDirichlet<M>::apply(std::map<std::pair<int, int>, double> &A_map,
//              std::span<double> b) {

//   using it_t = std::map<std::pair<int, int>, double>::iterator;

//   size_t N = space.NbDoF();

//   for (auto &[df, dof_data] : boundary_dofs) {
//     double g_x = dof_data.val;

//     // Subtract g(x) * i-th column from the RHS
//     for (size_t row = 0; row < N; ++row) {
//       auto it = A_map.find(std::make_pair(row, df));
//       if (it != A_map.end()) {
//         b[row] -= g_x * it->second;
//       }
//     }

//     // Zero out the i-th column in the matrix
//     for (size_t row = 0; row < N; ++row) {
//       auto it = A_map.find(std::make_pair(row, df));
//       if (it != A_map.end()) {
//         A_map.erase(it);
//       }
//     }

//     // Zero out the i-th row in the matrix
//     for (size_t col = 0; col < N; ++col) {
//       auto it = A_map.find(std::make_pair(df, col));
//       if (it != A_map.end()) {
//         A_map.erase(it);
//       }
//     }

//     // Set A[df, df] = 1
//     A_map[std::make_pair(df, df)] = 1.0;

//     // Set the i-th entry of the RHS to g(x)
//     b[df] = g_x;
//   }
// }

// template <typename M>
// void BoundaryDirichlet<M>::apply_homogeneous(std::map<std::pair<int, int>, double> &A_map,
//              std::span<double> b) {

//   using it_t = std::map<std::pair<int, int>, double>::iterator;

//   size_t N = space.NbDoF();

//   for (auto &[df, dof_data] : boundary_dofs) {

//     // Zero out the i-th column in the matrix
//     for (size_t row = 0; row < N; ++row) {
//       auto it = A_map.find(std::make_pair(row, df));
//       if (it != A_map.end()) {
//         A_map.erase(it);
//       }
//     }

//     // Zero out the i-th row in the matrix
//     for (size_t col = 0; col < N; ++col) {
//       auto it = A_map.find(std::make_pair(df, col));
//       if (it != A_map.end()) {
//         A_map.erase(it);
//       }
//     }

//     // Set A[df, df] = 1
//     A_map[std::make_pair(df, df)] = 1.0;

//     // Set the i-th entry of the RHS to zero
//     b[df] = 0.;
//   }
// }

// template <typename M>
// void BoundaryDirichlet<M>::apply_homogeneous(std::map<std::pair<int, int>, double> &A_map,
//                                                std::span<double> b) {
//     using it_t = std::map<std::pair<int, int>, double>::iterator;

//     size_t N = space.NbDoF();
//     size_t nz_rm_row = 0; // Track removed non-zero elements in rows
//     size_t nz_rm_col = 0; // Track removed non-zero elements in columns

//     it_t it_start = A_map.begin();

//     for (auto &[df, dof_data] : boundary_dofs) {
//         // Step 1: Remove all rows corresponding to the boundary DOF `df`
//         std::pair<int, int> min_val_row = std::make_pair(df, 0);
//         std::pair<int, int> max_val_row = std::make_pair(df + 1, 0);

//         auto it_begin_row = std::find_if(it_start, A_map.end(), [&min_val_row](auto &a) {
//             return a.first >= min_val_row;
//         });
//         auto it_end_row = std::find_if(it_begin_row, A_map.end(), [&max_val_row](auto &a) {
//             return a.first > max_val_row;
//         });

//         nz_rm_row += std::distance(it_begin_row, it_end_row);
//         A_map.erase(it_begin_row, it_end_row);

//         // Step 2: Remove column entries for boundary DOF `df`
//         for (auto it = A_map.begin(); it != A_map.end();) {
//             // Check if the column index matches `df`
//             if (it->first.second == df) {
//                 // Keep the diagonal element intact
//                 if (it->first.first != df) {
//                     it = A_map.erase(it); // Erase column entry
//                     ++nz_rm_col;          // Increment column removal counter
//                 } else {
//                     ++it; // Skip the diagonal entry
//                 }
//             } else {
//                 ++it;
//             }
//         }

//         // Step 3: Set the diagonal entry A(df, df) = 1
//         A_map[std::make_pair(df, df)] = 1.0;

//         // Step 4: Set the RHS value for the boundary DOF to `g(x)`
//         b[df] = dof_data.val;
//     }

//     // Debug information (optional)
//     // std::cout << nz_rm_row << " row elements removed" << std::endl;
//     // std::cout << nz_rm_col << " column elements removed" << std::endl;
//     // std::cout << A_map.size() << " elements remain in the matrix map" << std::endl;
// }

template <typename M>
void BoundaryDirichlet<M>::finalize_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                                                  std::span<double> b,
                                                  size_t dof_start)
{
    for (const auto& kv : boundary_dofs) {
        const int I_dof = static_cast<int>(kv.first + dof_start);
        const double g = b[static_cast<size_t>(I_dof)];

        for (auto it = A.begin(); it != A.end(); ) {
            const int r = it->first.first;
            const int c = it->first.second;
            if (c == I_dof && r != I_dof) {
                b[static_cast<size_t>(r)] -= it->second * g; // RHS shift
                it = A.erase(it);                            // zero A(r,I_dof)
            } else {
                ++it;
            }
        }
        // A(I_dof,I_dof)=1 already set; row I_dof already zeroed.
    }

}


#endif // BOUNDARYHPP