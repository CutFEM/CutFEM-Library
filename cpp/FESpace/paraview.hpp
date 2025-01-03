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
#ifndef PARAVIEW_HPP
#define PARAVIEW_HPP

#include "../common/cut_method.hpp"
#include "FESpace.hpp"
#include "macroElement.hpp"
#include "../num/util.hpp"
#include "expression.hpp"
#include "../algoim/quadrature_general.hpp"

// template<class F>
// class FEMFunction;
//
template <class F> class FunFEM;

template <class F> class FEMTimeFunction;

static double paraviewFormat(double x) {
    // return (fabs(x) < 1e-25)? util::fsign(x)*1e-25: x;
    return (fabs(x) < 1e-25) ? 0 : x;
}

/*
 *   New version of Parview
 *
 */

template <class M> class Paraview {
  public:
    typedef M Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename Mesh::Element Element;
    typedef typename Element::Rd Rd;
    typedef Partition<Element> PartitionT;
    typedef FunFEM<Mesh> Fun_h;
    typedef ExpressionVirtual Expression;

    const int nv_cut_element = Rd::d + 1;
    int nbDataFile           = 0;
    std::string outFile_;
    int ntCut = -1;
    int nt_cut, nt_notcut;

    struct ParaviewMesh {

        int ntCut_;
        int ntNotcut_;
        int nv_;

        int numCell_;
        int numCutCell_;

        int nvCutCell_;
        int nvCell_;

        std::map<int, int> element_status;
        std::vector<std::vector<Rd>> mesh_node;
        std::vector<std::pair<int, int>> idx_in_Vh; // idx_K -> (k_back, domain)
        std::vector<std::pair<int, int>> num_cell;  //(nb nodes, num_cell)

        ParaviewMesh()
            : numCell_(Element::ParaviewNumCell),
              numCutCell_((Rd::d == 2) ? Triangle2::ParaviewNumCell : Tet::ParaviewNumCell), nvCutCell_(Rd::d + 1),
              nvCell_(Element::nv) {}

        void check_and_resize_array(int next_idx) {
            int size0 = mesh_node.size();
            if (next_idx < size0)
                return;
            mesh_node.resize((int)(3 * size0 / 2));
            idx_in_Vh.resize((int)(3 * size0 / 2));
            num_cell.resize((int)(3 * size0 / 2));
        }

        void clearAndResize(int size0) {
            mesh_node.clear();
            idx_in_Vh.clear();
            num_cell.clear();
            mesh_node.resize(size0);
            idx_in_Vh.resize(size0);
            num_cell.resize(size0);
        }
        void shrinkToFit(int kk) {
            mesh_node.resize(kk);
            mesh_node.shrink_to_fit();
            idx_in_Vh.resize(kk);
            idx_in_Vh.shrink_to_fit();
            num_cell.resize(kk);
            num_cell.shrink_to_fit();
        }
        void build(const Mesh &Th) {

            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = Th.get_nb_element();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < Th.get_nb_element(); ++k) {

                idx_in_Vh[k] = std::make_pair(k, 0);
                num_cell[k]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[k].push_back(Th[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
            }
        }
        void build(const ActiveMesh<Mesh> &cutTh);
        void buildCut(const ActiveMesh<Mesh> &cutTh) {

            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            // std::vector<Rd> list_node;
            int kk    = 0;
            int size0 = 1.5 * cutTh.NbElement();
            clearAndResize(size0);
            double area = 0.;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);
                if (cutTh.isInactive(k, 0)) {
                    continue;
                }
                if (cutTh.isCut(k, 0)) {

                    const Cut_Part<Element> cutK(cutTh.get_cut_part(k, 0));

                    if (cutK.multi_interface()) {
                        assert(0);
                    }

                    for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

                        check_and_resize_array(kk);
                        for (int i = 0; i < nvCutCell_; ++i) {
                            Rd x(cutK.get_vertex(it, i));
                            mesh_node[kk].push_back(x);
                        }
                        num_cell[kk]  = std::make_pair(nvCutCell_, numCutCell_);
                        idx_in_Vh[kk] = std::make_pair(kb, domain);
                        nv_ += nvCutCell_;
                        kk++;
                        ntCut_++;
                    }

                } else {
                    // not cut
                    idx_in_Vh[kk] = std::make_pair(kb, domain);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[k][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildNoCut(const ActiveMesh<Mesh> &cutTh) {

            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = 1.5 * cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);
                // Sometimes the following three lines should be outcommented, to view the solution in all quadrature
                // points
                if (cutTh.isInactive(k, 0)) {
                    continue;
                }
                if (cutTh.isCut(k, 0)) {

                    const Cut_Part<Element> cutK(cutTh.get_cut_part(k, 0));

                    if (cutK.multi_interface()) {
                        assert(0);
                        //     // Here we need to create triangles to feel the multi cut
                        //     for (int l = 0; l < cutK.get_nb_element(); ++l) {
                        //         check_and_resize_array(kk);
                        //         CutElement<Element> K = cutK.get_element(l);
                        //         for (int i = 0; i < nvCutCell_; ++i) {
                        //             mesh_node[kk].push_back(K[i]);
                        //         }
                        //         num_cell[kk]  = std::make_pair(nvCutCell_, 7);
                        //         idx_in_Vh[kk] = std::make_pair(kb, domain);
                        //         nv_ += nvCutCell_;
                        //         ntCut_++;
                        //         kk++;
                        //     }
                    } else {
                        check_and_resize_array(kk);

                        cutK.get_list_node(list_node);
                        int nv_loc = list_node.size();

                        num_cell[kk]  = std::make_pair(nv_loc, 7);
                        idx_in_Vh[kk] = std::make_pair(kb, domain);
                        nv_ += nv_loc;
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(list_node[i]);
                        }
                        ntCut_++;
                        kk++;
                    }
                } else {
                    check_and_resize_array(kk);

                    idx_in_Vh[kk] = std::make_pair(kb, domain);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[k][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void build_active_mesh(const ActiveMesh<Mesh> &cutTh) {

            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
            shrinkToFit(kk + 1);
        }
        void build_face_stab(const ActiveMesh<Mesh> &cutTh, int domain) {

            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int kb = cutTh.idxElementInBackMesh(k);
                if ((cutTh.isStabilizeElement(k) && cutTh.get_domain_element(k) == domain) || domain == -1) {
                    const Element &K(cutTh[k]);
                    for (int e = 0; e < Element::nea; ++e) {
                        check_and_resize_array(kk);

                        int je = e;
                        int kn = cutTh.ElementAdj(k, je);

                        if (kn == -1)
                            continue;

                        int nv_loc    = Element::nva;
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        idx_in_Vh[kk] = std::make_pair(kb, domain);
                        nv_ += nv_loc;
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(K[Element::nvedge[e][i]]);
                        }
                        ntCut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }

        void build_macro_element(const MacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th_);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroInnerEdge(const MacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th_);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroOutterEdge(const MacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th_);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMeshNoStab(const MacroElement<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th_);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (macro.isSmall(k) || macro.isRootFat(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        void buildMeshNoStab(const ActiveMesh<Mesh> &cutTh) {
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                // if(domain != dom) continue;
                if (cutTh.isStabilizeElement(k))
                    continue;

                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }
        void buildMeshNoStabEdge(const ActiveMesh<Mesh> &cutTh) {
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement() * 3;
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                // if(domain != dom) continue;
                if (cutTh.isStabilizeElement(k))
                    continue;

                for (int e = 0; e < Mesh::Element::nea; ++e) {
                    int je = e;
                    int kn = cutTh.ElementAdj(k, je);
                    if (kn == -1)
                        continue;
                    if (cutTh.isStabilizeElement(kn))
                        continue;

                    check_and_resize_array(kk);

                    idx_in_Vh[kk] = std::make_pair(kb, domain);
                    num_cell[kk]  = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[k][Element::nvedge[e][i]]);
                    }
                    nv_ += nv_loc;

                    ntNotcut_++;
                    kk++;
                }
            }
        }

        void build_macro_element(const MacroElementSurface<M> &macro) {

            const Mesh &Th(*macro.interface.backMesh);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = Th.nbElements();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                // if(cutTh.get_domain_element(MK.get_index_root()) != dom)
                // continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int ki        = MK.get_index_element(k);
                    int kb        = macro.interface.idxElementOfFace(ki);
                    idx_in_Vh[kk] = std::make_pair(kb, 0);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(Th[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroOutterEdge(const MacroElementSurface<M> &macro) {

            const Mesh &Th(*macro.interface.backMesh);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = Th.nbElements();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                for (int k = 0; k < MK.size(); ++k) {
                    int ki = MK.get_index_element(k);
                    int kb = macro.interface.idxElementOfFace(ki);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = Th.ElementAdj(kb, je);

                        int kin;
                        if (macro.interface.isCut(kn)) {
                            kin = macro.interface.idxFaceOfElement(kn);
                        } else {
                            kin = -1;
                        }
                        if (MK.containElement(kin))
                            continue;

                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kb, 0);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(Th[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroInnerEdge(const MacroElementSurface<M> &macro) {

            const Mesh &Th(*macro.interface.backMesh);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = Th.nbElements();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int ki                   = edge.first;
                    int kb                   = macro.interface.idxElementOfFace(ki);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kb, 0);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(Th[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMeshNoStab(const MacroElementSurface<M> &macro) {
            const Interface<Mesh> &interface(macro.interface);
            const Mesh &Th(*macro.interface.backMesh);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = interface.nbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < interface.nbElement(); ++k) {

                if (macro.isSmall(k) || macro.isRootFat(k))
                    continue;
                // if(macro.isSmall(k)) continue;

                check_and_resize_array(kk);
                int kb = interface.idxElementOfFace(k);

                idx_in_Vh[kk] = std::make_pair(kb, 0);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(Th[kb][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        void build_macro_element(const TimeMacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void build_macro_element(const TimeMacroElement2<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void build_macro_element(const MacroElementPartition<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }

        void buildMacroInnerEdge(const TimeMacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroInnerEdge(const TimeMacroElement2<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroInnerEdge(const MacroElementPartition<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }

        void buildMacroOutterEdge(const TimeMacroElement<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroOutterEdge(const TimeMacroElement2<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroOutterEdge(const MacroElementPartition<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }

        void buildMeshNoStab(const TimeMacroElement<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (macro.isSmall(k) || macro.isRootFat(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }
        void buildMeshNoStab(const MacroElementPartition<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (macro.isSmall(k) || macro.isRootFat(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        void buildSmallElements(const TimeMacroElement<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }
        void buildSmallElements(const TimeMacroElement2<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }
        void buildSmallElements(const MacroElementPartition<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        void build_macro_element(const TimeMacroElementSurface<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroInnerEdge(const TimeMacroElementSurface<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildMacroOutterEdge(const TimeMacroElementSurface<M> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        void buildSmallElements(const TimeMacroElementSurface<M> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        void buildSmallElements(const MacroElementSurface<M> &macro) {
            assert(macro.Th_active);
            const ActiveMesh<Mesh> &cutTh(*(macro.Th_active));
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            int dom = 0;

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom) {
                    assert(0);
                    continue;
                }
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }

        template <typename L> void buildMacroInnerEdge(const AlgoimMacro<M, L> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.get_nb_inner_edge(); ++k) {
                    check_and_resize_array(kk);

                    std::pair<int, int> edge = MK.get_inner_edge(k);
                    int kb                   = edge.first;
                    int kbb                  = cutTh.idxElementInBackMesh(kb);
                    int ie                   = edge.second;
                    idx_in_Vh[kk]            = std::make_pair(kbb, dom);
                    num_cell[kk]             = std::make_pair(nv_loc, 3);
                    for (int i = 0; i < nv_loc; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                    }
                    nv_ += nv_loc;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }
        template <typename L> void buildMacroOutterEdge(const AlgoimMacro<M, L> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk     = 0;
            int nv_loc = M::Element::nva;
            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb  = MK.get_index_element(k);
                    int kbb = cutTh.idxElementInBackMesh(kb);
                    for (int ie = 0; ie < M::Element::nea; ++ie) {
                        int je = ie;
                        int kn = cutTh.ElementAdj(kb, je);
                        if (MK.containElement(kn))
                            continue;
                        check_and_resize_array(kk);

                        idx_in_Vh[kk] = std::make_pair(kbb, dom);
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(cutTh[kb][Element::nvedge[ie][i]]);
                        }
                        nv_ += nv_loc;
                        ntNotcut_++;
                        kk++;
                    }
                }
            }
            shrinkToFit(kk + 1);
        }
        template <typename L> void buildSmallElements(const AlgoimMacro<M, L> &macro, int dom) {
            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                if (domain != dom)
                    continue;
                if (!macro.isSmall(k))
                    continue;
                check_and_resize_array(kk);

                idx_in_Vh[kk] = std::make_pair(kb, domain);
                num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                for (int i = 0; i < nvCell_; ++i) {
                    mesh_node[kk].push_back(cutTh[k][i]);
                }
                nv_ += nvCell_;
                ntNotcut_++;
                kk++;
            }
        }
        template <typename L> void build_macro_element(const AlgoimMacro<M, L> &macro, int dom) {

            const ActiveMesh<Mesh> &cutTh(macro.Th);
            ntCut_    = 0;
            ntNotcut_ = 0;
            nv_       = 0;
            int size0 = cutTh.NbElement();
            clearAndResize(size0);
            std::vector<Rd> list_node;
            int kk = 0;

            for (auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
                check_and_resize_array(kk);

                const MElement &MK(it->second);
                if (cutTh.get_domain_element(MK.get_index_root()) != dom)
                    continue;

                for (int k = 0; k < MK.size(); ++k) {
                    int kb        = MK.get_index_element(k);
                    int kbb       = cutTh.idxElementInBackMesh(kb);
                    idx_in_Vh[kk] = std::make_pair(kbb, dom);
                    num_cell[kk]  = std::make_pair(nvCell_, numCell_);
                    for (int i = 0; i < nvCell_; ++i) {
                        mesh_node[kk].push_back(cutTh[kb][i]);
                    }
                    nv_ += nvCell_;
                    ntNotcut_++;
                    kk++;
                }
            }
            shrinkToFit(kk + 1);
        }

        template <typename L>
        void buildAlgoimQuadrature(const ActiveMesh<M> &cutTh, L &phi, const TimeSlab &In,
                                   const QuadratureFormular1d &qTime, const int itq, const int algoim_domain, const int side) {
            using mesh_t    = M;
            using fespace_t = GFESpace<mesh_t>;
            using FElement  = typename fespace_t::FElement;
            using Rd        = typename FElement::Rd;
            using Element   = typename mesh_t::Element;

            // algoim_domain = -1 -> integrate {phi < 0} \cap K
            // algoim_domain = 2 -> integrate {phi = 0} \cap K

            //assert(algoim_domain == -1 || algoim_domain == 2);
            if (algoim_domain == 0 || algoim_domain == 1)
                assert(side == 0 || side == 1); 

            const int quadrature_order = 5;
            const int domain           = 0; // only for one subdomain
            ntCut_                     = 0;
            ntNotcut_                  = 0;
            nv_                        = 0;

            clearAndResize(quadrature_order * quadrature_order * cutTh.NbElement());

            // std::vector<Rd> all_quadrature_nodes;

            GQuadraturePoint<R1> tq((qTime)[itq]);
            const double t = In.mapToPhysicalElement(tq);
            phi.t          = t;

            int kk = 0;
            for (int k = cutTh.first_element(); k < cutTh.last_element(); k += cutTh.next_element()) {
                if (domain != cutTh.get_domain_element(k))
                    continue;

                // if (cutTh.isInactive(k, itq) || !cutTh.isCut(k, itq))
                //     continue;
                if (cutTh.isInactive(k, itq))
                    continue;

                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);
                const Element &K(cutTh[k]);

                const auto &V0(K.at(0)); // vertex 0
                const auto &V2(K.at(2)); // vertex 2   diagonally opposed
                // R2 zero(0., 0.);
                // std::cout << "Paraview phi(0,0) = " << phi(zero) << "\n";
                // std::cout << "t = " << t << "\n";
                // std::cout << "kb = " << kb << "\n";
                // std::cout << "V0 = " << V0 << ", V2 = " << V2 << "\n";
                algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
                algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

                algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax),
                                                                 algoim_domain, side, quadrature_order);

                // if (algoim_domain == 2)
                //     assert((q.nodes.size() == quadrature_order) ||
                //     (q.nodes.size() == 2*quadrature_order) ||
                //     (q.nodes.size()) == 0);
                // else if (algoim_domain == -1)
                //     assert((q.nodes.size() == quadrature_order * quadrature_order) ||
                //            (q.nodes.size() == 2 * quadrature_order * quadrature_order)||
                //            (q.nodes.size() == 3 * quadrature_order * quadrature_order));

                // if (q.nodes.size() == 75) {
                //     std::cout << "k = " << k << ", itq = " << itq << ", q size: " << q.nodes.size() << "\n";
                //     std::cout << "kb = " << kb << "\n";
                //     std::cout << "V0 = " << V0 << ", V2 = " << V2 << "\n";
                // }
                for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
                    check_and_resize_array(kk);
                    mesh_node[kk].push_back(R2(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1)));
                    num_cell[kk]  = std::make_pair(1, 1);
                    idx_in_Vh[kk] = std::make_pair(kb, domain);
                    kk++;
                    ntNotcut_++;
                }

                nv_ += q.nodes.size();
            }
        }

        template <typename L> void buildAlgoimQuadrature(const ActiveMesh<M> &cutTh, L &phi, const int algoim_domain) {
            using mesh_t    = M;
            using fespace_t = GFESpace<mesh_t>;
            using FElement  = typename fespace_t::FElement;
            using Rd        = typename FElement::Rd;
            using Element   = typename mesh_t::Element;

            // algoim_domain = -1 -> integrate {phi < 0} \cap K
            // algoim_domain = 2 -> integrate {phi = 0} \cap K

            assert(algoim_domain == -1 || algoim_domain == 2);

            const int quadrature_order = 5;
            const int domain           = 0; // only for one subdomain
            ntCut_                     = 0;
            ntNotcut_                  = 0;
            nv_                        = 0;

            clearAndResize(quadrature_order * quadrature_order * cutTh.NbElement());

            // std::vector<Rd> all_quadrature_nodes;

            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); k += 1) {
                if (domain != cutTh.get_domain_element(k))
                    continue;

                // if (!cutTh.isCut(k, 0))
                //     continue;

                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);

                check_and_resize_array(kk);

                const Element &K(cutTh[k]);
                const auto &V0(K.at(0)); // vertex 0
                const auto &V2(K.at(2)); // vertex 2   diagonally opposed

                algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
                algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

                algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax),
                                                                 algoim_domain, -1, quadrature_order);

                if (algoim_domain == 2)
                    assert((q.nodes.size() == quadrature_order) || (q.nodes.size()) == 0);
                else if (algoim_domain == -1)
                    assert((q.nodes.size() == quadrature_order * quadrature_order) ||
                           (q.nodes.size() == 2 * quadrature_order * quadrature_order));

                for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
                    mesh_node.at(kk).push_back(R2(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1)));
                    num_cell.at(kk)  = std::make_pair(1, 1);
                    idx_in_Vh.at(kk) = std::make_pair(kb, domain);
                    kk++;
                    ntNotcut_++;
                }

                nv_ += q.nodes.size();
            }
        }

        template <typename Q> void buildQuadrature(const ActiveMesh<M> &cutTh, Q &qf) {
            using mesh_t    = M;
            using fespace_t = GFESpace<mesh_t>;
            using FElement  = typename fespace_t::FElement;
            using Rd        = typename FElement::Rd;
            using Element   = typename mesh_t::Element;

            const int n      = qf.getNbrOfQuads();
            const int domain = 0; // only for one subdomain
            ntCut_           = 0;
            ntNotcut_        = 0;
            nv_              = 0;

            clearAndResize(2 * n * cutTh.NbElement());

            // std::vector<Rd> all_quadrature_nodes;

            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); k += 1) {
                if (domain != cutTh.get_domain_element(k))
                    continue;

                if (cutTh.isCut(k, 0))
                    continue;
                const auto &K = cutTh[k];

                // if (cutTh.isCut(k, 0)) {

                // const Cut_Part<Element> cutK(cutTh.get_cut_part(k, 0));

                int domain = cutTh.get_domain_element(k);
                int kb     = cutTh.idxElementInBackMesh(k);
                check_and_resize_array(kk);

                // for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                    auto ip(qf[ipq]);

                    mesh_node.at(kk).push_back(K.mapToPhysicalElement(ip));
                    num_cell.at(kk)  = std::make_pair(1, 1);
                    idx_in_Vh.at(kk) = std::make_pair(kb, domain);
                    kk++;
                    ntNotcut_++;
                }

                //}

                nv_ += qf.getNbrOfQuads();
            }
        }

        bool isCut(int k) const {
            element_status.find(k);
            return (element_status.find(k) != element_status.end());
        }
        int numCell(int k) const {
            if (isCut(k))
                return numCutCell_;
            else
                return numCell_;
        }
        int nvCell(int k) const {
            if (isCut(k))
                return nvCutCell_;
            else
                return nvCell_;
        }
        int nbElement() const { return ntCut_ + ntNotcut_; }
        int nbNode() const { return nv_; }
        int num_stab_elems(const ActiveMesh<Mesh> &cutTh, int domain) {

            int nstab = 0; // number of stabilized edges

            int size0 = cutTh.NbElement();
            clearAndResize(size0);

            std::vector<Rd> list_node;
            int kk = 0;
            for (int k = 0; k < cutTh.NbElement(); ++k) {
                int kb = cutTh.idxElementInBackMesh(k);
                if ((cutTh.isStabilizeElement(k) && cutTh.get_domain_element(k) == domain) || domain == -1) {
                    const Element &K(cutTh[k]);
                    for (int e = 0; e < Element::nea; ++e) {
                        check_and_resize_array(kk);

                        int je = e;
                        int kn = cutTh.ElementAdj(k, je);

                        if (kn == -1)
                            continue;

                        int nv_loc    = Element::nva;
                        num_cell[kk]  = std::make_pair(nv_loc, 3);
                        idx_in_Vh[kk] = std::make_pair(kb, domain);

                        for (int i = 0; i < nv_loc; ++i) {
                            mesh_node[kk].push_back(K[Element::nvedge[e][i]]);
                        }
                        kk++;
                    }
                    nstab++;
                }
            }
            shrinkToFit(kk + 1);
            return nstab;
        }

        int sizeDataCell() const { return this->nbNode() + this->nbElement(); }
        Rd node(int k, int i) const { return mesh_node[k][i]; }

    } mesh_data;

    Paraview() {}
    Paraview(const ActiveMesh<Mesh> &cutTh, std::string name) {
        outFile_ = name;
        mesh_data.build(cutTh);
        this->writeFileMesh();
        this->writeFileCell();
    }
    Paraview(const Mesh &Th, std::string name) {
        outFile_ = name;
        mesh_data.build(Th);
        this->writeFileMesh();
        this->writeFileCell();
    }

    int get_nb_stab_elems(const ActiveMesh<Mesh> &cutTh, int domain) { return mesh_data.num_stab_elems(cutTh, domain); }
    void writeFaceStab(const ActiveMesh<Mesh> &cutTh, int domain, std::string name) {
        outFile_ = name;
        mesh_data.build_face_stab(cutTh, domain);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeActiveMesh(const ActiveMesh<Mesh> &cutTh, std::string name) {
        outFile_ = name;
        mesh_data.build_active_mesh(cutTh);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const MacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const MacroElementSurface<M> &macro, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const TimeMacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const TimeMacroElement2<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const MacroElementPartition<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroElement(const TimeMacroElementSurface<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const MacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const MacroElementSurface<M> &macro, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const TimeMacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const TimeMacroElement2<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const MacroElementPartition<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroInnerEdge(const TimeMacroElementSurface<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroOutterEdge(const MacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroOutterEdge(const MacroElementSurface<M> &macro, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroOutterEdge(const TimeMacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroOutterEdge(const TimeMacroElement2<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeMacroOutterEdge(const MacroElementPartition<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }

    void writeMacroOutterEdge(const TimeMacroElementSurface<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeNonStabMesh(const ActiveMesh<Mesh> &cutTh, std::string name) {
        outFile_ = name;
        mesh_data.buildMeshNoStab(cutTh);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeNonStabMeshEdge(const ActiveMesh<Mesh> &cutTh, std::string name) {
        outFile_ = name;
        mesh_data.buildMeshNoStabEdge(cutTh);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeNonStabMesh(const MacroElement<M> &macro, int domain, std::string name) {
        outFile_ = name;
        mesh_data.buildMeshNoStab(macro, domain);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeNonStabMesh(const MacroElementSurface<M> &macro, std::string name) {
        outFile_ = name;
        mesh_data.buildMeshNoStab(macro);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeNonStabMesh(const TimeMacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMeshNoStab(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeSmallElements(const TimeMacroElement<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeSmallElements(const TimeMacroElement2<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeSmallElements(const MacroElementPartition<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeSmallElements(const TimeMacroElementSurface<M> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    void writeSmallElements(const MacroElementSurface<M> &macro, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro);
        this->writeFileMesh();
        this->writeFileCell();
    }

    template <typename L> void writeMacroElement(const AlgoimMacro<M, L> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.build_macro_element(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    template <typename L> void writeMacroInnerEdge(const AlgoimMacro<M, L> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroInnerEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    template <typename L> void writeMacroOutterEdge(const AlgoimMacro<M, L> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildMacroOutterEdge(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }
    template <typename L> void writeSmallElements(const AlgoimMacro<M, L> &macro, int dom, std::string name) {
        outFile_ = name;
        mesh_data.buildSmallElements(macro, dom);
        this->writeFileMesh();
        this->writeFileCell();
    }

    template <typename L>
    void writeAlgoimQuadrature(const ActiveMesh<M> &cutTh, L &phi, const TimeSlab &In,
                               const QuadratureFormular1d &qTime, const int itq, const int algoim_domain, const int side,
                               std::string name) {
        outFile_ = name;
        mesh_data.template buildAlgoimQuadrature<L>(cutTh, phi, In, qTime, itq, algoim_domain, side);
        this->writeFileMesh();
        this->writeFileCell();
    }

    template <typename L>
    void writeAlgoimQuadrature(const ActiveMesh<M> &cutTh, L &phi, const int algoim_domain, std::string name) {
        outFile_ = name;
        mesh_data.template buildAlgoimQuadrature<L>(cutTh, phi, algoim_domain);
        this->writeFileMesh();
        this->writeFileCell();
    }

    template <typename Q> void writeAlgoimQuadrature(const ActiveMesh<M> &cutTh, Q &qf, std::string name) {
        outFile_ = name;
        mesh_data.template buildQuadrature<Q>(cutTh, qf);
        this->writeFileMesh();
        this->writeFileCell();
    }

    void writeFileMesh();
    void writeFileCell();

    void add(Fun_h &, std::string, int, int);
    void add(const ExpressionVirtual &fh, std::string);
    void add(const std::shared_ptr<ExpressionVirtual> &fh, std::string);

    void writeFileScalarData(const ExpressionVirtual &, std::string);
    void writeFileScalarData(const std::shared_ptr<ExpressionVirtual> &, std::string);
    void writeFileVectorData(Fun_h &, int, std::string);
};

// Writting the nodes
// // ---------------------------------------------------------------------------------------
// template <typename M, typename L> void Paraview<M>::writeAlgoimQuadrature(const ActiveMesh<M> &cutTh, L &phi, const
// double t, std::string name) {
//     // Prints the quadrature nodes of the interface integral, in the first time quadrature point of the given active
//     mesh.
//     // If the problem is stationary, just put t=0 (or any other number).
//     outFile_ = name;

//     using mesh_t        = M;
//     using fespace_t     = GFESpace<mesh_t>;
//     using FElement      = typename fespace_t::FElement;
//     using Rd            = typename FElement::Rd;
//     using Element       = typename mesh_t::Element;

//     std::ofstream point(outFile_.c_str(), std::ofstream::out);
//     point << "# vtk DataFile Version 1.0" << std::endl
//           << "unstructured Grid" << std::endl
//           << "ASCII" << std::endl
//           << "DATASET UNSTRUCTURED_GRID" << std::endl;

//     std::vector<R2> all_quadrature_nodes;

//     phi.t = t;
//     for (int k = cutTh.first_element(); k < cutTh.last_element(); k += cutTh.next_element()) {
//         if (cutTh.isCut(k, 0)) {
//             const Element &K(cutTh[k]);
//             const auto &V0(K.at(0)); // vertex 0
//             const auto &V2(K.at(2)); // vertex 2   diagonally opposed

//             algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
//             algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

//             algoim::QuadratureRule<2> q =
//                 algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, 5);
//             for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
//                 all_quadrature_nodes.push_back(R2(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1)));
//                 // Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
//                 // plot << mip << std::endl;
//             }
//         }
//     }

//     point << "POINTS " << all_quadrature_nodes.size() << " float " << std::endl;

//     // now print the coordinates

//     for (const auto & elem : all_quadrature_nodes) {
//         point << elem << "\t 0.0" << std::endl;
//     }

//     point.close();
// }

// Visualizing algoim quadrature
template <class M> void Paraview<M>::writeFileMesh() {

    std::ofstream point(outFile_.c_str(), std::ofstream::out);
    point << "# vtk DataFile Version 1.0" << std::endl
          << "unstructured Grid" << std::endl
          << "ASCII" << std::endl
          << "DATASET UNSTRUCTURED_GRID" << std::endl
          << "POINTS " << mesh_data.nbNode() << " float " << std::endl;

    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        for (int i = 0; i < mesh_data.mesh_node[k].size(); ++i) {
            if (Rd::d == 3) {
                point << mesh_data.node(k, i) << std::endl;
            } else {
                point << mesh_data.node(k, i) << "\t 0.0" << std::endl;
            }
        }
    }

    point.close();
}

// Writting the cells
// ---------------------------------------------------------------------------------------
template <class M> void Paraview<M>::writeFileCell() {

    std::ofstream cell(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
    cell << "CELLS " << mesh_data.nbElement() << " " << mesh_data.sizeDataCell() << std::endl;
    int l = 0;
    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        int nvc = mesh_data.num_cell[k].first;
        cell << nvc;
        for (int i = 0; i < nvc; ++i, ++l)
            cell << "  " << l;
        cell << std::endl;
    }

    cell << "CELL_TYPES " << mesh_data.nbElement() << std::endl;
    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        cell << mesh_data.num_cell[k].second << std::endl;
    }

    cell.close();
}

template <class M> void Paraview<M>::add(Fun_h &fh, std::string nameField, int begin_comp, int nb_comp) {

    if (nb_comp == 1) {
        ExpressionFunFEM<M> ui(fh, begin_comp, op_id);
        writeFileScalarData(ui, nameField);
    } else {
        writeFileVectorData(fh, begin_comp, nameField);
    }
}

template <class M> void Paraview<M>::add(const ExpressionVirtual &fh, std::string nameField) {

    writeFileScalarData(fh, nameField);
}

template <class M> void Paraview<M>::add(const std::shared_ptr<ExpressionVirtual> &fh, std::string nameField) {

    writeFileScalarData(fh, nameField);
}

template <class M> void Paraview<M>::writeFileScalarData(const ExpressionVirtual &fh, std::string name) {

    std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
    if (nbDataFile == 0)
        data << "POINT_DATA " << mesh_data.nbNode() << std::endl;
    data << "SCALARS " + name + " float" << std::endl;
    data << "LOOKUP_TABLE default" << std::endl;

    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        auto k_dom = mesh_data.idx_in_Vh[k];
        int kb     = k_dom.first;
        int domain = k_dom.second;
        for (int i = 0; i < mesh_data.mesh_node[k].size(); ++i) {
            R1 val = fh.evalOnBackMesh(kb, domain, mesh_data.node(k, i));
            data << paraviewFormat(val) << std::endl;
            ;
        }
    }
    data.close();
    nbDataFile++;
}

template <class M>
void Paraview<M>::writeFileScalarData(const std::shared_ptr<ExpressionVirtual> &fh, std::string name) {

    std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
    if (nbDataFile == 0)
        data << "POINT_DATA " << mesh_data.nbNode() << std::endl;
    data << "SCALARS " + name + " float" << std::endl;
    data << "LOOKUP_TABLE default" << std::endl;

    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        auto k_dom = mesh_data.idx_in_Vh[k];
        int kb     = k_dom.first;
        int domain = k_dom.second;
        for (int i = 0; i < mesh_data.mesh_node[k].size(); ++i) {
            R1 val = fh->evalOnBackMesh(kb, domain, mesh_data.node(k, i));
            data << paraviewFormat(val) << std::endl;
            ;
        }
    }
    data.close();
    nbDataFile++;
}

template <class M> void Paraview<M>::writeFileVectorData(Fun_h &fh, int c0, std::string name) {

    std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);

    if (nbDataFile == 0)
        data << "POINT_DATA " << mesh_data.nbNode() << std::endl;
    data << "VECTORS " + name + " float" << std::endl;

    for (int k = 0; k < mesh_data.nbElement(); ++k) {
        auto k_dom = mesh_data.idx_in_Vh[k];
        int kb     = k_dom.first;
        int domain = k_dom.second;
        int kf     = fh.idxElementFromBackMesh(kb, domain);

        for (int i = 0; i < mesh_data.mesh_node[k].size(); ++i) {

            for (int dd = 0; dd < Rd::d; ++dd) {
                // R1 val = fh.evalOnBackMesh(kb, domain, mesh_data.node(k,i),
                // c0+dd);
                R1 val = fh.eval(kf, mesh_data.node(k, i), c0 + dd);
                data << paraviewFormat(val) << "\t";
            }
            if (Rd::d == 2)
                data << "0.0";
            data << std::endl;
        }
    }
    data.close();
    nbDataFile++;
}

//
// template<class M>
// void Paraview<M>::add(Fun_h& fh, std::string nameField, int begin_comp, int
// nb_comp, MacroElement& macro){
//
//   if(nb_comp ==1) {
//     ExpressionFunFEM<M> ui(fh, begin_comp, op_id);
//     writeFileScalarData(ui, nameField, macro);
//   }
//   else
//     writeFileVectorData(fh, begin_comp, nameField, macro);
// }
//
//
// template<class M>
// void Paraview<M>::add(const ExpressionVirtual& fh, std::string nameField,
// MacroElement& macro){
//
//     writeFileScalarData(fh, nameField, macro);
// }
//
//
// template<class M>
// void Paraview<M>::writeFileVectorData(Fun_h& fh,int c0, std::string name,
// MacroElement& macro){
//
//   std::ofstream data(outFile_.c_str(), std::ofstream::out |
//   std::ofstream::app);
//
//   if(nbDataFile == 0) data << "POINT_DATA " << ntCut*Element::nv
// 			   << std::endl;
//   data << "VECTORS "+ name +" float" << std::endl;
//
//   if(levelSet) {
//     double loc_ls[Element::nv];
//     for(int k=0; k<Vh.NbElement(); ++k) {
//
//       const FElement& FK(Vh[k]);
//       int kback = Vh.idxElementInBackMesh(k);
//
//       levelSet->eval(loc_ls, Vh.Th(FK.T));
//       const int domain = FK.whichDomain();
//       const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
//       ElementSignEnum the_part = cutK.what_part(domain);
//
//
//       int kf = fh.idxElementFromBackMesh(kback, domain);
//       if(!macro.isRootFat(k)) {
//         kf = macro.getIndexRootElement(k);
//       }
//
//       if(the_part == NoElement) continue;
//
//       if(cutK.is_cut()) {
//         for(typename PartitionT::const_element_iterator it =
//         cutK.element_begin(the_part); it != cutK.element_end(the_part);
//         ++it){
//           for(int i=0;i<Element::nv;++i) {
//             Rd x(cutK.get_Vertex(it,i));
//             for(int  dd=0;dd<Rd::d;++dd){
//               R val = fh.eval(kf, x, c0+dd);
//               data << paraviewFormat(val) << "\t" ;
//             }
//             if (Rd::d==2) data << "0.0";
//             data << std::endl;
//           }
//         }
//       }
//       else
//       for(int i=0;i<Element::nv;++i) {
//         for(int  dd=0;dd<Rd::d;++dd){
//           R val = fh.eval(kf, FK.T[i], c0+dd);
//           data << paraviewFormat(val) << "\t" ;
//         }
//         if (Rd::d==2) data << "0.0";
//         data << std::endl;
//       }
//     }
//   }
//   else{
//
//     for(int k=0;k<Vh.NbElement();++k) {
//       const FElement& FK(Vh[k]);
//
//       int kback = Vh.idxElementInBackMesh(k);
//       int kf = fh.idxElementFromBackMesh(kback);
//
//       for(int i=0;i<Element::nv;++i){
//         for(int  dd=0;dd<Rd::d;++dd){
//           R val = fh.eval(kf, FK.T[i], c0+dd);
//           data << paraviewFormat(val) << "\t" ;
//         }
//         if (Rd::d==2) data << "0.0";
//         data << std::endl;
//       }
//     }
//   }
//   data.close();
//   nbDataFile++;
// }
//
//
// template<class M>
// void Paraview<M>::writeFileScalarData(const ExpressionVirtual& fh,
// std::string name, MacroElement& macro){
//
//   std::ofstream data(outFile_.c_str(), std::ofstream::out |
//   std::ofstream::app); if(nbDataFile == 0) data << "POINT_DATA " <<
//   ntCut*Element::nv << std::endl; data << "SCALARS "+ name +" float" <<
//   std::endl; data << "LOOKUP_TABLE default" << std::endl;
//
//   double loc_ls[Element::nv];
//
//   if(levelSet){
//     for(int k=0; k<Vh.NbElement(); ++k) {
//       const FElement& FK(Vh[k]);
//       levelSet->eval(loc_ls, Vh.Th(FK.T));
//
//       const int domain = FK.whichDomain();
//       const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
//       ElementSignEnum the_part = cutK.what_part(domain);
//       int kback = Vh.idxElementInBackMesh(k);
//
//       if(!macro.isRootFat(k)) {
//         int kk = macro.getIndexRootElement(k);
//         kback = Vh.idxElementInBackMesh(kk);
//       }
//
//       if(the_part == NoElement) continue;
//
//       if(cutK.is_cut()) {
//         for(typename PartitionT::const_element_iterator it =
//         cutK.element_begin(the_part); it != cutK.element_end(the_part);
//         ++it){
//           for(int i=0;i<Element::nv;++i) {
//             Rd x(cutK.get_Vertex(it,i));
//             R1 val = fh.evalOnBackMesh(kback, domain, x);
//             // if(!macro.isRootFat(k)) val = (2*domain-1)*1e7;
//             data << paraviewFormat(val) << std::endl;;
//           }
//         }
//       }
//       else{
//         for(int i=0;i<Element::nv;++i) {
//           R1 val = fh.evalOnBackMesh(kback, domain, FK.T[i]);
//           data << paraviewFormat(val.x) << std::endl;
//         }
//       }
//     }
//   }
//   else {
//     for(int k=0;k<Vh.NbElement();++k) {
//       const FElement& FK(Vh[k]);
//       int kback = Vh.idxElementInBackMesh(k);
//       for(int i=0;i<Element::nv;++i){
//         R1 val = fh.evalOnBackMesh(kback, -1, FK.T[i]);
//         data << paraviewFormat(val.x) << std::endl;
//       }
//     }
//   }
//
//   data.close();
//   nbDataFile++;
// }
//

#endif
