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
#ifndef _MACRO_ELEMENT_HPP
#define _MACRO_ELEMENT_HPP

#include "../common/base_interface.hpp"
#include "../common/cut_mesh.hpp"
#include "../common/global.hpp"
#include "../algoim/quadrature_general.hpp"
#include "../FESpace/FESpace.hpp"

#include <set>

class Extension;
class GMacro;

struct MElement {
  public:
    int idx_root_element;
    std::vector<int> idx_element;
    std::vector<std::pair<int, int>> inner_edge;  // idx_Element, idx_edge
    std::vector<std::pair<int, int>> outter_edge; // idx_Element, idx_edge

    GMacro *macro_ptr = nullptr;

    double area_root_  = 0;
    double area_total_ = 0;

    MElement(int idx_root = -1, GMacro *macroo_ptr = nullptr) : macro_ptr(macroo_ptr) {
        idx_root_element = idx_root;
        idx_element.push_back(idx_root);
    }
    MElement(int idx_root, double s, GMacro *macroo_ptr) : macro_ptr(macroo_ptr) {
        idx_root_element = idx_root;
        idx_element.push_back(idx_root);
        area_root_ = s;
        area_total_ += s;
    }
    int get_index_root() const { return idx_root_element; }
    int size() const { return idx_element.size(); }
    void add(int idx_K, std::pair<int, int> ke, double s) {
        idx_element.push_back(idx_K);
        // int kk =(idx_root_element < idx_K)? idx_root_element : idx_K;
        inner_edge.push_back(ke);
        area_total_ += s;
    }
    void add(int idx_K, double s) {
        idx_element.push_back(idx_K);
        area_total_ += s;
    }
    int get_index_element(int k) const { return idx_element[k]; }
    int get_nb_inner_edge() const { return inner_edge.size(); }
    std::pair<int, int> get_inner_edge(int k) const { return inner_edge[k]; }
    bool containElement(int k) const {
        for (int i = 0; i < idx_element.size(); ++i) {
            if (idx_element[i] == k)
                return true;
        }
        return false;
    }
    void print() const {
        std::cout << " Root idx   : \t " << idx_root_element << std::endl;
        std::cout << " small element : \n";
        for (auto it = inner_edge.begin(); it != inner_edge.end(); ++it) {
            std::cout << it->first << " with edge " << it->second << std::endl;
        }
    }

    const GMacro &getMacroSpace() const {
        if (!macro_ptr) {
            std::runtime_error("Macro pointer is not set");
        }
        return *macro_ptr;
    }
};

struct SmallElement {
    int index;
    int index_root;
    int chain_position;
    int idx_edge_to_root;
    double area;

    SmallElement(int idx = -1, int idx_root = -1)
        : index(idx), index_root(idx_root), chain_position(0), idx_edge_to_root(-1) {}
    void setRoot(int i) { index_root = i; }
    void setChainPosition(int i) { chain_position = i; }
    void setEdgeDirection(int i) { idx_edge_to_root = i; }
};

class GMacro {
  public:
    const int small = -1;
    const int good = 0, extension = 1, exhaust = 2;

    std::map<int, MElement> macro_element;                  // idx_root -> idx_macroElement
    std::map<int, SmallElement> small_element;              // idx_element -> idx_small_element
    std::map<int, MElement *> idx_element_to_macro_element; // idx_element -> macroElement_ptr

    double tol;

    GMacro() {}

    const auto begin() const { return macro_element.begin(); }
    const auto end() const { return macro_element.end(); }
    auto begin() { return macro_element.begin(); }
    auto end() { return macro_element.end(); }

    int getIndexRootElement(int k) const {
        auto it = small_element.find(k);
        if (it == small_element.end()) {
            return k;
        } else {
            return it->second.index_root;
        }
    }

    const SmallElement &getSmallElement(int k) const {
        auto it = small_element.find(k);
        assert(it != small_element.end());
        return it->second;
    }

    int nb_macro_element() const { return macro_element.size(); }

    bool isRootFat(int k) const { return (macro_element.find(k) != macro_element.end()); }
    bool isSmall(int k) const { return (small_element.find(k) != small_element.end()); }
    bool isInMacro(int k) const { return (idx_element_to_macro_element.find(k) != idx_element_to_macro_element.end()); }

    const MElement &getMacroElementOfElement(int k) const {
        auto it = idx_element_to_macro_element.find(k);
        assert(it != idx_element_to_macro_element.end());
        return *it->second;
    }
    //  int getMacroElementOfElement(int k) const {
    //    if(!isInMacro(k)) return -1;
    //    // if(isSmall(k)) return getSmallElement(k).index_root;
    //  }
};

template <typename Mesh> class MacroElement : public GMacro {
  public:
    const ActiveMesh<Mesh> &Th_;
    R tol_;
    int nb_element_0, nb_element_1;

    MacroElement(const ActiveMesh<Mesh> &th, const double C);

    double get_area(int k) const {
        if (isSmall(k)) {
            return getSmallElement(k).area;
        } else if (isRootFat(k)) {
            const auto it(macro_element.find(k));
            return it->second.area_root_;
        } else
            assert(0);
    }

  private:
    void findSmallElement();
    void createMacroElement();
    void createMacroElementInside();
    void findPathToInside(int, std::vector<std::pair<int, int>> &);
    void setOutterEdgeMacroElement();
    void setListElementToMacroElement();
    friend class Extension;
};

template <typename Mesh> MacroElement<Mesh>::MacroElement(const ActiveMesh<Mesh> &th, const double C) : Th_(th) {
    double h     = Th_[0].lenEdge(0);
    double meas  = Th_[0].measure();
    nb_element_0 = 0;
    nb_element_1 = 0;
    tol          = C * meas;

    findSmallElement();
    if (globalVariable::verbose > 0) {
        LOG_INFO << " ---  INFO MACRO ELEMENT  --- " << logger::endl;
        LOG_INFO << " Tolerance small element  :\t" << tol << logger::endl;
        LOG_INFO << " Small element in Omega 1 :\t" << nb_element_0 << logger::endl;
        LOG_INFO << " Small element in Omega 2 :\t" << nb_element_1 << logger::endl;
    }
    createMacroElement();
    setOutterEdgeMacroElement();
    setListElementToMacroElement();
    // createMacroElementInside();
}

template <typename Mesh> void MacroElement<Mesh>::findSmallElement() {
    for (int k = 0; k < Th_.get_nb_element(); k += 1) {
        if (!Th_.isCut(k, 0))
            continue;

        const typename Mesh::Element &K(Th_[k]);

        const Cut_Part<typename Mesh::Element> cutK(Th_.get_cut_part(k, 0));
        const int domain = Th_.get_domain_element(k);
        double areaCut   = cutK.measure();

        if (areaCut < tol) {
            small_element[k]      = SmallElement(k);
            small_element[k].area = areaCut;
            if (domain == 0) {
                nb_element_0++;
            } else {
                nb_element_1++;
            }
        }
    }
}

template <typename Mesh> void MacroElement<Mesh>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th_.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);
        ;
        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // lLOOP OVER FACES
            for (int ifac = 0; ifac < 3; ++ifac) {

                int ifacn = ifac;
                int kn    = Th_.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if (small_or_fat_K[kn] == small)
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;
                if (it != macro_element.end()) { // already exist
                    it->second.add(k, std::make_pair(kk, ie), Ks.area);
                } else {

                    const Cut_Part<typename Mesh::Element> cutK(Th_.get_cut_part(root_id, 0));
                    double areaCut = cutK.measure();

                    macro_element[root_id] = MElement(root_id, areaCut, this);
                    macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }
}

template <typename Mesh> void MacroElement<Mesh>::setOutterEdgeMacroElement() {
    // loop over macro elements
    for (auto &[idx_root, MK] : macro_element) {
        // loop over element contained in the macro element
        for (int k = 0; k < MK.size(); ++k) {
            int ki = MK.get_index_element(k);
            // loop over each edges of the element
            for (int ie = 0; ie < Mesh::Element::nea; ++ie) {
                auto [kn, je] = Th_.elementAdjacent(ki, ie);
                // if no neighbor or inner edge we do not add the edge as outter edge
                if (kn == -1 || MK.containElement(kn)) {
                    continue;
                }

                // else it is added to the list of outter edge
                MK.outter_edge.push_back(std::make_pair(ki, ie));
            }
        }
    }
}

template <typename Mesh> void MacroElement<Mesh>::setListElementToMacroElement() {
    // loop over macro elements
    for (auto &[idx_root, MK] : macro_element) {
        for (int k = 0; k < MK.size(); ++k) {
            int ki                           = MK.get_index_element(k);
            idx_element_to_macro_element[ki] = &MK;
        }
    }
}

template <typename Mesh> class MacroElementSurface : public GMacro {

  public:
    const Interface<Mesh> &interface;
    const ActiveMesh<Mesh> *Th_active;

    // MacroElementSurface(const Interface<Mesh> &, const double);
    MacroElementSurface(const Interface<Mesh> &, const double, const ActiveMesh<Mesh> *Th_active_ = nullptr);
    void findSmallElement();
    virtual void findRootElement();
    int checkDirection(const int, const int, int &);
};

template <typename Mesh>
MacroElementSurface<Mesh>::MacroElementSurface(const Interface<Mesh> &gh, const double C, const ActiveMesh<Mesh> *Th_active_) : interface(gh), Th_active(Th_active_) {
    double h = (*interface.backMesh)[0].lenEdge(0);
    tol      = C * h;

    // std::cout << " tolerance macro surface\t" << tol << std::endl;
    // findSmallElement();
    // std::cout << " Found " << small_element.size() << " small elements " << std::endl;
    // findRootElement();
}

template <typename Mesh> void MacroElementSurface<Mesh>::findSmallElement() {
    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const typename Interface<Mesh>::Face &face = interface[iface];
        const int kb                               = interface.idxElementOfFace(iface);
        const R meas                               = interface.measure(face);

        if (meas > tol)
            continue;

        small_element[iface]      = SmallElement(iface);
        small_element[iface].area = meas;
    }
}

template <typename Mesh> void MacroElementSurface<Mesh>::findRootElement() {
    const Mesh &Th(*interface.backMesh);
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {

        int iface                                  = it->first;
        const typename Interface<Mesh>::Face &face = interface[iface];
        const int kb                               = interface.idxElementOfFace(iface);
        KN<int> root_v(2), chain_position_v(2), ie_v(2);

        for (int i = 0; i < 2; ++i) {
            int chain_position = 0;

            int idx_node = face[i];
            int ie       = interface.edge_of_node_[idx_node];
            int root_K   = checkDirection(kb, ie, chain_position);

            root_v(i)           = root_K;
            chain_position_v(i) = chain_position;
            ie_v(i)             = ie;
        }

        // Choose the root with the lowest mesh index
        int i = (chain_position_v(0) <= chain_position_v(1)) ? 0 : 1;
        if (root_v(i) == -1) {
            i = (i == 0);
        }
        int root = interface.face_of_element_.find(root_v(i))->second;
        it->second.setRoot(root);
        it->second.setChainPosition(chain_position_v(i));
        it->second.setEdgeDirection(ie_v(i));
    }

    for (auto it = small_element.begin(); it != small_element.end(); ++it) {

        int idx_root = it->second.index_root;
        auto pRoot   = macro_element.find(idx_root);
        int k_loc    = it->first;
        int k        = interface.idxElementOfFace(k_loc);
        int ie       = it->second.idx_edge_to_root;
        int je       = ie;
        int kn       = Th.ElementAdj(k, je);
        assert(kn != -1);
        int kn_loc = interface.idxFaceOfElement(kn);
        assert(kn != -1);
        int kk = (k < kn) ? k_loc : kn_loc;
        int ke = (k < kn) ? ie : je;
        if (pRoot != macro_element.end()) { // already exist
            pRoot->second.add(k_loc, std::make_pair(kk, ke), 0);
        } else {
            macro_element[idx_root] = MElement(idx_root, this);
            macro_element[idx_root].add(k_loc, std::make_pair(kk, ke), 0.);
        }
    }
}

template <typename Mesh> int MacroElementSurface<Mesh>::checkDirection(const int k, const int ie, int &chain_position) {
    assert(chain_position < interface.nbElement());
    const Mesh &Th(*interface.backMesh);
    int e1 = ie;
    int kn = Th.ElementAdj(k, e1);

    if (kn == -1) {
        return -1;
    }

    chain_position++;
    int kn_loc = interface.idxFaceOfElement(kn);
    if (small_element.find(kn_loc) == small_element.end()) { // fat element
        return kn;
    }
    // find which edge to look at
    int iface                                  = interface.face_of_element_.find(kn)->second;
    const typename Interface<Mesh>::Face &face = interface[iface];
    int ie_next =
        (interface.edge_of_node_[face[0]] == e1) ? interface.edge_of_node_[face[1]] : interface.edge_of_node_[face[0]];

    return checkDirection(kn, ie_next, chain_position);
}



template <typename Mesh> class MacroElementSurfaceExtension : public MacroElementSurface<Mesh> {

  public:

    // MacroElementSurface(const Interface<Mesh> &, const double);
    MacroElementSurfaceExtension(const Interface<Mesh> &, const double, const ActiveMesh<Mesh> *Th_active_ = nullptr);
    void findRootElement() override; 
};

template <typename Mesh>
MacroElementSurfaceExtension<Mesh>::MacroElementSurfaceExtension(const Interface<Mesh> &gh, const double C, const ActiveMesh<Mesh> *Th_active_) : MacroElementSurface<Mesh>(gh, C, Th_active_) {}

template <typename Mesh> void MacroElementSurfaceExtension<Mesh>::findRootElement() {
    const Mesh &Th(*(this->interface).backMesh);
    
    // Create a map to keep track of the number of small elements connected to each large element
    std::unordered_map<int, int> largeElementCount;

    for (auto it = this->small_element.begin(); it != this->small_element.end(); ++it) {
        int iface                                  = it->first;
        const typename Interface<Mesh>::Face &face = this->interface[iface];
        const int kb                               = this->interface.idxElementOfFace(iface);
        KN<int> root_v(2), chain_position_v(2), ie_v(2);

        // Find potential roots
        for (int i = 0; i < 2; ++i) {
            int chain_position = 0;

            int idx_node = face[i];
            int ie       = this->interface.edge_of_node_[idx_node];
            int root_K   = this->checkDirection(kb, ie, chain_position);    // index of root (-1) if no root found in that direction

            if ()
            root_v(i)           = root_K;
            chain_position_v(i) = chain_position;
            ie_v(i)             = ie;
        }


        for (int i=0; i<2; ++i) {
            std::cout << "root_v(i) = " << root_v(i) << "\t chain_pos(i)" << chain_position_v(i) << "\t ie_v(i)" << ie_v(i) << std::endl;
        }

        // getchar();

        // Find the best root considering the count of connected small elements
        int selected_root = -1;
        int local_face_idx_selected_root = -1;
        int min_count = std::numeric_limits<int>::max();    // Initialize min_count to maximum value
        
        for (int i = 0; i < 2; ++i) {
            if (root_v(i) != -1) {
                int count = largeElementCount[root_v(i)];

                if (count < min_count) {
                    min_count = count;
                    selected_root = root_v(i);
                    local_face_idx_selected_root = i;
                }
            }
            else {
                selected_root = -1;
                std::cout << "No root found for small element " << iface << std::endl;
                break;
            }
        }
        
        //! Error is in the function below
        if (selected_root != -1) {
            assert(local_face_idx_selected_root != -1);
            int root = this->interface.face_of_element_.find(selected_root)->second;
            largeElementCount[selected_root]++;
            it->second.setRoot(root);
            it->second.setChainPosition(chain_position_v(local_face_idx_selected_root));
            it->second.setEdgeDirection(ie_v(local_face_idx_selected_root));
        }
        
    }

    for (auto it = this->small_element.begin(); it != this->small_element.end(); ++it) {
        int idx_root = it->second.index_root;
        auto pRoot   = this->macro_element.find(idx_root);
        int k_loc    = it->first;
        int k        = this->interface.idxElementOfFace(k_loc);
        int ie       = it->second.idx_edge_to_root;
        int je       = ie;
        int kn       = Th.ElementAdj(k, je);
        assert(kn != -1);
        int kn_loc = this->interface.idxFaceOfElement(kn);
        assert(kn != -1);
        int kk = (k < kn) ? k_loc : kn_loc;
        int ke = (k < kn) ? ie : je;
        if (pRoot != this->macro_element.end()) { // already exist
            pRoot->second.add(k_loc, std::make_pair(kk, ke), 0);
        } else {
            this->macro_element[idx_root] = MElement(idx_root, this);
            this->macro_element[idx_root].add(k_loc, std::make_pair(kk, ke), 0.);
        }
    }


    // print out the largeElementCount
    for (auto it = largeElementCount.begin(); it != largeElementCount.end(); ++it) {
        std::cout << "Large element " << it->first << " has " << it->second << " connected small elements" << std::endl;
    }
}


template <typename Mesh> class TimeMacroElementSurface : public GMacro {

  public:
    const ActiveMesh<Mesh> &Th;
    const TimeInterface<Mesh> &interface;
    R tol;

    TimeMacroElementSurface(const ActiveMesh<Mesh> &, const TimeInterface<Mesh> &, const QuadratureFormular1d &,
                            const double);
    void findSmallElement();
    void createMacroElement();
    int checkDirection(const int, const int, int &);

    int number_of_inner_edges();

  private:
    const QuadratureFormular1d &qTime;
};

template <typename Mesh>
TimeMacroElementSurface<Mesh>::TimeMacroElementSurface(const ActiveMesh<Mesh> &Th_, const TimeInterface<Mesh> &gh,
                                                       const QuadratureFormular1d &qTime_, const double C)
    : Th(Th_), interface(gh), qTime(qTime_) {
    double h = (*(*interface(0)).backMesh)[0].lenEdge(0);

    tol = C * h;

    std::cout << " tolerance macro surface\t" << tol << std::endl;
    findSmallElement();
    std::cout << " Found " << small_element.size() << " small elements " << std::endl;
    createMacroElement();
}

template <typename Mesh> void TimeMacroElementSurface<Mesh>::findSmallElement() {

    for (int k = 0; k < Th.get_nb_element(); k += 1) {

        // if (!Th.isStabilizeElement(k))
        //    continue;

        bool is_large           = false; // is element large in any quadrature point?
        bool is_inactive        = false; // is element inactive in any quadrature point?
        bool is_small           = false;
        int numb_times_inactive = 0;

        for (int itq = 0; itq < qTime.n; ++itq) {

            if (Th.isCut(k, itq)) {

                int kb    = Th.idxElementInBackMesh(k);
                int iface = interface(itq)->idxFaceOfElement(kb);

                const typename Interface<Mesh>::Face &face = (*interface(itq))[iface];
                const R meas                               = (*interface(itq)).measure(face);
                // std::cout << "k = " << k << "\n";
                // std::cout << "meas = " << meas << "\n";
                // std::cout << "tol = " << tol << "\n";

                if (meas > tol) {
                    is_large = true;
                } else
                    is_small = true;

            } else if (Th.isInactive(k, itq))
                is_inactive = true;
            // else if (Th.isInactive(k,itq)) ++numb_times_inactive;
        }

        if (!is_large || is_inactive) {
            // if (is_small || is_inactive) {
            // if (!is_large) {
            //  if (!is_large || numb_times_inactive>=2) {

            small_element[k] = SmallElement(k);
        }
    }
}

template <typename Mesh> void TimeMacroElementSurface<Mesh>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);
        ;
        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // lLOOP OVER FACES
            for (int ifac = 0; ifac < 3; ++ifac) {

                int ifacn = ifac;
                int kn    = Th.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if (small_or_fat_K[kn] == small)
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;

                if (it != macro_element.end()) { // already exist
                    it->second.add(k, std::make_pair(kk, ie), Ks.area);
                } else {

                    const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
                    double areaCut = cutK.measure();

                    macro_element[root_id] = MElement(root_id, areaCut, this);
                    macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }
}

template <typename Mesh> int TimeMacroElementSurface<Mesh>::number_of_inner_edges() {
    int num_of_inner_edges = 0;
    for (auto me = this->macro_element.begin(); me != this->macro_element.end(); ++me) {
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
            num_of_inner_edges += 1;
        }
    }
    return num_of_inner_edges;
}

/**
 * @brief General class for macro element partition of bulk domain.
 *
 * @tparam Mesh
 * @param Th Active mesh
 * @param tol Tolerance determining if an element is small or not
 * @note This class is used to generate a macro element partition of the active mesh
 * in both time dependent and stationary problems. Since the quadrature rule in time
 * is given by the active mesh object, the input arguments only need be the active mesh
 * and the tolerance parameter in both cases.
 */
template <typename Mesh> class MacroElementPartition : public GMacro {

  public:
    const ActiveMesh<Mesh> &Th;
    double tol, C;

    int nb_element_0,
        nb_element_1; // number of small elements in outer and inner domain respectively w.r.t level-set function
                      // sign

    const int number_of_faces      = Th[0].ne; // number of faces of the mesh element
    int number_of_stabilized_edges = 0;

    MacroElementPartition(const ActiveMesh<Mesh> &, const double);

    double get_area(int k) const;

  private:
    void findSmallElement();
    void createMacroElement();
    void setInnerEdges();
};

template <typename Mesh> double MacroElementPartition<Mesh>::get_area(int k) const {
    if (isSmall(k)) {
        return getSmallElement(k).area;
    } else if (isRootFat(k)) {
        const auto it(macro_element.find(k));
        return it->second.area_root_;
    } else
        assert(0);
}

template <typename Mesh>
MacroElementPartition<Mesh>::MacroElementPartition(const ActiveMesh<Mesh> &Th_, const double C_) : Th(Th_), C(C_) {

    //double h       = Th[0].lenEdge(1); // catheter of triangle
    double measure = Th[0].measure();  // measure = h^2/2

    nb_element_0 = 0;
    nb_element_1 = 0;
    tol          = C * measure;

    std::cout << "Tolerance ratio: \t" << C << std::endl;
    findSmallElement();
    std::cout << nb_element_0 << " \t small elements in Omega 1 " << std::endl;
    std::cout << nb_element_1 << " \t small elements in Omega 2 " << std::endl;
    createMacroElement();
    std::cout << "Macro element created\n";
    setInnerEdges();
    std::cout << "Inner edges set\n";
}

template <typename Mesh> void MacroElementPartition<Mesh>::findSmallElement() {

    // Iterate over all elements in the active mesh (over the whole time-slab)

    for (int k = 0; k < Th.get_nb_element(); k += 1) {

        if (!Th.isStabilizeElement(k))
            continue; // if the element is not cut or if it doesn't change domain it doesn't need stabilization

        const typename Mesh::Element &K(Th[k]);

        const int domain = Th.get_domain_element(k);

        // Iterate over the quadrature points in the time-slab In

        bool is_large           = false; // is element large in any quadrature point?
        bool is_small           = false; // is element large in any quadrature point?
        bool is_inactive        = false; // is element inactive in any quadrature point?
        int numb_times_inactive = 0;

        for (int itq = 0; itq < Th.nb_quadrature_time_; ++itq) {

            Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(k, itq));
            double areaCut = cutK.measure();
            double part    = areaCut / K.measure();

            if ((areaCut > tol) && (!Th.isInactive(k, itq))) {

                // element can be cut or not cut, but must be active in itq
                assert(0 < part && part < 1+1e-10); 

                is_large = true;
                
            } else if ((areaCut <= tol) && (!Th.isInactive(k, itq))) {

                // element must be cut and active in itq
                assert(0 <= part && part < C); // make sure cut area is not equal to element area if cut

                is_small = true;
            }

            else {
                // element lies entirely outside the domain in itq
                
                assert((Th.isInactive(k, itq)));
                assert(part == 1);  // if an element is inactive in itq, its size will be set to the full element size (unintuitively)

                is_inactive = true;
        
            }

        }

        if (is_small || is_inactive) { 
            
            small_element[k] = SmallElement(k);
            //small_element[k].area = areaCut;
            
            if (domain == 0)
                nb_element_0++;
            else
                nb_element_1++;
        }
    }
}

template <typename Mesh> void MacroElementPartition<Mesh>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);

        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // LOOP OVER FACES
            for (int ifac = 0; ifac < number_of_faces; ++ifac) {

                int ifacn = ifac;
                int kn    = Th.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if (small_or_fat_K[kn] == small)
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;

                if (it != macro_element.end()) { // already exist

                    it->second.add(k, Ks.area);

                } else {

                    const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
                    double areaCut = cutK.measure();

                    macro_element[root_id] = MElement(root_id, areaCut, this);

                    macro_element[root_id].add(k, Ks.area);     // add element to the macroelement
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }
}

template <typename Mesh> void MacroElementPartition<Mesh>::setInnerEdges() {
    // loop over macro elements
    for (auto &[idx_root, MK] : macro_element) {
        // loop over element contained in the macro element
        for (int k = 0; k < MK.size(); ++k) {
            int ki = MK.get_index_element(k);
            // loop over each edges of the element
            for (int ie = 0; ie < Mesh::Element::nea; ++ie) {
                auto [kn, je] = Th.elementAdjacent(ki, ie);
                // if neighbor element belongs to the macro element, we mark as inner edge
                if (MK.containElement(kn)) {

                    // Add only elements with higher index, otherwise the same element
                    // will be added twice
                    if (ki < kn) {
                        MK.inner_edge.push_back(std::make_pair(ki, ie));
                        number_of_stabilized_edges++;
                    }
                }
            }
        }
    }
}

template <typename Mesh> class TimeMacroElement : public GMacro {

  public:
    const ActiveMesh<Mesh> &Th;
    R tol;
    int nb_element_0, nb_element_1;

    TimeMacroElement(const ActiveMesh<Mesh> &Th_, const QuadratureFormular1d &qTime_, const double C_);

    double get_area(int k) const {
        if (isSmall(k)) {
            return getSmallElement(k).area;
        } else if (isRootFat(k)) {
            const auto it(macro_element.find(k));
            return it->second.area_root_;
        } else
            assert(0);
    }

    int number_of_inner_edges();

  private:
    const QuadratureFormular1d &qTime;
    void findSmallElement();
    void createMacroElement();
};

template <typename Mesh>
TimeMacroElement<Mesh>::TimeMacroElement(const ActiveMesh<Mesh> &Th_, const QuadratureFormular1d &qTime_,
                                         const double C_)
    : Th(Th_), qTime(qTime_) {

    double h       = Th[0].lenEdge(1); // catheter of triangle
    double measure = Th[0].measure();  // measure = h^2/2

    nb_element_0 = 0;
    nb_element_1 = 0;
    tol          = 2 * C_ * measure;

    std::cout << "tolerance \t" << tol << std::endl;
    findSmallElement();
    std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
    createMacroElement();
    std::cout << " Macro element created" << std::endl;
}

template <typename Mesh> void TimeMacroElement<Mesh>::findSmallElement() {

    // Iterate over all elements in the active mesh (over the whole time-slab)

    for (int k = 0; k < Th.get_nb_element(); k += 1) {

        if (!Th.isStabilizeElement(k))
            continue; // if the element is not cut or if it doesn't change domain
                      // it doesn't need stabilization

        const typename Mesh::Element &K(Th[k]);

        const int domain = Th.get_domain_element(k);

        // Iterate over the quadrature points in the time-slab In

        bool is_large     = false; // is element large in any quadrature point?
        bool is_small     = true;  // is element small in any quadrature point?
        bool is_inactive  = false; // is element inactive in any quadrature point?
        bool is_never_cut = true;
        int times_small   = 0; // how many times the element is small

        for (int itq = 0; itq < qTime.n; ++itq) {

            Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(k, itq));
            double areaCut = cutK.measure();

            if (Th.isCut(k, itq))
                is_never_cut = false;

            // if element is large and active in itq
            if ((areaCut > tol) && (!Th.isInactive(k, itq))) {
                is_large = true;
                is_small = false;
            }
            // if element is small and active in itq
            else if ((areaCut < tol) && (!Th.isInactive(k, itq))) {
                is_small = true;
                times_small += 1;
            }

            if (Th.isInactive(k, itq))
                is_inactive = true;
        }

        if (!is_large || is_inactive) {
            // if (is_small || is_inactive) {
            // if (times_small >= qTime.n-1 || is_inactive) {
            //  if ((is_inactive && is_never_cut) || !is_large) {
            small_element[k] = SmallElement(k);
            // small_element[k].area = areaCut;
            if (domain == 0)
                nb_element_0++;
            else
                nb_element_1++;
        }
    }
}

template <typename Mesh> void TimeMacroElement<Mesh>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);
        ;
        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // lLOOP OVER FACES
            for (int ifac = 0; ifac < 3; ++ifac) {

                int ifacn = ifac;
                int kn    = Th.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if (small_or_fat_K[kn] == small)
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;

                if (it != macro_element.end()) { // already exist
                    it->second.add(k, std::make_pair(kk, ie), Ks.area);
                } else {

                    const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
                    double areaCut = cutK.measure();

                    macro_element[root_id] = MElement(root_id, areaCut, this);
                    macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }
}

template <typename Mesh> int TimeMacroElement<Mesh>::number_of_inner_edges() {
    int num_of_inner_edges = 0;
    for (auto me = this->macro_element.begin(); me != this->macro_element.end(); ++me) {
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
            num_of_inner_edges += 1;
        }
    }
    return num_of_inner_edges;
}

template <typename Mesh> class TimeMacroElement2 : public GMacro {

  public:
    const ActiveMesh<Mesh> &Th;
    R tol;
    int nb_element_0, nb_element_1;

    TimeMacroElement2(const ActiveMesh<Mesh> &, const QuadratureFormular1d &, const double);

    double get_area(int k) const {
        if (isSmall(k)) {
            return getSmallElement(k).area;
        } else if (isRootFat(k)) {
            const auto it(macro_element.find(k));
            return it->second.area_root_;
        } else
            assert(0);
    }

    int number_of_inner_edges();

  private:
    const QuadratureFormular1d &qTime;
    void findSmallElement();
    void createMacroElement();
};

template <typename Mesh>
TimeMacroElement2<Mesh>::TimeMacroElement2(const ActiveMesh<Mesh> &Th_, const QuadratureFormular1d &qTime_,
                                           const double C_)
    : Th(Th_), qTime(qTime_) {

    double h       = Th[0].lenEdge(1); // catheter of triangle
    double measure = Th[0].measure();  // measure = h^2/2

    nb_element_0 = 0;
    nb_element_1 = 0;
    tol          = 2 * C_ * measure;

    std::cout << "tolerance \t" << tol << std::endl;
    findSmallElement();
    std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
    createMacroElement();
    std::cout << " Macro element created" << std::endl;
}

template <typename Mesh> void TimeMacroElement2<Mesh>::findSmallElement() {

    // Iterate over all elements in the active mesh (over the whole time-slab)

    for (int k = 0; k < Th.get_nb_element(); k += 1) {

        if (!Th.isStabilizeElement(k))
            continue; // if the element is not cut or if it doesn't change domain
                      // it doesn't need stabilization

        const typename Mesh::Element &K(Th[k]);

        const int domain = Th.get_domain_element(k);

        // Iterate over the quadrature points in the time-slab In

        bool is_large           = false; // is element large in any quadrature point?
        bool is_inactive        = false; // is element inactive in any quadrature point?
        int numb_times_inactive = 0;

        for (int itq = 0; itq < qTime.n; ++itq) {

            Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(k, itq));
            double areaCut = cutK.measure();

            if ((areaCut > tol) && (!Th.isInactive(k, itq)))
                is_large = true;
            if (Th.isInactive(k, itq))
                is_inactive = true;
            if (Th.isInactive(k, itq))
                ++numb_times_inactive;
        }

        // if (!is_large || is_inactive) {
        if (!is_large || numb_times_inactive >= 2) {
            small_element[k] = SmallElement(k);
            // small_element[k].area = areaCut;
            if (domain == 0)
                nb_element_0++;
            else
                nb_element_1++;
        }
    }
}

template <typename Mesh> void TimeMacroElement2<Mesh>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);
        ;
        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // lLOOP OVER FACES
            for (int ifac = 0; ifac < 3; ++ifac) {

                int ifacn = ifac;
                int kn    = Th.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if ((small_or_fat_K[kn] == small))
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;

                if (it != macro_element.end()) { // already exist
                    it->second.add(k, std::make_pair(kk, ie), Ks.area);
                } else {

                    const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
                    double areaCut = cutK.measure();

                    macro_element[root_id] = MElement(root_id, areaCut, this);
                    macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }
}

template <typename Mesh> int TimeMacroElement2<Mesh>::number_of_inner_edges() {
    int num_of_inner_edges = 0;
    for (auto me = this->macro_element.begin(); me != this->macro_element.end(); ++me) {
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
            num_of_inner_edges += 1;
        }
    }
    return num_of_inner_edges;
}




template <typename Mesh, typename L> class AlgoimMacro : public GMacro {
  public:
    AlgoimMacro(const ActiveMesh<Mesh> &, const double, L &, const TimeSlab &, const QuadratureFormular1d &);

    const int get_number_of_stabilized_edges() { return number_of_stabilized_edges; }

    const ActiveMesh<Mesh> &Th;

    virtual void findSmallElement();
    void createMacroElement();
    void setInnerEdges();

  protected:
    L phi;
    R tol;
    int nb_element_0, nb_element_1;
    double measure_K;
    int number_of_stabilized_edges;
    int number_of_faces; // number of faces of the mesh element
    const QuadratureFormular1d &qTime;
    const TimeSlab &In;

    
};

template <typename Mesh, typename L>
AlgoimMacro<Mesh, L>::AlgoimMacro(const ActiveMesh<Mesh> &Th_, const double C_, L &phi_, const TimeSlab &In_,
                                  const QuadratureFormular1d &qTime_)
    : Th(Th_), phi(phi_), In(In_), qTime(qTime_), number_of_faces(Th[0].ne), number_of_stabilized_edges(0) {

    nb_element_0 = 0;
    nb_element_1 = 0;
    measure_K    = Th[0].measure();
    tol          = C_ * measure_K;

    // this->findSmallElement(&In_, &qTime_);
    // std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    // std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
    // createMacroElement();
    // setInnerEdges();
    // std::cout << "Macro element created\n";
}

template <typename Mesh, typename L>
void AlgoimMacro<Mesh, L>::findSmallElement() {

    std::cout << "Tolerance: \t" << tol << std::endl;
    std::cout << "|K|: \t" << measure_K << std::endl;

    // Iterate over all elements in the active mesh (over the whole time-slab)
    for (int k = 0; k < Th.get_nb_element(); k += 1) {

        if (!Th.isStabilizeElement(k))
            continue; // if the element is not cut or if it doesn't change domain it doesn't need stabilization

        const typename Mesh::Element &K(Th[k]);

        const int domain = Th.get_domain_element(k);

        // Iterate over the quadrature points in the time-slab In

        bool is_large           = false; // is element large in any quadrature point?
        bool is_small           = false; // is element large in any quadrature point?
        bool is_inactive        = false; // is element inactive in any quadrature point?
        int numb_times_inactive = 0;

        for (int itq = 0; itq < Th.nb_quadrature_time_; ++itq) {

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get current time

            GQuadraturePoint<R1> tq(qTime.at(itq));
            const double t = In.mapToPhysicalElement(tq);
            phi.t          = t;

            // Get quadrature rule for the intersection between the element K and the negative part of the level set
            // function, using quadrature order 5
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, 5);

            double cut_area = q.sumWeights();

            double part = cut_area / measure_K;

            if ((cut_area > tol) && (!Th.isInactive(k, itq))) {
                // std::cout << "LARGE: kb: " << Th.idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: " << k
                //           << ", area_cut: " << cut_area << ", |K|: " << measure_K << ", cut part %: " << part <<
                //           "\n";

                assert(0 < cut_area / measure_K &&
                       cut_area / measure_K <= 1 + 1e-10); // make sure cut area is not equal to element area if cut
                is_large = true;

            } else if ((cut_area <= tol) && (!Th.isInactive(k, itq))) {
                // assert(std::fabs(cut_area - measure_K) >
                //        1e-10); // make sure cut area is not equal to element area if cut
                assert(0 <= cut_area / measure_K &&
                       cut_area / measure_K < 1); // make sure cut area is not equal to element area if cut
                is_small = true;
                // std::cout << "SMALL: kb: " << Th.idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: " << k
                //           << ", area_cut: " << cut_area << ", |K|: " << measure_K << ", cut part %: " << part <<
                //           "\n";
            }

            if (Th.isInactive(k, itq)) {

                assert(cut_area == 0); // make sure element is not cut
                is_inactive = true;
                // std::cout << "INACTIVE. kb: " << Th.idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: "
                //<< k
                // << "\n";
            }

            // if (Th.isInactive(k, itq))
            //     ++numb_times_inactive;
        }

        // if (!is_large || is_inactive) { // method 1
        if (is_small || is_inactive) { // method 2
            // if (!is_large || numb_times_inactive >= 2) {
            (this->small_element)[k] = SmallElement(k);
            // small_element[k].area = areaCut;
            if (domain == 0)
                nb_element_0++;
            else
                nb_element_1++;
        }
    }

    std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
}

template <typename Mesh, typename L> void AlgoimMacro<Mesh, L>::createMacroElement() {

    std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
    std::vector<int> small_or_fat_K(Th.get_nb_element());
    std::vector<std::pair<int, int>> big_element_found;

    for (int i = 0; i < small_or_fat_K.size(); ++i)
        small_or_fat_K[i] = i;
    int ii = 0;
    for (auto it = small_element.begin(); it != small_element.end(); ++it) {
        idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);

        small_or_fat_K[it->second.index] = small;
    }
    int pos = 0;
    while (idx_small_K_temp.size() > 0) {
        int nb_of_small_K_left = idx_small_K_temp.size();
        pos += 1;
        big_element_found.clear();
        for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
            // LOOP OVER SMALL ELEMENTS LEFT

            int k      = idx_small_K_temp[i].first;
            int idx_Ks = idx_small_K_temp[i].second;
            SmallElement &Ks(small_element[idx_Ks]);

            // LOOP OVER FACES
            for (int ifac = 0; ifac < number_of_faces; ++ifac) {

                int ifacn = ifac;
                int kn    = Th.ElementAdj(k, ifacn);
                if (kn == -1)
                    continue;

                if (small_or_fat_K[kn] == small)
                    continue;

                // set position of the small element
                Ks.setChainPosition(pos);
                Ks.setRoot(small_or_fat_K[kn]);
                big_element_found.push_back(std::make_pair(k, kn));

                // find the correonding macro element
                int root_id = small_or_fat_K[kn];
                auto it     = macro_element.find(root_id);
                // for unique edge
                int ie      = (k < kn) ? ifac : ifacn;
                int kk      = (k < kn) ? k : kn;

                if (it != macro_element.end()) { // already exist
                    // it->second.add(k, std::make_pair(kk, ie), Ks.area);
                    it->second.add(k, Ks.area); // add element to macro element
                } else {

                    // const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
                    // double areaCut = cutK.measure();
                    double areaCut = 0.;

                    macro_element[root_id] = MElement(root_id, areaCut, this);
                    // macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
                    macro_element[root_id].add(k, Ks.area); // add element to macro element
                }

                // remove small element from the list
                idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
                break;
            }
        }

        for (int j = 0; j < big_element_found.size(); ++j) {
            int k             = big_element_found[j].first;
            int kn            = big_element_found[j].second;
            small_or_fat_K[k] = small_or_fat_K[kn];
        }
    }

    //std::cout << "Macro element successfully created!\n";
}

template <typename Mesh, typename L> void AlgoimMacro<Mesh, L>::setInnerEdges() {
    // loop over macro elements
    for (auto &[idx_root, MK] : macro_element) {
        // loop over element contained in the macro element
        for (int k = 0; k < MK.size(); ++k) {
            int ki = MK.get_index_element(k);
            // loop over each edges of the element
            for (int ie = 0; ie < Mesh::Element::nea; ++ie) {
                auto [kn, je] = Th.elementAdjacent(ki, ie);
                // if no neighbor or inner edge we do not add the edge as outter edge
                if (MK.containElement(kn)) {

                    // Add only elements with higher index, otherwise the same element
                    // will be added twice
                    if (ki < kn) {
                        MK.inner_edge.push_back(std::make_pair(ki, ie));
                        number_of_stabilized_edges++;
                    }
                }
            }
        }
    }
    std::cout << "Number of stabilized edges: " << number_of_stabilized_edges << "\n";
}






template <typename Mesh, typename L> class AlgoimMacroSurface : public AlgoimMacro<Mesh, L> {
  public:
    AlgoimMacroSurface(const ActiveMesh<Mesh> &, const double, L &, const TimeSlab &,
                const QuadratureFormular1d &);

    void findSmallElement() override;
    
};


template <typename Mesh, typename L>
AlgoimMacroSurface<Mesh, L>::AlgoimMacroSurface(const ActiveMesh<Mesh> &Th_, const double C_, L &phi_, const TimeSlab &In_,
                                  const QuadratureFormular1d &qTime_)
    : AlgoimMacro<Mesh, L>(Th_, C_, phi_, In_, qTime_) {

    this->nb_element_0 = 0;
    // measure_K = length of diagonal of 2D quadrilateral
    this->measure_K = std::sqrt(std::pow((this->Th)[0].mesureBord(0), 2) + std::pow((this->Th)[0].mesureBord(1), 2));
    //measure_K    = Th[0].measure();
    this->tol          = C_ * this->measure_K;


    // findSmallElement();
    // std::cout << this->nb_element_0 << " \t in Gamma" << std::endl;
    // this->createMacroElement();
    // this->setInnerEdges();
    // std::cout << "Macro element created\n";
}


template <typename Mesh, typename L>
void AlgoimMacroSurface<Mesh, L>::findSmallElement() {

    // std::cout << "Tolerance: \t" << tol << std::endl;
    // std::cout << "|E|: \t" << this->measure_K << std::endl;

    // Iterate over all elements in the active mesh (over the whole time-slab)
    for (int k = 0; k < (this->Th).get_nb_element(); k += 1) {

        const typename Mesh::Element &K((this->Th)[k]);

        const int domain = (this->Th).get_domain_element(k);

        // Iterate over the quadrature points in the time-slab In

        bool is_large           = false; // is element large in any quadrature point?
        bool is_small           = false; // is element large in any quadrature point?
        bool is_inactive        = false; // is element inactive in any quadrature point?
        int numb_times_inactive = 0;

        for (int itq = 0; itq < (this->Th).nb_quadrature_time_; ++itq) {

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get current time
            GQuadraturePoint<R1> tq((this->qTime).at(itq));
            const double t = (this->In).mapToPhysicalElement(tq);
            (this->phi).t          = t;

            // Get quadrature rule for the intersection between the element K and the negative part of the level set
            // function, using quadrature order 5
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(this->phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, 5);

            double cut_segment = q.sumWeights();

            double part = cut_segment / (this->measure_K);

            if ((cut_segment > (this->tol)) && (!(this->Th).isInactive(k, itq))) {
                // std::cout << "LARGE: kb: " << (this->Th).idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: " << k
                //           << ", area_cut: " << cut_segment << ", |K|: " << this->measure_K << ", cut part %: " << part <<
                //           "\n";

                assert(0 < cut_segment/(this->measure_K) && cut_segment/(this->measure_K) <= 2.); // cut segment can actually be larger than diagonal of element
                is_large = true;
                
            } else if ((cut_segment <= this->tol) && (!(this->Th).isInactive(k, itq))) {
                // assert(std::fabs(cut_area - measure_K) >
                //        1e-10); // make sure cut area is not equal to element area if cut
                assert(0 <= cut_segment/(this->measure_K) && cut_segment/(this->measure_K) < 1); // make sure cut area is not equal to element area if cut
                is_small = true;
                // std::cout << "SMALL: kb: " << Th.idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: " << k
                //           << ", area_cut: " << cut_area << ", |K|: " << measure_K << ", cut part %: " << part <<
                //           "\n";
            }

            if ((this->Th).isInactive(k, itq)) {
                
                assert(cut_segment == 0); // make sure element is not cut
                is_inactive = true;
                // std::cout << "INACTIVE. kb: " << Th.idx_in_background_mesh_[0][k] << ", itq: " << itq << ", k: "
                //<< k
                // << "\n";
            }

            // if (Th.isInactive(k, itq))
            //     ++numb_times_inactive;
        }

        // if (!is_large || is_inactive) { // method 1
        if (is_small || is_inactive) { // method 2
            // if (!is_large || numb_times_inactive >= 2) {
            (this->small_element)[k] = SmallElement(k);
            // small_element[k].area = areaCut;
            (this->nb_element_0)++;
       
        }
    }

    // std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
    // std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
}


// template <typename Mesh, typename L> void AlgoimMacroSurface<Mesh, L>::createMacroElement() {

//     std::vector<std::pair<int, int>> idx_small_K_temp(small_element.size());
//     std::vector<int> small_or_fat_K(Th.get_nb_element());
//     std::vector<std::pair<int, int>> big_element_found;

//     for (int i = 0; i < small_or_fat_K.size(); ++i)
//         small_or_fat_K[i] = i;
//     int ii = 0;
//     for (auto it = small_element.begin(); it != small_element.end(); ++it) {
//         idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);

//         small_or_fat_K[it->second.index] = small;
//     }
//     int pos = 0;
//     while (idx_small_K_temp.size() > 0) {
//         int nb_of_small_K_left = idx_small_K_temp.size();
//         pos += 1;
//         big_element_found.clear();
//         for (int i = nb_of_small_K_left - 1; i >= 0; --i) {
//             // LOOP OVER SMALL ELEMENTS LEFT

//             int k      = idx_small_K_temp[i].first;
//             int idx_Ks = idx_small_K_temp[i].second;
//             SmallElement &Ks(small_element[idx_Ks]);

//             // LOOP OVER FACES
//             for (int ifac = 0; ifac < number_of_faces; ++ifac) {

//                 int ifacn = ifac;
//                 int kn    = Th.ElementAdj(k, ifacn);
//                 if (kn == -1)
//                     continue;

//                 if (small_or_fat_K[kn] == small)
//                     continue;

//                 // set position of the small element
//                 Ks.setChainPosition(pos);
//                 Ks.setRoot(small_or_fat_K[kn]);
//                 big_element_found.push_back(std::make_pair(k, kn));

//                 // find the correonding macro element
//                 int root_id = small_or_fat_K[kn];
//                 auto it     = macro_element.find(root_id);
//                 // for unique edge
//                 int ie      = (k < kn) ? ifac : ifacn;
//                 int kk      = (k < kn) ? k : kn;

//                 if (it != macro_element.end()) { // already exist
//                     // it->second.add(k, std::make_pair(kk, ie), Ks.area);
//                     it->second.add(k, Ks.area); // add element to macro element
//                 } else {

//                     // const Cut_Part<typename Mesh::Element> cutK(Th.get_cut_part(root_id, 0));
//                     // double areaCut = cutK.measure();
//                     double areaCut = 0.;

//                     macro_element[root_id] = MElement(root_id, areaCut, this);
//                     // macro_element[root_id].add(k, std::make_pair(kk, ie), Ks.area);
//                     macro_element[root_id].add(k, Ks.area); // add element to macro element
//                 }

//                 // remove small element from the list
//                 idx_small_K_temp.erase(idx_small_K_temp.begin() + i);
//                 break;
//             }
//         }

//         for (int j = 0; j < big_element_found.size(); ++j) {
//             int k             = big_element_found[j].first;
//             int kn            = big_element_found[j].second;
//             small_or_fat_K[k] = small_or_fat_K[kn];
//         }
//     }
// }

// template <typename Mesh, typename L> void AlgoimMacroSurface<Mesh, L>::setInnerEdges() {
//     // loop over macro elements
//     for (auto &[idx_root, MK] : macro_element) {
//         // loop over element contained in the macro element
//         for (int k = 0; k < MK.size(); ++k) {
//             int ki = MK.get_index_element(k);
//             // loop over each edges of the element
//             for (int ie = 0; ie < Mesh::Element::nea; ++ie) {
//                 auto [kn, je] = Th.elementAdjacent(ki, ie);
//                 // if no neighbor or inner edge we do not add the edge as outter edge
//                 if (MK.containElement(kn)) {

//                     // Add only elements with higher index, otherwise the same element
//                     // will be added twice
//                     if (ki < kn) {
//                         MK.inner_edge.push_back(std::make_pair(ki, ie));
//                     }
//                 }
//             }
//         }
//     }
// }


#endif

// template<typename Mesh>
// void MacroElementCL<Mesh>::createMacroElementInside(){
//   // loop over the macro Element
//   for(auto it = macro_element.begin(); it != macro_element.end();++it) {
//     const MElement& MK(it->second);
//     int idx_root = MK.get_index_root();
//
//     // if it is NOT cut => we don't do anything
//     if(!Th_.isCut(idx_root)) continue;
//
//     // else we need to find a way to a non cut element
//     std::vector<std::pair<int,int>> path;
//     findPathToInside(idx_root, path);
//
//
//   }
//
//
// }
//
// template<typename Mesh>
// void MacroElementCL<Mesh>::findPathToInside(int
// idx_root,std::vector<std::pair<int,int>>& path){
//
//
//   std::vector<std::pair<int,int>> optimal_path = path;
//   // check element around
//   bool find_a_way_out = false;
//   // Always best to find non cut elements
//   for(int ie=0;ie<Mesh::Element::nea;++ie) {
//
//     int je = ie;
//     int kn = Th_.ElementAdj(idx_root, je);
//     // check if exists
//     if(kn == -1) continue;
//
//     // if it is not cut then its ok
//     if(!Th_.isCut(kn)) {
//       path.push_back(std::make_pair(idx_root, ie));
//       return ;
//     }
//   }
//
//
//   // if only cut elements
//   int nb_of_path_found = 0;
//   for(int ie=0;ie<Mesh::Element::nea;++ie) {
//     std::vector<std::pair<int,int>> found_path = path;
//
//     int je = ie;
//     int kn = Th_.ElementAdj(idx_root, je);
//
//     // check if exists
//     if(kn == -1) continue;
//     // check if it is a small element
//     if(this->isSmall(kn)) continue;
//     // check if it is a root of another element
//     // MAYBE IT IS OK BUT LETS START WITH NOT
//     if(this->isRootFat(kn)) continue;
//
//     found_path.push_back(std::make_pair(idx_root, ie));
//     findPathToInside(kn, found_path);
//
//     if(found_path.size() < optimal_path.size() || nb_of_path_found == 0){
//       optimal_path = found_path;
//       find_a_way_out = true;
//     }
//     nb_of_path_found++;
//   }
//   if(!find_a_way_out) {
//     std::cout << " Root element cannot find his way to a inside element" <<
//     std::endl; assert(0);
//   }
//   path = optimal_path;
// }
