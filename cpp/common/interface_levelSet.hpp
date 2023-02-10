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

template <typename M> class InterfaceLevelSet : public Interface<M> {

    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    static const int nve = Rd::d;
    typedef FaceInterface<nve> Face;
    typedef SortArray<Ubyte, Element::Rd::d + 1> ElementIdx;

    KN<byte> ls_sign;
    KN<double> ls_;

  public:
    template <typename Fct> InterfaceLevelSet(const Mesh &MM, const Fct &lss, int label = 0);

    SignElement<Element> get_SignElement(int k) const;
    Partition<Element> get_partition(int k) const;
    Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k, int ifac) const;
   

    bool isCutFace(int k, int ifac) const;

    void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                       std::list<int> &erased_element, int sign_part) const;

    R measure(const Face &f) const;

  private:
    void make_patch(int label);

    const Face make_face(const typename RefPatch<Element>::FaceIdx &ref_tri, const typename Mesh::Element &K,
                         const double lset[Element::nv], int label);

    Rd make_normal(const typename Mesh::Element &K, const double lset[Element::nv]);

    // Rd get_intersection_node(int k, const Rd A, const Rd B) const;

    Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const;
    
};

#include "interface_levelSet.tpp"

#endif
