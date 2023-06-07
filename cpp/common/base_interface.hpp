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

#ifndef COMMON_BASE_INTERFACE_HPP_
#define COMMON_BASE_INTERFACE_HPP_

#include <iostream>
#include <cassert>
#include <bitset>
#include <memory>
#include "../concept/function.hpp"
#include "RNM.hpp"
#include "Label.hpp"
#include "../num/util.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "cut_method.hpp"

#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif
#include "../FESpace/QuadratureFormular.hpp"

template <int N> struct FaceInterface {};
template <> struct FaceInterface<2> : public SortArray<Uint, 2>, public Label {
    typedef SortArray<Uint, 2> FaceIdx;

    FaceInterface(const Uint &a0, const Uint &a1, int l = 0) : FaceIdx(a0, a1), Label(l) {}
    FaceInterface(Uint *a, int l = 0) : FaceIdx(a), Label(l) {}
    FaceInterface() : FaceIdx(), Label(0) {}
};
template <> struct FaceInterface<3> : public SortArray<Uint, 3>, public Label {
    typedef SortArray<Uint, 3> FaceIdx;

    FaceInterface(const Uint &a0, const Uint &a1, const Uint &a2, int l = 0) : FaceIdx(a0, a1, a2), Label(l) {}
    FaceInterface(Uint *a, int l = 0) : FaceIdx(a), Label(l) {}
    FaceInterface() : FaceIdx(), Label(0) {}
};

template <typeMesh M> class InterfaceLevelSet;

template <typeMesh M> class Interface {

  public:
    using mesh_t = M;
    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;

    static const int nve = Rd::d;
    typedef FaceInterface<nve> Face;
    typedef SortArray<Ubyte, Element::Rd::d + 1> ElementIdx;

    typedef const Face *const_face_iterator;
    typedef Face *face_iterator;
    typedef const Rd *const_vertex_iterator;
    typedef Rd *vertex_iterator;

    const Mesh *backMesh;

    std::vector<Face> faces_;
    std::vector<Rd> vertices_;
    // std::vector<Rd> outward_normal_;
    std::vector<Uint> element_of_face_;
    std::map<int, int> face_of_element_;
    std::vector<Uint> edge_of_node_;

  public:
    Interface(const Mesh &MM) : backMesh(&MM) {}

    virtual Rd normal(int k, std::span<double> x = std::span<double>()) const = 0;
    virtual SignElement<Element> get_SignElement(int k) const                 = 0;

    virtual Partition<Element> get_partition(int k) const                        = 0;
    virtual Partition<typename Element::Face> get_partition_face(const typename Element::Face &face, int k,
                                                                 int ifac) const = 0;

    virtual void cut_partition(Physical_Partition<Element> &local_partition, std::vector<ElementIdx> &new_element_idx,
                               std::list<int> &erased_element, int sign_part) const   = 0;
    virtual R measure(const Face &f) const                                            = 0;
    virtual Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x) const = 0;
    virtual bool isCutFace(int k, int ifac) const                                     = 0;
    virtual bool isCut(int k) const                                                   = 0;

    Rd operator()(const int k, const int i) const { return vertices_[faces_[k][i]]; }
    const Rd &operator()(const int i) const { return vertices_.at(i); }
    const Face &operator[](const int k) const { return faces_.at(k); }
    const Face &getFace(const int k) const { return faces_.at(k); }

    Uint idxElementOfFace(const int k) const { return element_of_face_.at(k); }
    Uint idxFaceOfElement(const int k) const {
        const auto it = face_of_element_.find(k);
        assert(it != face_of_element_.end());
        return it->second;
    }

    Uint nbElement() const { return size(); }

    const Element &get_element(int k) const { return (*backMesh)[k]; }
    const Mesh &get_mesh() const {
        assert(backMesh);
        return *backMesh;
    }
    const_face_iterator face_begin() const { return (faces_.begin()).base(); }
    const_face_iterator face_end() const { return (faces_.end()).base(); }

    const_vertex_iterator vertex_begin() const { return (vertices_.begin()).base(); }
    const_vertex_iterator vertex_end() const { return (vertices_.end()).base(); }

    virtual size_t size() const = 0;

#ifdef USE_MPI
    size_t first_element() const { return MPIcf::first_element(size()); }
    size_t next_element() const { return MPIcf::next_element(size()); }
    size_t last_element() const { return MPIcf::last_element(size()); }

#else
    size_t first_element() const { return 0; }
    size_t next_element() const { return 1; }
    size_t last_element() const { return size(); }
#endif

    virtual R measure(int i) const = 0;

    virtual ~Interface() {}

  private:
    Interface(const Interface &);      // pas de construction par copie
    void operator=(const Interface &); // pas affectation par copy
};

#endif
