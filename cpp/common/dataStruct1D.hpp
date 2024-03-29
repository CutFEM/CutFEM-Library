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

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef COMMON_DATA_STRUCT_1D_HPP
#define COMMON_DATA_STRUCT_1D_HPP

#include "GenericElement.hpp"

typedef double R;
typedef GenericVertex<R1> Vertex1;

struct DataPoint1 {
    static const int NbOfVertices          = 1;
    static const int NbOfEdges             = 0;
    static const int NbOfFaces             = 0;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = 1;
    static const int NbOfVertexOnHyperFace = 1;
    static const int NbOfRef               = 0;
    static const int NbOfVerticesCut       = 0;
    static const int ParaviewNumCell       = 1;
    static const int nva                   = 0;
    static const int NvOnFace              = 1;
    static const int NbSignPattern         = 3;
    static const int NbNtCut               = 1;
    static const int NbNtPatch             = 1;

    typedef Vertex1 V;
    typedef V::Rd Rd;
    typedef R0 Face;
    static R mesure(V *pv[NbOfVertices]) { return 1.; }

    typedef R0 RdHatBord;
    typedef R0 RdHat;
    static RdHat PBord(const int *nvb, const RdHatBord &P) { return R0(); }
    // static RdHat PBord()  { return R0() ;}
};
class Node1 : public GenericElement<DataPoint1> {
  public:
    Node1(){}; // constructor empty for array
    Rd operator()(const RdHat &Phat) const {
        Rd r = (*(Rd *)vertices[0]);
        return r;
    }
    // std::array<int, 2> index_adjacent_element_;
    // void set_adjacent_element(int k1, int k2) {
    //   index_adjacent_element_ = {k1, k2};
    // }
    // int get_indes_adjacent_element(int i) const {
    //   return index_adjacent_element_[i];
    // }
};

struct DataSeg1 {
    static const int NbOfVertices          = 2;
    static const int NbOfVerticesCut       = 1;
    static const int NbOfFaces             = 0;
    static const int NbOfEdges             = 1;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = NbOfVertices;
    static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
    static const int ParaviewNumCell       = 3;
    static const int NbOfRef               = 2;
    static const int NvOnFace              = 1;
    static const int NbSignPattern         = 9;
    static const int NbNtCut               = 1;
    static const int NbNtPatch             = 1;

    typedef Vertex1 V;
    typedef V::Rd Rd;
    typedef Node1 Face;
    static R mesure(V *pv[NbOfVertices]) { return pv[1]->X() - pv[0]->X(); }
    typedef R0 RdHatBord;
    typedef R1 RdHat;
    static RdHat PBord(const int *nvb, const RdHatBord &P) { return R1(*nvb); }
    // static RdHat PBord(const int * nvb)  { return R1(*nvb) ;}
};
class Seg1 : public GenericElement<DataSeg1> {
  public:
    Seg1(){}; // constructor empty for array

    void Gradlambda(R1 *GradL) const {
        GradL[1] = 1. / measure();
        GradL[0] = -GradL[1];
    }

    R1 toKref(const R1 &P) const {

        const R &A = *vertices[0];
        const R &B = *vertices[1];

        R mes = fabs(B - A);
        return R1((P.X() - A) / mes);
    }
    R1 toReferenceElement(const R1 &P) const {
        const R &A = *vertices[0];
        const R &B = *vertices[1];

        R mes = fabs(B - A);
        return R1((P.X() - A) / mes);
    }
    Rd operator()(const RdHat &Phat) const {
        Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
        for (int i = 1; i < nv; ++i)
            r += Phat[i - 1] * (*(Rd *)vertices[i]);
        return r;
    }
};
class BoundaryPoint1 : public GenericElement<DataPoint1> {
  public:
    BoundaryPoint1(){}; // constructor empty for array

    Rd operator()(const RdHat &Phat) const {
        Rd r = (*(Rd *)vertices[0]);
        return r;
    }
};

template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvedge;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvface;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvhyperFace;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::edgeOfFace;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::faceOfEdge;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg1>::commonVertOfEdges;

#endif
