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

#ifndef COMMON_DATA_STRUCT_2D_HPP
#define COMMON_DATA_STRUCT_2D_HPP
#include <cassert>
#include <array>

#include "GenericElement.hpp"

typedef GenericVertex<R2> Vertex2;

struct DataPoint2 {
    static const int NbOfVertices          = 1;
    static const int NbOfEdges             = 0;
    static const int NbOfFaces             = 0;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = 1;
    static const int NbOfVertexOnHyperFace = 1;
    static const int NbOfRef               = 0;
    static const int ParaviewNumCell       = 1;
    static const int NbOfVerticesCut       = 0;
    static const int nva                   = 0;
    static const int NvOnFace              = 1;
    static const int NbSignPattern         = 3;
    static const int NbNtCut               = 1;
    static const int NbNtPatch             = 1;

    typedef Vertex2 V;
    typedef V::Rd Rd;
    typedef R0 Face;
    static double mesure(V *pv[NbOfVertices]) { return 1.; }

    typedef R0 RdHatBord;
    typedef R0 RdHat;
    static RdHat PBord(const int *nvb, const RdHatBord &P) { return R0(); }
    // static RdHat PBord()  { return R0() ;}
};

class Node2 : public GenericElement<DataPoint2> {
  public:
    Node2(){}; // constructor empty for array
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

struct DataSeg2 {
    static const int NbOfVertices          = 2;
    static const int NbOfEdges             = 1;
    static const int NbOfFaces             = 0;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = NbOfVertices;
    static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
    static const int NbOfRef               = 2;
    static const int NbOfVerticesCut       = 1;
    static const int ParaviewNumCell       = 3;
    static const int nva                   = 1;
    static const int NvOnFace              = 1;
    static const int NbSignPattern         = 9;
    static const int NbNtCut               = 1;
    static const int NbNtPatch             = 1;

    typedef Vertex2 V;
    typedef V::Rd Rd;
    static double mesure(V *pv[NbOfVertices]) { return R2(*pv[0], *pv[1]).norme(); }
    typedef R1 RdHat;
    typedef R0 RdHatBord;
    typedef Node2 Face;
    static RdHat PBord(const int *nvb, const RdHatBord &P) { return RdHat(*nvb); }
    // static RdHat PBord(const int * nvb)  { return RdHat(*nvb) ;}

    // static const int (* const nvface)[3];// = nvfaceSeg ;
    // static const int (* const nvedge)[2];//  = nvedgeSeg;
};

class Edge2 : public GenericElement<DataSeg2> {
  public:
    Edge2(){}; // constructor empty for array
    Rd operator()(const RdHat &Phat) const {
        Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
        for (int i = 1; i < nv; ++i)
            r += Phat[i - 1] * (*(Rd *)vertices[i]);
        return r;
    }
    // std::array<int, 2> index_adjacent_element_;
    // void set_adjacent_element(int k1, int k2) {
    //   index_adjacent_element_ = {k1, k2};
    // }
    // int get_index_adjacent_element(int i) const {
    //   return index_adjacent_element_[i];
    // }
};

class BoundaryEdge2 : public GenericElement<DataSeg2> {
  public:
    BoundaryEdge2(){}; // constructor empty for array
    Rd operator()(const RdHat &Phat) const {
        Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
        for (int i = 1; i < nv; ++i)
            r += Phat[i - 1] * (*(Rd *)vertices[i]);
        return r;
    }
};

struct DataTriangle2 {
    static const int NbOfVertices          = 3;
    static const int NbOfFaces             = 1;
    static const int NbOfEdges             = 3;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = NbOfVertices;
    static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
    static const int NbOfVerticesCut       = 2;
    static const int NbOfRef               = 4;
    static const int ParaviewNumCell       = 5;
    static const int NvOnFace              = 3;
    static const int NbSignPattern         = 27;
    static const int NbNtCut               = 1;
    static const int NbNtPatch             = 1;

    typedef Vertex2 V;
    typedef V::Rd Rd;
    typedef Edge2 Face;
    static double mesure(V *pv[NbOfVertices]) { return det(*pv[0], *pv[1], *pv[2]) * 0.5; }
    typedef R2 RdHat;
    typedef R1 RdHatBord;
    static RdHat PBord(const int *nvb, const RdHatBord &P) {
        return RdHat::KHat[nvb[0]] * (1 - P.X()) + R2::KHat[nvb[1]] * (P.X());
    }
};

class Triangle2 : public GenericElement<DataTriangle2> {
  public:
    typedef Edge2 Face;
    typedef Triangle2 TypeCutElement;
    Triangle2(){}; // constructor empty for array
    Triangle2(Vertex *v0, int *iv, int r = 0, double mss = globalVariable::UnSetMesure) {
        this->set(v0, iv, r, mss);
    }; // constructor empty for array

    /**
     * @brief Computes grad(phi_i) for i=1,2
     *
     * @param i Edge index
     * @return R2 If i is a horizontal edge, return [0, +- hx]/2A. If i is a vertical edge, return [+- hy, 0]/2A.
     * @note The sign depends on the orientation of the triangle
     */
    R2 H(int i) const {
        assert(i >= 0 && i < 3);
        R2 E = Edge(i); // R2 vector from vertex 0 to vertex 1 of edge i
        // Edge(0) = Vector from edge (1, 2)
        // Edge(1) = Physical coordinates of edge (2, 0)
        // Edge(2) = Physical coordinates of edge (0, 1)

        // std::cout << "i = " << i << std::endl;
        // std::cout << "Edge(i) = " << E << std::endl;
        // std::cout << "Edge(0) = " << Edge(0) << std::endl;
        // std::cout << "E.perp() = " << E.perp() << std::endl;
        // getchar();

        // E = (x, y)
        // E.perp() = (-y, x)

        return E.perp() / (2. * this->measure());
    }

    /**
     * @brief Compute gradients of the basis functions on the physical element
     *
     * @param GradL Empty array of empty R2 objects.
     * @return Fill GradL such that it becomes
     * GradL = (grad(phi_0), grad(phi_1), grad(phi_2))
     */
    void Gradlambda(R2 *GradL) const {

        GradL[1] = H(1);                 // grad(phi_1) = [+-h_y, 0]/2A
        GradL[2] = H(2);                 // grad(phi_2) = [0, -+h_x]/2A
        GradL[0] = -GradL[1] - GradL[2]; // grad(phi_0) = [-+h_y, +-h_x]/2A
    }

    R2 toKref(const R2 &P) const {
        double l[3];
        const R2 &A = *vertices[0];
        const R2 &B = *vertices[1];
        const R2 &C = *vertices[2];

        R2 PA(P, A), PB(P, B), PC(P, C);
        l[0] = 0.5 / mes * ((PB ^ PC));
        l[1] = 0.5 / mes * ((PC ^ PA));
        l[2] = 1 - l[0] - l[1];
        return l[0] * R2::KHat[0] + l[1] * R2::KHat[1] + l[2] * R2::KHat[2];
    }
    R2 mapToReferenceElement(const R2 &P) const {
        double l[3];
        const R2 &A = *vertices[0];
        const R2 &B = *vertices[1];
        const R2 &C = *vertices[2];

        R2 PA(P, A), PB(P, B), PC(P, C);
        l[0] = 0.5 / mes * ((PB ^ PC));
        l[1] = 0.5 / mes * ((PC ^ PA));
        l[2] = 1 - l[0] - l[1];
        return l[0] * R2::KHat[0] + l[1] * R2::KHat[1] + l[2] * R2::KHat[2];
    }
    Rd operator()(const RdHat &Phat) const {
        Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
        for (int i = 1; i < nv; ++i)
            r += Phat[i - 1] * (*(Rd *)vertices[i]);
        return r;
    }

    R2 centroid() const { return 1. / 3 * ((*vertices[0]) + (*vertices[1]) + (*vertices[2])); }

    R2 toKref(const R1 &P, int i) const;
    R2 mapToReferenceElement(const R1 &P, int i) const;
    double mesureBord(int i) const;
};

struct DataQuad2 {
    static const int NbOfVertices          = 4;
    static const int NbOfFaces             = 1;
    static const int NbOfEdges             = 4;
    static const int NbOfTet               = 0;
    static const int NbOfAdjElem           = 4;
    static const int NbOfVertexOnHyperFace = 2;
    static const int NbOfVerticesCut       = 2;
    static const int NvOnFace              = 4;
    static const int NbOfRef               = 4;
    static const int ParaviewNumCell       = 9;
    static const int NbSignPattern         = 81;
    static const int NbNtCut               = 2;
    static const int NbNtPatch             = 2;

    typedef Vertex2 V;
    typedef V::Rd Rd;
    typedef Edge2 Face;
    static double mesure(V *pv[NbOfVertices]) { return det(*pv[0], *pv[1], *pv[2]); }
    typedef R2 RdHat;
    typedef R1 RdHatBord;
    static RdHat PBord(const int *nvb, const RdHatBord &P) {
        return RdHat::KHat[nvb[0]] * (1 - P.X()) + R2::KHat[nvb[1]] * (P.X());
    }
};
/**
 * @brief Class for 2D quadrilateral mesh elements
 *
 */
class Quad2 : public GenericElement<DataQuad2> {
  public:
    typedef Edge2 Face;
    typedef Triangle2 TypeCutElement;

    Quad2(){}; // constructor empty for array
    Quad2(Vertex *v0, int *iv, int r = 0, double mss = globalVariable::UnSetMesure) {
        this->set(v0, iv, r, mss);
    }; // constructor empty for array

    /**
     * @brief Map coordinates in reference element to physical element
     *
     * @param Phat Coordinates in reference element.
     * @note We map a point Phat in the reference element to its
     * corresponding point P in the physical element.
     * @return Rd Coordinates in physical element.
     */
    Rd operator()(const RdHat &Phat) const override {
        // Phat = (Phat.x, Phat.y)
        // P = (x, y) where
        // x = x0 + h_x * xhat
        // y = y0 + h_y * yhat

        R x = vertices[0]->x + (vertices[1]->x - vertices[0]->x) * Phat.x;
        R y = vertices[0]->y + (vertices[3]->y - vertices[0]->y) * Phat.y;

        return Rd(x, y);
    }

    // R2 H(int i) const { assert(i>=0 && i <3);
    //   R2 E=Edge(i);return E.perp()/(2.*this->mesure());
    // } // heigth

    /**
     * @brief Compute gradients of the basis functions on the physical 2D quadrilateral
     *
     * @param GradL Empty array of empty R2 objects.
     * @return Fill GradL such that it becomes
     * GradL = (grad(phi_0), grad(phi_1), grad(phi_2), grad(phi_3))
     */
    // void Gradlambda(R2 *GradL, const R2 &P) const {
    void Gradlambda(R2 *GradL) const {
        // We are working on the physical elements, but when we compute the basis functions
        // and its gradients, we use the uniformity required between the physical element
        // and the reference element, so we equate e.g. (y-y0)/hy = yhat, (y1 - y)/hy = 1-yhat.
        // In this function, however, P = (x, y) so we map it to the reference element first.

        double hx = (vertices[1]->x - vertices[0]->x);
        double hy = (vertices[3]->y - vertices[0]->y);

        // R2 PHat(this->(P));    // P will be physical coordinate, so we want to map it to reference coordinates PHat
        R2 PHat(0.5, 0.5);

        // Assert reference element
        assert(0. <= PHat.x && PHat.x <= 1.);
        assert(0. <= PHat.y && PHat.y <= 1.);

        GradL[0] = R2(-1. * (1. - PHat.y) / hx, -1. * (1 - PHat.x) / hy);
        GradL[1] = R2((1. - PHat.y) / hx, -PHat.x / hy);
        GradL[2] = R2(PHat.y / hx, PHat.x / hy);
        GradL[3] = R2(-PHat.y / hx, (1. - PHat.x) / hy);

        // GradL[0] = R2(1./4,  1./4);
        // GradL[1] = R2(1./4 , 1./4);
        // GradL[2] = R2(1./4,  1./4);
        // GradL[3] = R2(1./4,  1./4);

        // GradL[0] = R2(0., 0.);
        // GradL[1] = R2(0., 0.);
        // GradL[2] = R2(0., 0.);
        // GradL[3] = R2(0., 0.);
    }

    R2 toKref(const R2 &P) const {

        // same as mapToReferenceElement
        const R2 &A = *vertices[0];
        const R2 &B = *vertices[1];
        const R2 &C = *vertices[3];

        R2 AB(A, B), AC(A, C);
        R2 AP(A, P);
        double c1 = (AP, AB) / AB.norme2();
        double c2 = (AP, AC) / AC.norme2();
        return R2(c1, c2);
    }

    /**
     * @brief Map point in physical element to reference element
     *
     * @param P Point in physical element
     * @return R2 Corresponding element in reference element
     */
    R2 mapToReferenceElement(const R2 &P) const {
        const R2 &A = *vertices[0];
        const R2 &B = *vertices[1];
        const R2 &C = *vertices[3];

        R2 AB(A, B), AC(A, C);
        R2 AP(A, P);
        double c1 = (AP, AB) / AB.norme2();
        double c2 = (AP, AC) / AC.norme2();
        return R2(c1, c2);
    }

    /**
     * @brief Map point on unit interval to point on face of reference element
     *
     * @param P Point between 0 and 1
     * @param i Edge index of reference element
     * @return R2 Coordinates of P on edge i in the reference element
     */
    R2 mapToReferenceElement(const R1 &P, int i) const;

    // R2 centroid() const {
    //   return 1./3*((*vertices[0])+(*vertices[1])+(*vertices[2]));
    // }

    // R2 toKref(const R1& P, int i) const;
    // double mesureBord(int i) const ;

    double mesureBord(int i) const;
};

// Forward declaration
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvedge;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvface;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvhyperFace;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::edgeOfFace;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::faceOfEdge;
template <> const std::vector<std::vector<int>> GenericElement<DataSeg2>::commonVertOfEdges;

template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvedge;
template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvface;
template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvhyperFace;
template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::edgeOfFace;
template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::faceOfEdge;
template <> const std::vector<std::vector<int>> GenericElement<DataTriangle2>::commonVertOfEdges;

template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvedge;
template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvface;
template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvhyperFace;
template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::edgeOfFace;
template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::faceOfEdge;
template <> const std::vector<std::vector<int>> GenericElement<DataQuad2>::commonVertOfEdges;

#endif
