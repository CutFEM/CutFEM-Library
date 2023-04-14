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

#include "dataStruct2D.hpp"

const std::vector<R2> R2::KHat = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};
const std::vector<R2> R2::QuadHat = {R2(0., 0.), R2(1., 0.), R2(1., 1.), R2(0., 1.)};

/// @brief Static mumber for Point2
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint2>::nvedge = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint2>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint2>::nvhyperFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint2>::edgeOfFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint2>::faceOfEdge = {
    {}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataPoint2>::commonVertOfEdges = {{}};

/// @brief Static mumber for Seg2
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvedge = {{0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg2>::nvhyperFace = {
    {0}, {1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg2>::edgeOfFace = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg2>::faceOfEdge = {{}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataSeg2>::commonVertOfEdges = {{}};

/// @brief Static mumber for Triangle2
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvedge = {
    {1, 2}, {2, 0}, {0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvface = {
    {0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle2>::nvhyperFace =
    {{1, 2}, {2, 0}, {0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle2>::edgeOfFace =
    {{0, 1, 2}};
template <>
const std::vector<std::vector<int>> GenericElement<DataTriangle2>::faceOfEdge{
    {0}, {0}, {0}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataTriangle2>::commonVertOfEdges = {
        {-1, 2, 1}, {2, -1, 0}, {1, 0, -1}};

/// @brief Static mumber for Quad2
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvedge = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvface = {
    {0, 1, 2, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad2>::nvhyperFace = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad2>::edgeOfFace = {
    {0, 1, 2, 3}};
template <>
const std::vector<std::vector<int>> GenericElement<DataQuad2>::faceOfEdge{
    {0}, {0}, {0}, {0}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataQuad2>::commonVertOfEdges = {
        {-1, 1, -1, 0}, {1, -1, 2, -1}, {-1, 2, -1, 3}, {0, -1, 3, -1}};

/// Some specialization
R2 Triangle2::toKref(const R1 &P, int i) const {
   return (1 - P.X()) * R2::KHat[nvedge[i][0]] + P.X() * R2::KHat[nvedge[i][1]];
}
R2 Triangle2::mapToReferenceElement(const R1 &P, int i) const {
    
    //? When this functions is called from laplace_beltrami_triangle.cpp,
    //? all values i=0,1,2 are called, which means that all edges are called.

    // nvedge[0] = {1, 2}
    // nvedge[1] = {2, 0}
    // nvedge[2] = {0, 1}

    // R2::KHat[0] = {0, 0}
    // R2::KHat[1] = {1, 0}
    // R2::KHat[2] = {0, 1}

    // Example: i = 0
    // p_0^ = (1-P)*{1, 0} + P*{0, 1} = {1-P, P} = point on the hypothenuse
    

    return (1 - P.X()) * R2::KHat[nvedge[i][0]] + P.X() * R2::KHat[nvedge[i][1]];
}

double Triangle2::mesureBord(int i) const {
   return (at(nvedge[i][0]) - at(nvedge[i][1])).norm();
}

double Quad2::mesureBord(int i) const {
   return (at(nvedge[i][0]) - at(nvedge[i][1])).norm();
}

R2 Quad2::mapToReferenceElement(const R1 &P, int i) const {

    // P.X() = quadrature point on line [0,1]
    // map to point on face in 2d reference quadrilateral

    //? When this function is called from laplace_beltrami.cpp, 
    //? i is only 1 or 2, never 0 or 3?
    //? This means that only edge 1 and 2 are iterated through,
    //? and isn't that a problem?

    //? If I switch places of e1 and e2 in "void BaseFEM<M>::addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1, const std::pair<int, int> &e2, const TimeSlab *In, int itq, double cst_time)"
    //? then I get the same behavior but 0 and 3 are called, but never 1 or 2.
    
    //std::cout << "i = " << i << "\n";

    // nvedge[0] = {0, 1}
    // nvedge[1] = {1, 2}
    // nvedge[2] = {2, 3}
    // nvedge[3] = {3, 0}

    // R2::QuadHat[0] = {0, 0}
    // R2::QuadHat[1] = {1, 0}
    // R2::QuadHat[2] = {1, 1}
    // R2::QuadHat[3] = {0, 1}

    return (1 - P.x) * R2::QuadHat[nvedge[i][0]] + P.x * R2::QuadHat[nvedge[i][1]];
}

