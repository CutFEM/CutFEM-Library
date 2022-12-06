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
   return (1 - P.X()) * R2::KHat[nvedge[i][0]] + P.X() * R2::KHat[nvedge[i][1]];
}

double Triangle2::mesureBord(int i) const {
   return (at(nvedge[i][0]) - at(nvedge[i][1])).norm();
}
