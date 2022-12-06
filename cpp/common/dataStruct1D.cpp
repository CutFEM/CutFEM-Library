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

#include "dataStruct1D.hpp"

const std::vector<R1> R1::KHat = {R1(0.), R1(1.)};

/// @brief Static mumber for Point1
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvedge = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvface = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::nvhyperFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::edgeOfFace = {
    {}};
template <>
const std::vector<std::vector<int>> GenericElement<DataPoint1>::faceOfEdge = {
    {}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataPoint1>::commonVertOfEdges = {{}};

/// @brief Static mumber for Seg1
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvedge = {{0, 1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvface{{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::nvhyperFace = {
    {0}, {1}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::edgeOfFace = {{}};
template <>
const std::vector<std::vector<int>> GenericElement<DataSeg1>::faceOfEdge = {{}};
template <>
const std::vector<std::vector<int>>
    GenericElement<DataSeg1>::commonVertOfEdges = {{}};
