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
#include "paraview.hpp"

template <>
void Paraview<Mesh2>::ParaviewMesh::build(const ActiveMesh<Mesh> &Vh) {
   return buildNoCut(Vh);
}
template <>
void Paraview<MeshQuad2>::ParaviewMesh::build(const ActiveMesh<Mesh> &Vh) {
   return buildNoCut(Vh);
}
template <>
void Paraview<MeshHexa>::ParaviewMesh::build(const ActiveMesh<Mesh> &Vh) {
   return buildCut(Vh);
}
template <>
void Paraview<Mesh3>::ParaviewMesh::build(const ActiveMesh<Mesh> &Vh) {
   return buildCut(Vh);
}

// template<> void Paraview<Mesh2>::ParaviewMesh::build(const FESpace & Vh,
// Fun_h* levelSet) {
//   return buildNoCut(Vh, levelSet);
// }
// template<> void Paraview<MeshQuad2>::ParaviewMesh::build(const FESpace & Vh,
// Fun_h* levelSet) {
//   return buildNoCut(Vh, levelSet);
// }
// template<> void Paraview<MeshHexa>::ParaviewMesh::build(const FESpace & Vh,
// Fun_h* levelSet) {
//   return buildCut(Vh, levelSet);
// }
// template<> void Paraview<Mesh3>::ParaviewMesh::build(const FESpace & Vh,
// Fun_h* levelSet) {
//   return buildCut(Vh, levelSet);
// }
