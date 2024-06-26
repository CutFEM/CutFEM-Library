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
#include "problem.hpp"

template <>
const QuadratureFormular1d &
QuadratureOfProblem<Mesh2>::get_quadrature_formular_dK() const {
   return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular1d &
QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_dK() const {
   return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<Mesh3>::get_quadrature_formular_dK() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<MeshHexa>::get_quadrature_formular_dK() const {
   return *QF_Quad(order_space_element_quadrature_);
}

template <>
const QuadratureFormular2d &
QuadratureOfProblem<Mesh2>::get_quadrature_formular_K() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_K() const {
   return *QF_Quad(order_space_element_quadrature_);
}
template <>
const QuadratureFormular3d &
QuadratureOfProblem<Mesh3>::get_quadrature_formular_K() const {
   return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular3d &
QuadratureOfProblem<MeshHexa>::get_quadrature_formular_K() const {
   return *QF_Hexa(order_space_element_quadrature_);
}

// QUADRATURE FOR CUT ELEMENT AND FACE
template <>
const QuadratureFormular3d &
QuadratureOfProblem<Mesh3>::get_quadrature_formular_cutK() const {
   return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular3d &
QuadratureOfProblem<MeshHexa>::get_quadrature_formular_cutK() const {
   return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<Mesh3>::get_quadrature_formular_cutFace() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<MeshHexa>::get_quadrature_formular_cutFace() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}

template <>
const QuadratureFormular2d &
QuadratureOfProblem<Mesh2>::get_quadrature_formular_cutK() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular2d &
QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_cutK() const {
   return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular1d &
QuadratureOfProblem<Mesh2>::get_quadrature_formular_cutFace() const {
   return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template <>
const QuadratureFormular1d &
QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_cutFace() const {
   return *QF_Simplex<R1>(order_space_element_quadrature_);
}
