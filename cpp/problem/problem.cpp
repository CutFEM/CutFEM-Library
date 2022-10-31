#include "problem.hpp"



template<>
const QuadratureFormular1d& QuadratureOfProblem<Mesh2>::get_quadrature_formular_dK() const {
  return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular1d& QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_dK() const {
  return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d& QuadratureOfProblem<Mesh3>::get_quadrature_formular_dK() const {
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d& QuadratureOfProblem<MeshHexa>::get_quadrature_formular_dK() const {
  return *QF_Quad(order_space_element_quadrature_);
}



template<>
const QuadratureFormular2d& QuadratureOfProblem<Mesh2>::get_quadrature_formular_K() const {
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d& QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_K() const {
  return *QF_Quad(order_space_element_quadrature_);
}
template<>
const QuadratureFormular3d& QuadratureOfProblem<Mesh3>::get_quadrature_formular_K() const {
  return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular3d& QuadratureOfProblem<MeshHexa>::get_quadrature_formular_K() const {
  return *QF_Hexa(order_space_element_quadrature_);
}



// QUADRATURE FOR CUT ELEMENT AND FACE
template<>
const QuadratureFormular3d&  QuadratureOfProblem<Mesh3>::get_quadrature_formular_cutK() const{
  return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular3d&  QuadratureOfProblem<MeshHexa>::get_quadrature_formular_cutK() const{
  return *QF_Simplex<R3>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d&  QuadratureOfProblem<Mesh3>::get_quadrature_formular_cutFace() const{
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d&  QuadratureOfProblem<MeshHexa>::get_quadrature_formular_cutFace() const{
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}

template<>
const QuadratureFormular2d&  QuadratureOfProblem<Mesh2>::get_quadrature_formular_cutK() const{
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular2d&  QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_cutK() const{
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular1d&  QuadratureOfProblem<Mesh2>::get_quadrature_formular_cutFace() const{
  return *QF_Simplex<R1>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular1d&  QuadratureOfProblem<MeshQuad2>::get_quadrature_formular_cutFace() const{
  return *QF_Simplex<R1>(order_space_element_quadrature_);
}
