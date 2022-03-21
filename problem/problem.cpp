#include "problem.hpp"



template<>
const QuadratureFormular2d& QuadratureOfProblem<Mesh2>::get_quadrature_formular_K() const {
  return *QF_Simplex<R2>(order_space_element_quadrature_);
}
template<>
const QuadratureFormular3d& QuadratureOfProblem<Mesh3>::get_quadrature_formular_K() const {
  return *QF_Simplex<R3>(order_space_element_quadrature_);
}

template<>template<>
const QuadratureFormular3d&  QuadratureOfProblem<Mesh3>::get_quadrature_formular_cutK<Tet>() const{
  return *QF_Simplex<R3>(order_space_element_quadrature_);
}
