#include "transformation.hpp"



template<> void Transformation::compute_inverse<2>() {
  this->det = F_t(0,0)*F_t(1,1) - F_t(0,1)*F_t(1,0);
  invF_t(0,0) =  F_t(1,1);
  invF_t(1,1) =  F_t(0,0);
  invF_t(0,1) = -F_t(0,1);
  invF_t(1,0) = -F_t(1,0);
  invF_t /= det;
}
template<> void Transformation::compute_inverse<3>() {
  this->det = F_t(0,0)*F_t(1,1) - F_t(0,1)*F_t(1,0);
  invF_t(0,0) =  F_t(1,1);
  invF_t(1,1) =  F_t(0,0);
  invF_t(0,1) = -F_t(0,1);
  invF_t(1,0) = -F_t(1,0);
  invF_t /= det;
  assert(0);
}

template<> void Linear_Transformation<Quad2>::init() {
  for(int i=0;i<D;++i) {
    for(int j=0;j<D;++j) {
      F_t(i,j) = A[j+1][i] - A[0][i];
    }
  }
  this->compute_inverse<2>();
}
template<> void Linear_Transformation<Hexa>::init() {
  for(int i=0;i<D;++i) {
    for(int j=0;j<D;++j) {
      F_t(i,j) = A[j+1][i] - A[0][i];
    }
  }
  this->compute_inverse<3>();
}

template<> Linear_Transformation<Quad2>::Linear_Transformation(const Quad2& T) : Transformation(2) {
  A[0] = T[0];
  A[1] = T[1];
  A[2] = T[3];
  init();
}
template<> Linear_Transformation<Hexa>::Linear_Transformation(const Hexa& T)   : Transformation(3) {
  A[0] = T[0];
  A[1] = T[1];
  A[2] = T[3];
  A[3] = T[4];
  init();
}
