#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

#include "../common/RNM.hpp"
#include "../common/dataStruct2D.hpp"
#include "../common/dataStruct3D.hpp"



class Transformation {

protected:
  KNM<double> F_t, invF_t;
  R det;

  Transformation(int N) : F_t(N,N), invF_t(N,N) {}


  template<int d> void compute_inverse();

};


template<typename E>
class Linear_Transformation : public Transformation {

  typedef typename E::Rd Rd;
  static const int D = Rd::d;
  static const int nv = E::nv;

  Rd A[Rd::d+1];

public:
  Linear_Transformation(const E& T) ;
  void transform_gradient(const Rd phi_hat[nv], KN_<double>& f0i, int d_i);

private:
  void init();

};


// d_i derivative of the ith bf is the  d_i line of
// the matrix invJ multiply by the ith bf_hat
// the last input is op_dx, op_dy, op_dz so we have to do -1;
template<typename E>
void Linear_Transformation<E>::transform_gradient(const Rd phi_hat[nv], KN_<double>& f0i, int d_i){
  assert(1<=d_i && d_i<=D);
  int op = d_i-1;
  for(int i=0; i<nv;++i){
    f0i[i] = 0.;
    for(int j=0;j<D;++j) {
      f0i[i] += invF_t(op,j)*phi_hat[i][j];
    }
  }


}



#endif
