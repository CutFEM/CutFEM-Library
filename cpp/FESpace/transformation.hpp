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
#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

#include "../common/RNM.hpp"
#include "../common/dataStruct2D.hpp"
#include "../common/dataStruct3D.hpp"
#include <map>

template <typename E> class Transformation {

 protected:
   struct Memory {
      // number of element to remember
      const int n = 2;
   };
   static const int N = E::Rd::d;

   struct LocalTransformation {
      KNM<double> DF, invF_t;
      R detDF;
      LocalTransformation() : DF(N, N), invF_t(N, N) {}
   };

   std::map<const E *, LocalTransformation> mapK;
   LocalTransformation *transformation;

   Transformation() {}

   // template<int d> void compute_inverse();
   void initialize(const E &K){}; // new local transfo or give the existing one
   void transform_phi(const E &K, KNMK_<double> &bfMat);
   virtual void init(const E &K) = 0;
};

template <typename E> class PiolaContravariant : public Transformation<E> {

 public:
   PiolaContravariant() {}

 private:
   void init(const E &K) {

      this->initialize(K);
      for (int i = 0; i < this->N; ++i) {
         for (int j = 0; j < this->N; ++j) {
            this->transformation->DF(i, j) = K[j + 1][i] - K[0][i];
         }
      }
      this->transformation->detDF =
          this->transformation->DF(0, 0) * this->transformation->DF(1, 1) -
          this->transformation->DF(0, 1) * this->transformation->DF(1, 0);
   }
};

// template<typename E>
// class Linear_Transformation : public Transformation {
//
//   typedef typename E::Rd Rd;
//   static const int D = Rd::d;
//   static const int nv = E::nv;
//
//   Rd A[Rd::d+1];
//
// public:
//   Linear_Transformation(const E& T) ;
//   void transform_gradient(const Rd phi_hat[nv], KN_<double>& f0i, int d_i);
//
// private:
//   void init();
//
// };
//
//
// // d_i derivative of the ith bf is the  d_i line of
// // the matrix invJ multiply by the ith bf_hat
// // the last input is op_dx, op_dy, op_dz so we have to do -1;
// template<typename E>
// void Linear_Transformation<E>::transform_gradient(const Rd dphi_hat[nv],
// KN_<double>& f0i, int d_i){
//   assert(1<=d_i && d_i<=D);
//   int op = d_i-1;
//   for(int i=0; i<nv;++i){
//     f0i[i] = 0.;
//     for(int j=0;j<D;++j) {
//       f0i[i] += invF_t(op,j)*dphi_hat[i][j];
//     }
//   }
// }
//

#endif
