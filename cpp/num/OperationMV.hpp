#ifndef OPERATION_MV_HPP_
#define OPERATION_MV_HPP_

#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
//#include "lapacke.h"

#include "../util/util.hpp"
// #include "../solver/solver.hpp"
#include "../FESpace/FESpace.hpp"



class GOperationMV {

public :
  static inline int D(int i){
    int op[3] = {op_dx, op_dy,op_dz};
    return op[i];
  }



  static inline int Dx(int i){
    int op[3] = {op_dxx, op_dxy, op_dxz};
    return op[i];
  }
  static inline int Dy(int i){
    int op[3] = {op_dxy, op_dyy, op_dyz};
    return op[i];
  }
  static inline int Dz(int i){
    int op[3] = {op_dxz, op_dyz, op_dzz};
    return op[i];
  }



  // static inline KNM<R> multiply_MatMat(const KNM<R> & A, const KNM<R>& B){
  //   assert(A.M() == B.N());
  //   KNM<R> x(A.N(), B.M()); x = 0.0;
  //   for(int i=0;i<A.N();++i) {
  //     for(int j=0;j<B.M();++j) {
  // 	for(int k=0;k<A.M();++k)
  // 	  x(i,j) += A(i,k)*B(k,j);
  //     }
  //   }
  //   return x;
  // }

  // static inline void projection(Rnm & c, const R* a) {
  //   for(int i = 0; i < c.N(); ++i)
  //     for(int j = 0; j < c.M(); ++j)
  // 	c(i,j) = (i==j) - a[i]*a[j];
  // }

  // static inline KNMK<R> transform_dBF(const KNM<R> & A, const KNMK<R>& w) {
  //   KNMK<R> x(w.N(), w.M(), w.K()); x = 0.0;

  //   for(int e=0; e<w.N();++e)
  //     for(int i=0;i<A.N();++i)
  // 	for(int j=0;j<A.M();++j)
  // 	  x(e,0,D(i)) += A(i,j)*w(e,0,D(j));
  //   return x;
  // }





};


class OperationMV1 : public GOperationMV {

public:
  static inline int D(int i){
    int op[3] = {op_dx, op_dy,op_dz};
    return op[i];
  }



};
class OperationMV2 : public GOperationMV {

  // static constexpr int op2Der[16] = {op_id, op_dx, op_dxx, op_dxxx,
  // 		    op_dy, op_dxy, op_dxxy, op_dxxy,
  // 		    op_dyy, op_dxyy, op_dxyy, op_dxyy,
  // 		    op_dyyy, op_dyyy, op_dyyy, op_dyyy};
public :
  static inline int D(int i, int j){
    int op2Der[16] = {op_id, op_dx, op_dxx, op_dxxx,
		    op_dy, op_dxy, op_dxxy, op_dxxy,
		    op_dyy, op_dxyy, op_dxyy, op_dxyy,
		    op_dyyy, op_dyyy, op_dyyy, op_dyyy};
    return op2Der[i + j*4];
  }

  static inline int D(int i){
    int op[3] = {op_dx, op_dy,op_dz};
    return op[i];
  }


  static inline R2 multiply(const KNM<R> & A, const R* b){
    assert(A.M() == A.N());
    R2 x(0.0,0.0);
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
	x[i] += A(i,j)*b[j];
      }
    }
    return x;
  }

  static inline R2 gradSurface(const KNMK<R>& w, const R2 n, int i) {
    R2 gradS(0,0);
    const R2 n1(1-n.x*n.x  ,-n.x*n.y);
    const R2 n2(-n.x*n.y , 1-n.y*n.y);
    gradS.x = innerProduct(w,n1,i);
    gradS.y = innerProduct(w,n2,i);
    return gradS;
   }

  static inline R2 gradSurface(const R2& w, const R2 n) {
    R2 gradS(0,0);
    const R2 n1(1-n.x*n.x  ,-n.x*n.y);
    const R2 n2(-n.x*n.y , 1-n.y*n.y);
    gradS.x = (w,n1);
    gradS.y = (w,n2);
    return gradS;
   }

   static inline R innerProduct(const KNMK<R>& w, const R2& b, int c,int e) {
     R sum = 0;
     for(int i=0; i<2;++i) sum += w(e,c,D(i))*b[i];
     return sum;
   }

  static inline R innerProduct(const KNMK<R>& w, const R2& b, int e) {
    R sum = 0;
    for(int i=0; i<2;++i) sum += w(e,0,D(i))*b[i];
    return sum;
  }

  static inline R innerProduct(const KNMK<R>& w, const KN<R>& b, int e) {
    R sum = 0;
    for(int i=0; i<2;++i) sum += w(e,0,D(i))*b(i);
    return sum;
  }

  static inline R innerProduct(const KNMK<R>& w, const KN<R>& b, int c,int e) {
    R sum = 0;
    for(int i=0; i<2;++i) sum += w(e,c,D(i))*b(i);
    return sum;
  }

  static inline R innerProduct(const KN<R>& a, const KN<R>& b) {
    R sum = 0;
    for(int i=0; i<2;++i) sum += a(i)*b(i);
    return sum;
  }

  static inline R divergenceSurface(const Rnm& Dval, const R2 n) {
    R divS = 0.;
    const R2 n1(1-n.x*n.x  ,-n.x*n.y);
    const R2 n2(-n.x*n.y , 1-n.y*n.y);
    const R2 du1(Dval(0,0), Dval(0,1));
    const R2 du2(Dval(1,0), Dval(1,1));
    divS = (n1, du1) + (n2, du2);
    return divS;
  }

  static KNM<R> inverse(const KNM<R> & A) {
    assert(A.M() == 2 && A.N() == 2);
    KNM<R> B(2,2);
    B(0,0) = A(1,1);
    B(1,1) = A(0,0);
    B(0,1) = -A(0,1);
    B(1,0) = -A(1,0);
    R det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    B /= det;
    return B;
  }

  static inline R determinant(const KNM<R>& a) {
    R sum = 0;
    sum = a(0,0)*a(1,1) - a(1,0)*a(0,1);
    return std::fabs(sum);
  }

  static inline R crossTerm(const KNMK<R>& w, int ci, int i, int j) {
    return 2 * w(i,ci,D(ci))       * w(j,ci,D(ci)) +
               w(i,ci,D((ci+1)%2)) * w(j,ci,D((ci+1)%2));
  }

  static inline R Dxn(const KNMK<R>& w, R2 normal,  int ci, int i, int j) {
    return 2 * w(i,ci,op_id) * w(j,ci,D(ci))  * normal[ci] +
               w(i,ci,op_id) * w(j,ci,D((ci+1)%2)) * normal[(ci+1)%2];
  }

  static inline R vDnu(const KNMK<R>& w, R2 normal,  int ci, int i, int j) {
    return 2 * w(i,ci,op_id) * w(j,ci,D(ci))  * normal[ci] +
               w(i,ci,op_id) * w(j,ci,D((ci+1)%2)) * normal[(ci+1)%2];
  }

  // static inline R innerProduct(const KNMK<R>& w, const KN<R>& b, int e) {
  //   R sum = 0;
  //   for(int i=0; i<2;++i) sum += w(e,0,D(i))*b[i];
  //   return sum;
  // }



  static inline R normal2Der(const KNMK<R>& w, const R2& normal, int e) {
    assert(0);
    // need to check component in w
    R2 s(0.0,0.0);
    for(int i=0; i<2;++i)  {
      s.x += w(e,i,Dx(i))*normal[i];
      s.y += w(e,i,Dy(i))*normal[i];
    }
    // R sum =0;
    // for(int i=2,j=0;i>=0;--i,++j){
    //   sum += factorial(2)/(factorial(j)*factorial(2-j))
    // 	*pow(normal[0],i)*pow(normal[1],j)*w(e,D(i,j));
    // }

    // std::cout << (s,normal) << "\t" << sum << std::endl;
    // getchar();
    return (s,normal);
  }

  static inline R normal2Der(const KNMK<R>& w, const KNM<R>& DF, const R2& normal, int e) {
    assert(0);
    // need to check component in w
    int M[2][2] = {{op_dxx, op_dxy},{op_dxy, op_dyy}};
    R sum = 0;

    for(int i=0; i<2;++i)
    for(int j=0;j<2;++j)
  	for(int m=0;m<2;++m)
  	  for(int n=0; n<2;++n)
  	    sum += DF(i,m)*DF(j,n)*normal[i]*normal[j]*w(e,i,M[m][n]);

    return sum;
  }

  // static inline R normal3Der(const KNMK<R>& w, const R2& normal, int e) {

  //   R sum =0;
  //   for(int i=3,j=0;i>=0;--i,++j)
  //     sum +=  factorial(3)/(factorial(j)*factorial(3-j))*
  // 	pow(normal[0],i)*pow(normal[1],j)*w(e,0,D(i,j));

  //   return sum;
  // }


  // static inline R innerProduct(const KNMK<R>& w1, const KNMK<R>& w2, int e1, int e2) {
  //   R sum = 0;
  //   for(int i=0; i<2;++i) sum += w1(e1,0,D(i))*w2(e2,0,D(i));
  //   return sum;
  // }

  // static inline R contractProduct(const KNMK<R>& w, const KNM<R>& Pg, int d,int e) {
  //   R sum = 0;
  //   for(int i=0; i<2;++i) {
  //     for(int j=0;j<2;++j) {
  // 	  sum += Pg(d,i) *  w(e,0,D(j)) * Pg(i,j);
  //     }
  //   }
  //   return sum;
  // }




};


class OperationMV3 : public GOperationMV {

//   static constexpr int op2Der[64] = {op_id, op_dx, op_dxx, op_dxxx,
// 				     op_dy, op_dxy, op_dxxy, op_dxxy,
// 				     op_dyy, op_dxyy, op_dxyy, op_dxyy,
// 				     op_dyyy, op_dyyy, op_dyyy, op_dyyy,
// 				     op_dz, op_dxz, op_dxxz, op_dxxz, // z=1, y=0
// 				     op_dyz, op_dxyz, op_dxyz, op_dxyz, // y=1
// 				     op_dyyz, op_dyyz, op_dyyz, op_dyyz, // y=2
// 				     op_dyyz, op_dyyz, op_dyyz, op_dyyz, // y=2
// 				     op_dzz, op_dxzz, op_dxzz, op_dxzz, // z=2, y=0
// 				     op_dyzz, op_dyzz, op_dyzz, op_dyzz, // z=2, y=1
// 				     op_dyzz, op_dyzz, op_dyzz, op_dyzz, // z=2, y=2
// 				     op_dyzz, op_dyzz, op_dyzz, op_dyzz, // z=2, y=2
// 				     op_dzzz, op_dzzz, op_dzzz, op_dzzz, // z = 3
// 				     op_dzzz, op_dzzz, op_dzzz, op_dzzz, // z = 3
// 				     op_dzzz, op_dzzz, op_dzzz, op_dzzz, // z = 3
// 				     op_dzzz, op_dzzz, op_dzzz, op_dzzz  // z = 3
//   };

public :
//   static inline int D(int i, int j, int k){
//     return op2Der[i + j*4 + k*16];
//   }

  static inline R3 multiply(const KNM<R> & A, const R* b){
    assert(A.M() == A.N());
    R3 x(0.0,0.0,0.0);
    for(int i=0;i<A.N();++i) {
      for(int j=0;j<A.M();++j) {
	x[i] += A(i,j)*b[j];
      }
    }
    return x;
  }

  static inline R determinant(const KNM<R>& a) {
    R sum = 0;
    sum = a(0,0)*( a(1,1)*a(2,2) - a(1,2)*a(2,1) )
        - a(0,1)*( a(1,0)*a(2,2) - a(2,0)*a(1,2) )
        + a(0,2)*( a(1,0)*a(2,1) - a(2,0)*a(1,1));
    return sum;
  }


  static KNM<R> inverse(const KNM<R> & A) {
    KNM<R> B(3,3);
    B = A;

    assert(0);

    // int N = 3;
    // int *IPIV = new int[4];
    // LAPACKE_dgetrf(LAPACK_ROW_MAJOR,N,N,B,N,IPIV);
    // LAPACKE_dgetri(LAPACK_ROW_MAJOR,N,B,N,IPIV);
    // delete IPIV;

    return B;
  }

  static inline R innerProduct(const KNMK<R>& w, const R3 b, int c,int e) {
    R sum = 0;
    for(int i=0; i<3;++i) sum += w(e,c,D(i))*b[i];
    return sum;
  }

  static inline R innerProduct(const KNMK<R>& w, const R3 b, int e) {
    R sum = 0;
    for(int i=0; i<3;++i) sum += w(e,0,D(i))*b[i];
    return sum;
  }

  static inline R innerProduct(const KNMK<R>& w, const KN<R>& b, int e) {
    R sum = 0;
    for(int i=0; i<3;++i) sum += w(e,0,D(i))*b(i);
    return sum;
  }
  static inline R innerProduct(const KNMK<R>& w, const KN<R>& b,int c, int e) {
    R sum = 0;
    for(int i=0; i<3;++i) sum += w(e,c,D(i))*b(i);
    return sum;
  }
  static inline R innerProduct(const KN<R>& a, const KN<R>& b) {
    R sum = 0;
    for(int i=0; i<3;++i) sum += a(i)*b(i);
    return sum;
  }

  static inline R3 gradSurface(const KNMK<R>& w, const R3 n, int i) {
    R3 gradS(0,0,0);
    // R3 grad(w(i,0,op_dx),w(i,0,op_dy) ,w(i,0,op_dz));
    const R3 n1(1-n.x*n.x  ,-n.x*n.y, -n.x*n.z);
    const R3 n2(-n.x*n.y , 1-n.y*n.y, -n.y*n.z);
    const R3 n3(-n.x*n.z , -n.y*n.z,   1-n.z*n.z);
    gradS.x = innerProduct(w,n1,i);
    gradS.y = innerProduct(w,n2,i);
    gradS.z = innerProduct(w,n3,i);
    return gradS;
   }

  static inline R3 gradSurface(const R3& w, const R3 n) {
    R3 gradS(0,0,0);
    const R3 n1(1-n.x*n.x  ,-n.x*n.y, -n.x*n.z);
    const R3 n2(-n.x*n.y , 1-n.y*n.y, -n.y*n.z);
    const R3 n3(-n.x*n.z , -n.y*n.z,   1-n.z*n.z);
    gradS.x = (w,n1);
    gradS.y = (w,n2);
    gradS.z = (w,n3);
    return gradS;
   }

  static inline R divergenceSurface(const Rnm& Dval, const R3 n) {
    R divS = 0.;
    const R3 n1(1-n.x*n.x  ,-n.x*n.y, -n.x*n.z);
    const R3 n2(-n.x*n.y , 1-n.y*n.y, -n.y*n.z);
    const R3 n3(-n.x*n.z , -n.y*n.z,   1-n.z*n.z);
    const R3 du1(Dval(0,0), Dval(0,1), Dval(0,2));
    const R3 du2(Dval(1,0), Dval(1,1), Dval(1,2));
    const R3 du3(Dval(2,0), Dval(2,1), Dval(2,2));
    divS = (n1, du1) + (n2, du2) + (n3, du3);
    return divS;
    }

    static inline R crossTerm(const KNMK<R>& w, int ci, int i, int j) {
    return 2 * w(i,ci,D(ci))       * w(j,ci,D(ci)) +
               w(i,ci,D((ci+1)%3)) * w(j,ci,D((ci+1)%3)) +
               w(i,ci,D((ci+2)%3)) * w(j,ci,D((ci+2)%3))
      ;
  }

  static inline R Dxn(const KNMK<R>& w, R3 normal,  int ci, int i, int j) {
    return 2 * w(i,ci,op_id) * w(j,ci,D(ci))  * normal[ci] +
               w(i,ci,op_id) * w(j,ci,D((ci+1)%3)) * normal[(ci+1)%3] +
               w(i,ci,op_id) * w(j,ci,D((ci+2)%3)) * normal[(ci+2)%3];
  }


  static inline R normal2Der(const KNMK<R>& w, const R3& normal, int e) {
    R3 s(0.0,0.0,0.0);
    for(int i=0; i<3;++i)  {
      s.x += w(e,i,Dx(i))*normal[i];
      s.y += w(e,i,Dy(i))*normal[i];
      s.z += w(e,i,Dz(i))*normal[i];
    }
    return (s,normal);
  }

  static inline R normal2Der(const KNMK<R>& w,
  			     const KNM<R>& DF, const R3& normal, int e) {
               assert(0);
               // need to fix component in w
    int M[3][3] = {{op_dxx, op_dxy, op_dxz},
  		   {op_dxy, op_dyy, op_dyz},
  		   {op_dxz, op_dyx, op_dzz}};
    R sum = 0;

    for(int i=0; i<3;++i)
    for(int j=0;j<3;++j)
    for(int m=0;m<3;++m)
    for(int n=0; n<3;++n)
    sum += DF(i,m)*DF(j,n)*normal[i]*normal[j]*w(e,i,M[m][n]);

    return sum;
  }



  // static inline R innerProduct(const KNMK<R>& w, const KN<R>& b, int e) {
  //   R sum = 0;
  //   for(int i=0; i<3;++i) sum += w(e,0,D(i))*b[i];
  //   return sum;
  // }

  // static inline R innerProduct(const KNMK<R>& w1, const KNMK<R>& w2, int e1, int e2) {
  //   R sum = 0;
  //   for(int i=0; i<3;++i) sum += w1(e1,0,D(i))*w2(e2,0,D(i));
  //   return sum;
  // }

  // static inline R contractProduct(const KNMK<R>& w, const KNM<R>& Pg, int d,int e) {
  //   R sum = 0;
  //   for(int i=0; i<3;++i) {
  //     for(int j=0;j<3;++j) {
  // 	  sum += Pg(d,i) *  w(e,0,D(j)) * Pg(i,j);
  //     }
  //   }
  //   return sum;
  // }



  // static inline R normal3Der(const KNMK<R>& w, const R3& normal, int e) {

  //   R sum =0;
  //   // for(int i=3,j=0;i>=0;--i,++j)
  //   //   sum +=  factorial(3)/(factorial(j)*factorial(3-j))*
  //   // 	pow(normal[0],i)*pow(normal[1],j)*w(e,0,D(i,j));
  //   // std::cout << " third normal derivative not implemented in 3D " << std::endl;
  //   return sum;

  // }




};



template<int D> struct TypeOperationMV {};
template<> struct TypeOperationMV<1> { typedef OperationMV1 OperationMV;};
template<> struct TypeOperationMV<2> { typedef OperationMV2 OperationMV;};
template<> struct TypeOperationMV<3> { typedef OperationMV3 OperationMV;};




#endif
