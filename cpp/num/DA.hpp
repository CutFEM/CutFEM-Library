#ifndef COMMON_DA_HPP_
#define COMMON_DA_HPP_

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <array>

//==========================//
//      Base types          //
//==========================//

template <class U> struct isbase {
   static const bool test = false;
};
template <> struct isbase<bool> {
   static const bool test = true;
};
template <> struct isbase<int> {
   static const bool test = true;
};
template <> struct isbase<double> {
   static const bool test = true;
};
template <> struct isbase<float> {
   static const bool test = true;
};

template <class r_t, bool test = isbase<r_t>::test> struct accessDA_ {
   static inline int get_nb_D() {
      return 1 + accessDA_<typename r_t::v_t>::get_nb_D();
   }
   static inline const double get_val(const r_t &r_) {
      return accessDA_<typename r_t::v_t>::get_val(r_.val);
   }
   static inline void set_val(r_t &r_, double a) {
      accessDA_<typename r_t::v_t>::set_val(r_.val, a);
   }
   static inline void set_D1(r_t &r_, double a, int c) {
      accessDA_<typename r_t::v_t>::set_val(r_.d[c], a);
   }
   static inline void set_derivative(r_t &val, const r_t &di, int i) {
      val.d[i] = di.val;
      accessDA_<typename r_t::v_t>::set_derivative(val.val, di.val, i);
   }
};
template <class r_t> struct accessDA_<r_t, true> {
   static inline int get_nb_D() { return 1; }
   static inline const r_t get_val(const r_t &r_) { return r_; }
   static inline void set_val(r_t &r_, double a) { r_ = a; }
   static inline void set_derivative(r_t &val, const r_t &d, int i) {}
};

template <class r_t> inline double get_Dn(const r_t &r_, int i) {
   return accessDA_<typename r_t::v_t>::get_val(r_.d[i]);
}

template <typename r_t, typename... Args>
inline double get_Dn(const r_t &r_, int first, Args... args) {
   return get_Dn(r_.d[first], args...);
}

template <class V_t, int N = 1> struct DA {
   typedef std::array<V_t, N> d_t;
   typedef V_t v_t;
   typedef DA<V_t, N> this_t;
   const int dim = N;

   v_t val;
   d_t d;
   int nb_D = 1;

   DA() : val(0) {
      for (int i = 0; i < N; ++i)
         d[i] = 0;
   }
   DA(double a) {
      accessDA_<v_t>::set_val(val, a);
      for (int i = 0; i < N; ++i)
         d[i] = 0;
      nb_D = accessDA_<v_t>::get_nb_D();
   }

   DA(double a, int c) {
      accessDA_<v_t>::set_val(val, a);
      nb_D = accessDA_<v_t>::get_nb_D();
      if (N <= c) {
         std::cout << "component of DA construction too high " << std::endl;
         std::exit(EXIT_FAILURE);
      }
      for (int i = 0; i < N; ++i)
         accessDA_<v_t>::set_val(d[i], (i == c));
      for (int i = 0; i < N; ++i)
         accessDA_<V_t>::set_derivative(val, d[i], i);
   }

   DA &operator=(const DA &a) {
      val = a.val;
      d   = a.d;
      return *this;
   }

   DA(const DA &a) : val(a.val), d(a.d) { nb_D = accessDA_<v_t>::get_nb_D(); }

   // DA(const DA &a,const int k) : val(a.val) {
   //   for(int i=0; i<N;++i)
   //     d[i]=(i==k);
   // }
   // DA(const double &a,const int k) : val(a) {
   //   for(int i=0; i<N;++i)
   //     d[i]=(i==k);
   // }
   //
   // DA& operator+ (){return *this;};
   DA &operator=(double);
   // DA& operator= (const DA& a) {
   //   val = a.val;
   //   for(int i=0; i<N;++i) d[i] = a.d[i];
   //   return *this;
   // }
   //
   DA &operator-=(double);
   DA &operator-=(const DA &);
   DA &operator+=(double);
   DA &operator+=(const DA &);
   DA &operator*=(double);
   DA &operator*=(const DA &);
   DA &operator/=(double);
   DA &operator/=(const DA &);
};

// template <class V_t,int N>
// DA<V_t,N> operator + (double x, const DA<V_t,N>& y) {
//   DA<V_t,N> r(y); r.val += x;  return r;
// }
//
// template <class V_t,int N>
// DA<V_t,N> operator + (const DA<V_t,N>& y, double x){
//      DA<V_t,N> r(y);  r.val += x;  return r;
//    }
//
// template <class V_t,int N>
// DA<V_t,N> operator-(const DA<V_t,N>  & x,const DA<V_t,N>  & y){
//   DA<V_t,N> r(x.val-y.val);
//   for(int i=0; i<N;++i)
//     r.d[i] = x.d[i]-y.d[i];
//   return r;
// }
//
// template <class V_t,int N>
// DA<V_t,N> operator - (double x, const DA<V_t,N>& y){
//   DA<V_t,N> r;
//   for(int i=0; i<N;++i) r.d[i]  = (-1)* y.d[i];
//   r.val = x - y.val;   return r;
// }
//
// template <class V_t,int N>
// DA<V_t,N> operator - (const DA<V_t,N>& x, double y){
//   DA<V_t,N>  r(x);  r.val -= y;  return r;
// }
//
template <class V_t, int N>
DA<V_t, N> operator+(const DA<V_t, N> &x, const DA<V_t, N> &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) +
                accessDA_<DA<V_t, N>>::get_val(y));
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i] + y.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}
template <class V_t, int N>
DA<V_t, N> operator+(const DA<V_t, N> &x, const double y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) + y);
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}
template <class V_t, int N>
DA<V_t, N> operator+(const double y, const DA<V_t, N> &x) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) + y);
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator*(const DA<V_t, N> &x, const DA<V_t, N> &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) *
                accessDA_<DA<V_t, N>>::get_val(y));
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i] * y.val + x.val * y.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

// template <class V_t,int N>
// DA<V_t,N> operator + (double x, const DA<V_t,N>& y)
// {   DA<V_t,N> r(y); r.val += x;  return r;}

// template <class V_t,int N>
// DA<V_t,N> operator + (const DA<V_t,N>& y, double x)
// {   DA<V_t,N> r(y);  r.val += x;  return r;}

template <class V_t, int N>
DA<V_t, N> operator-(const DA<V_t, N> &x, const DA<V_t, N> &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) -
                accessDA_<DA<V_t, N>>::get_val(y));
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i] - y.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator-(double x, const DA<V_t, N> &y) {
   DA<V_t, N> r;
   for (int i = 0; i < N; ++i)
      r.d[i] = (-1) * y.d[i];
   r.val = x - y.val;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N> DA<V_t, N> &DA<V_t, N>::operator-=(double y) {
   val -= y;
   return *this;
}

template <class V_t, int N>
DA<V_t, N> &DA<V_t, N>::operator-=(const DA<V_t, N> &y) {
   val -= y.val;
   for (int i = 0; i < N; ++i)
      d[i] -= y.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(val, d[i], i);
   return *this;
}

template <class V_t, int N>
DA<V_t, N> operator-(const DA<V_t, N> &x, double y) {
   DA<V_t, N> r(x);
   r.val -= y;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator*(const DA<V_t, N> &x, const double &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) * y);
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i] * y;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator*(const double &y, const DA<V_t, N> &x) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) * y);
   for (int i = 0; i < N; ++i)
      r.d[i] = y * x.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

// template <class V_t,int N>
// DA<V_t,N> sqrt(const DA<V_t,N>  & x)
// {
//   DA<V_t,N> r(sqrt(x.val));
//   for(int i=0; i<N;++i)
//     r.d[i] = x.d[i]*0.5/r.val;
//   return r;
// }

template <class V_t, int N>
DA<V_t, N> operator/(const DA<V_t, N> &x, const DA<V_t, N> &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) /
                accessDA_<DA<V_t, N>>::get_val(y));
   for (int i = 0; i < N; ++i)
      r.d[i] = (x.d[i] * y.val - x.val * y.d[i]) / (y.val * y.val);
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator/(const DA<V_t, N> &x, const double &y) {
   DA<V_t, N> r(accessDA_<DA<V_t, N>>::get_val(x) / y);
   for (int i = 0; i < N; ++i)
      r.d[i] = x.d[i] / y;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N>
DA<V_t, N> operator/(const double &x, const DA<V_t, N> &y) {
   DA<V_t, N> r(x / accessDA_<DA<V_t, N>>::get_val(y));
   for (int i = 0; i < N; ++i)
      r.d[i] = (-x * y.d[i]) / (y.val * y.val);
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(r.val, r.d[i], i);
   return r;
}

template <class V_t, int N> DA<V_t, N> &DA<V_t, N>::operator=(double y) {
   val = y;
   for (int i = 0; i < N; ++i)
      d[i] = 0;
   return *this;
}

template <class V_t, int N> DA<V_t, N> &DA<V_t, N>::operator+=(double y) {
   val += y;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(val, d[i], i);
   return *this;
}

template <class V_t, int N>
DA<V_t, N> &DA<V_t, N>::operator+=(const DA<V_t, N> &y) {
   val += y.val;
   for (int i = 0; i < N; ++i)
      d[i] += y.d[i];
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(val, d[i], i);
   return *this;
}

template <class V_t, int N> DA<V_t, N> &DA<V_t, N>::operator*=(double y) {
   val *= y;
   for (int i = 0; i < N; ++i)
      d[i] *= y;
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(val, d[i], i);
   return *this;
}
template <class V_t, int N>
DA<V_t, N> &DA<V_t, N>::operator*=(const DA<V_t, N> &y) {
   for (int i = 0; i < N; ++i)
      accessDA_<V_t>::set_derivative(val, d[i], i);
   return *this = *this * y;
}

// template <class V_t,int N>
// DA<V_t,N>& DA<V_t,N>::operator /= (const DA<V_t,N>& y)
// {  return *this = *this / y;}

// template <class V_t,int N>
// DA<V_t,N>& DA<V_t,N>::operator /= (double y)
// { const double inv = 1.0 / y;
//     val *= inv;
//     for(int i=0; i<N;++i) d[i] *= inv;
//     return *this;
// }

// template <class V_t,int N>
// DA<V_t,N> pow (const DA<V_t,N>& x,const int y)
//   {
//     DA<V_t,N> r(1);
//     if(y>=0) for(int i=0;i<y;i++) r*=x;
//     else for(int i=0; i<-y;i++) r/=x;
//     return r;
//   }

// template<int N>
// R laplacian(DA< DA<double,N>,N> & x) {
//   R l = 0;
//   for(int i=0; i<N;++i) l += x.d[i].d[i];

//   return l;
// }

// // template<class V_t,int N, typename A>
// // DA<V_t,N> normalProjection(const DA<V_t,N>& x,const A normal) {
// //   for(int i=0;i<N;++i)
// //     for(int j=0;j<N;++j)

// // }

// template <class V_t,int N>
// DA<V_t,N> sin (const DA<V_t,N>& x)
// {   DA<V_t,N> r; r.val=sin(x.val);
//   for(int i=0; i<N;++i) r.d[i]=x.d[i]*cos(x.val);
//   return r;
// }

// template <class V_t,int N>
// DA<V_t,N> cos (const DA<V_t,N>& x)
// {   DA<V_t,N> r; r.val=cos(x.val);
//   for(int i=0; i<N;++i) r.d[i]=(-1)*x.d[i]*sin(x.val);
//   return r;
// }

// template <class V_t,int N>
// DA<V_t,N> atan (const DA<V_t,N>& x)
// {   DA<V_t,N> r; r.val=atan(x.val);
//   for(int i=0; i<N;++i) r.d[i] = 1 / (1 + x.val*x.val) * x.d[i];
//   return r;
// }

// //------------------------------------------------------
// //------------------------------------------------------
// static R2 gradient(const double* x, DA<double,2> (*fun)(DA<double,2>
// P[2]) ) {
//   DA<double,2> X[2];
//   for(int i=0;i<2;++i) {X[i].val = x[i];X[i].d[i] = 1;}

//   DA<double,2> dx = fun(X);
//   R2 sum;
//   for(int i=0;i<2;++i) sum[i] = dx.d[i];
//   return sum;
// }

// static R3 gradient(const double* x, DA<double,3> (*fun)(DA<double,3>
// P[3]) ) {
//   DA<double,3> X[3];
//   for(int i=0;i<3;++i) {X[i].val = x[i];X[i].d[i] = 1;}

//   DA<double,3> dx = fun(X);
//   R3 sum;
//   for(int i=0;i<3;++i) sum[i] = dx.d[i];
//   return sum;
// }

// static R2 gradient(const R2& x,
// 		   DA<double,2> (*fun)(DA<double,2> P[2], const double),
// 		   const double t ) {
//   DA<double,2> X[2];
//   for(int i=0;i<2;++i) {X[i].val = x[i];X[i].d[i] = 1;}

//   DA<double,2> dx = fun(X, t);
//   R2 sum;
//   for(int i=0;i<2;++i) sum[i] = dx.d[i];
//   return sum;
// }

// static R3 gradient(const R3& x,
// 		   DA<double,3> (*fun)(DA<double,3> P[3], const double),
// 		   const double t ) {
//   DA<double,3> X[3];
//   for(int i=0;i<3;++i) {X[i].val = x[i];X[i].d[i] = 1;}

//   DA<double,3> dx = fun(X, t);
//   R3 sum;
//   for(int i=0;i<3;++i) sum[i] = dx.d[i];
//   return sum;
// }

// static void eval(const R2& x, DA<double,2> (*fun)(DA<double,2> P[2]),
// 		 double& val, R2& grad ) {
//   DA<double,2> X[2];
//   for(int i=0;i<2;++i) {X[i].val = x[i]; X[i].d[i] = 1;}

//   DA<double,2> dx = fun(X);
//   val = dx.val;
//   for(int i=0;i<2;++i) grad[i] = dx.d[i];
// }
// static void eval(const R3& x, DA<double,3> (*fun)(DA<double,3> P[3]),
// 		 double& val, R3& grad ) {
//   DA<double,3> X[3];
//   for(int i=0;i<3;++i) {X[i].val = x[i]; X[i].d[i] = 1;}

//   DA<double,3> dx = fun(X);
//   val = dx.val;
//   for(int i=0;i<3;++i) grad[i] = dx.d[i];
// }

typedef DA<double, 1> Diff_R1;
typedef DA<double, 2> Diff_R2;
typedef DA<double, 3> Diff_R3;
typedef DA<Diff_R1, 1> DDiff_R1;
typedef DA<Diff_R2, 2> DDiff_R2;
typedef DA<Diff_R3, 3> DDiff_R3;
typedef DA<DDiff_R1, 1> DDDiff_R1;
typedef DA<DDiff_R2, 2> DDDiff_R2;
typedef DA<DDiff_R3, 3> DDDiff_R3;
typedef DA<DDDiff_R1, 1> DDDDiff_R1;
typedef DA<DDDiff_R2, 2> DDDDiff_R2;
typedef DA<DDDiff_R3, 3> DDDDiff_R3;

template <int N>
std::ostream &operator<<(std::ostream &f, const DA<double, N> &a) {
   typedef double V_t;
   f << "val = " << accessDA_<DA<V_t, N>>::get_val(a) << std::endl;
   f << "D1 = [  ";
   for (int i = 0; i < N; ++i)
      f << get_Dn(a, i) << "  ";
   f << "] " << std::endl;
   return f;
}
template <int N>
std::ostream &operator<<(std::ostream &f, const DA<DA<double, N>, N> &a) {
   typedef DA<double, N> V_t;
   f << "val = " << accessDA_<DA<V_t, N>>::get_val(a) << std::endl;
   f << "D1 = [  ";
   for (int i = 0; i < N; ++i)
      f << get_Dn(a, i) << "  ";
   f << "] " << std::endl;
   if (a.nb_D > 1) {
      f << "D2 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [ ";
         for (int j = 0; j < N; ++j) {
            f << get_Dn(a, i, j) << "  ";
         }
         f << "] ";
      }
      f << " ] " << std::endl;
   }

   return f;
}
template <int N>
std::ostream &operator<<(std::ostream &f,
                         const DA<DA<DA<double, N>, N>, N> &a) {
   typedef DA<DA<double, N>, N> V_t;
   f << "val = " << accessDA_<DA<V_t, N>>::get_val(a) << std::endl;
   f << "D1 = [  ";
   for (int i = 0; i < N; ++i)
      f << get_Dn(a, i) << "  ";
   f << "] " << std::endl;
   if (a.nb_D > 1) {
      f << "D2 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [ ";
         for (int j = 0; j < N; ++j) {
            f << get_Dn(a, i, j) << "  ";
         }
         f << "] ";
      }
      f << " ] " << std::endl;
   }
   if (a.nb_D > 2) {
      f << "D3 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [  ";
         for (int j = 0; j < N; ++j) {
            f << " [ ";
            for (int k = 0; k < N; ++k) {
               f << get_Dn(a, i, j, k) << "  ";
            }
            f << " ] ";
         }
         f << " ] ";
      }
      f << " ] " << std::endl;
   }
   return f;
}
template <int N>
std::ostream &operator<<(std::ostream &f,
                         const DA<DA<DA<DA<double, N>, N>, N>, N> &a) {
   typedef DA<DA<DA<double, N>, N>, N> V_t;
   f << "val = " << accessDA_<DA<V_t, N>>::get_val(a) << std::endl;
   f << "D1 = [  ";
   for (int i = 0; i < N; ++i)
      f << get_Dn(a, i) << "  ";
   f << "] " << std::endl;
   if (a.nb_D > 1) {
      f << "D2 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [ ";
         for (int j = 0; j < N; ++j) {
            f << get_Dn(a, i, j) << "  ";
         }
         f << "] ";
      }
      f << " ] " << std::endl;
   }
   if (a.nb_D > 2) {
      f << "D3 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [  ";
         for (int j = 0; j < N; ++j) {
            f << " [ ";
            for (int k = 0; k < N; ++k) {
               f << get_Dn(a, i, j, k) << "  ";
            }
            f << " ] ";
         }
         f << " ] ";
      }
      f << " ] " << std::endl;
   }
   if (a.nb_D > 3) {
      f << "D4 = [ ";
      for (int i = 0; i < N; ++i) {
         f << " [  ";
         for (int j = 0; j < N; ++j) {
            f << " [ ";
            for (int k = 0; k < N; ++k) {
               f << " [ ";
               for (int l = 0; l < N; ++l) {
                  f << get_Dn(a, i, j, k, l) << "  ";
               }
               f << " ] ";
            }
            f << " ] ";
         }
         f << " ] ";
      }
      f << " ] " << std::endl;
   }
   return f;
}

#endif
