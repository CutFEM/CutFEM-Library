
#ifndef RNM_tpl_
#define RNM_tpl_

#include "RNM.hpp"

//   version du 22 nov 1999
//   Voila une debut de class vecteur + matrice
//   les class avec termine avec un _ ne travail que sur
//   un pointeur existant => pas de new et delete de pointeur
//   la class correspondant sans de _ genere les pointeurs
//
//   avec ses classes on peut prendre une ligne
//   ou une colonne d'une matrice
// -----------------------

namespace RNM {

inline double norm(double x) { return x * x; }
inline double norm(float x) { return x * x; }
inline long norm(long x) { return x * x; }
inline int norm(int x) { return x * x; }

} // namespace RNM

template <class R> void MatMul(KNM_<R> &ab, KNM_<R> &a, KNM_<R> &b) {
   // attention ne marche que si les adresses ne sont pas les memes
   long N = a.shapei.n;
   long M = a.shapej.n;
   K_throwassert(a.shapej.n == a.shapei.n);
   K_throwassert(a.shapei.n == ab.shapei.n);
   K_throwassert(b.shapej.n == ab.shapej.n);
   K_throwassert(b.v != ab.v);
   K_throwassert(a.v != ab.v);
   KN_<R> ai = a(0);
   for (long i = 0; i < N; i++, ++ai) {
      KN_<R> bj = b[0];
      for (long j = 0; j < M; j++, ++bj)
         ab(i, j) = (ai, bj);
   }
}

inline std::ostream &operator<<(std::ostream &f, const ShapeOfArray &s) {
   f << s.n;
   const int i10 = 10;
   int prec      = f.precision();
   if (prec < i10)
      f.precision(i10);
   if (s.step != 1)
      f << ":" << s.step;
   if (s.step * s.n != s.next)
      f << " n: " << std::setw(3) << s.next;
   f << ",";
   if (prec < i10)
      f.precision(prec);
   return f;
};

template <class R>
std::ostream &
operator<<(std::ostream &f,
           const KN_<R> &v) { // f <<  " KN_ : " << (ShapeOfArray) v << "
                              // "   <<  (R *) v << " :\n\t"  ;
   f << v.N() << "\t\n\t";
   const int i10 = 10;
   int prec      = f.precision();
   if (prec < i10)
      f.precision(i10);
   for (long i = 0; i < v.N(); i++)
      f << std::setw(3) << v[i] << ((i % 5) == 4 ? "\n\t" : "\t");
   if (prec < i10)
      f.precision(prec);
   return f;
};

template <class R> std::istream &operator>>(std::istream &f, KN_<R> &v) {
   int n;
   char c;
   f >> n;
   ffassert(f.good());
   ffassert(n == v.N());
   while (f.get(c) && (c != '\n' && c != '\r'))
      ; // eat until control (new line

   for (int i = 0; i < n; i++) {
      f >> v[i];
      ffassert(f.good());
   } // modif FH  main 2006
   return f;
}

template <class R> std::istream &operator>>(std::istream &f, KN<R> &v) {
   int n;
   char c;
   f >> n;
   if (v.unset())
      v.init(n);
   std::cout << n << " == " << v.N() << std::endl;
   ffassert(n == v.N());
   while (f.get(c) && (c != '\n' && c != '\r'))
      ; // eat until control (new line

   for (int i = 0; i < n; i++) {
      f >> v[i];
      ffassert(f.good());
   } // modif FH  main 2006
   return f;
}

template <class R>
std::ostream &
operator<<(std::ostream &f,
           const KNM_<R> &v) { // f << " KNM_ "<<v.N()<<"x"<<v.M()<< ": "
                               // << (ShapeOfArray) v
                               //<< " i "  << v.shapei
                               //  << " j "  << v.shapej
                               //  << " " << &v(0,0) << " :\n\t";
   const int i10 = 10;
   int prec      = f.precision();
   if (prec < i10)
      f.precision(i10);
   f << v.N() << ' ' << v.M() /*<< "  n" << v.next<<" :"<< v.shapei.next << ","
                                 << v.shapej.next */
     << "\t\n\t";
   for (long i = 0; i < v.N(); i++) {
      for (long j = 0; j < v.M(); j++)
         f << " " << std::setw(3) << v(i, j);
      f << "\n\t";
   }
   if (prec < i10)
      f.precision(prec);
   return f;
};

template <class R> R KN_<R>::operator,(const KN_<R> &u) const {
   K_throwassert(u.n == n);
   R s = 0;
   R *l(v);
   R *r(u.v);
   for (long i = 0; i < n; i++, l += step, r += u.step)
      s += *l * *r;
   return s;
}

template <class R> R KN_<R>::min() const {
   R minv = v[index(0)];
   for (long i = 1; i < n; i++)
      minv = std::min(minv, v[index(i)]);
   return minv;
}
template <class R> R KN_<R>::max() const {
   R maxv = v[index(0)];
   for (long i = 1; i < n; i++)
      maxv = std::max(maxv, v[index(i)]);
   return maxv;
}

template <class R> R KN_<R>::sum() const {
   R s = v[index(0)];
   for (long i = 1; i < n; i++)
      s += v[index(i)];
   //    std::cout << " sum = " << s << std::endl;
   return s;
}

template <class R> double KN_<R>::norm() const {
   double s = 0.;
   for (long i = 0; i < n; i++)
      s += RNM::norm(v[index(i)]);
   return s;
}

template <class R> double KN_<R>::l2() const {
   double s = 0.;
   for (long i = 0; i < n; i++)
      s += RNM::norm(v[index(i)]);
   return sqrt(s);
}
template <class R> double KN_<R>::l1() const {
   double s = 0.;
   for (long i = 0; i < n; i++)
      s += std::abs(v[index(i)]);
   return (s);
}
template <class R> double KN_<R>::linfty() const {
   double s = 0.;
   for (long i = 0; i < n; i++)
      s = std::max((double)std::abs(v[index(i)]), s);
   return (s);
}
template <class R> double KN_<R>::lp(double p) const {
   if (p == 1.)
      return l1();
   else if (p == 2.)
      return l2();
   else if (p > 1.e10)
      return linfty();
   else {
      double s = 0.;
      for (long i = 0; i < n; i++)
         s += pow(std::max((double)std::abs(v[index(i)]), s), p);
      return pow(s, 1. / p);
   }
}

template <class R> template <class T> long KN_<R>::last(const T &a) const {
   for (long i = n; i-- > 0;)
      if (a(v[index(i)]))
         return i;
   return -1;
}

template <class R> template <class T> long KN_<R>::first(const T &a) const {
   for (long i = 0; i < n; i++)
      if (a(v[index(i)]))
         return i;
   return n;
}

template <class R> void KN_<R>::map(R (*f)(R)) {
   for (long i = 0; i < n; i++) {
      R &x(v[index(i)]);
      x = f(x);
   }
}

template <class R> void KN_<R>::map(R (*f)(const R &)) {
   for (long i = 0; i < n; i++) {
      R &x(v[index(i)]);
      x = f(x);
   }
}

///////////////// definition des operateurs d'affectation
////////////////////////////
#define oper =
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper +=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper -=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper *=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper /=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"

#endif
