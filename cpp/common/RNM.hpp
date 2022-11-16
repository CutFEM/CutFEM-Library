#ifndef KNM_H_
#define KNM_H_

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <fstream>
#include <list>

// using namespace std;

// inline void Check_Kn(const char *str, const char *file, int line) {
//    std::std::cerr << "CHECK_KN: " << str << " in file: " << file << ", line "
//    << line
//              << std::endl;
//    assert(0);
//    abort();
// }

// #define K_bigassert(i)                                                         \
//    if (!(i))                                                                   \
//       Check_Kn(#i, __FILE__, __LINE__);
// #define assert( 0&&i) Check_Kn(i, __FILE__, __LINE__);
// #ifdef CHECK_KN
// #define assert(i) \
//    if (!(i)) \
//       Check_Kn(#i, __FILE__, __LINE__);
// #else
// #define assert(i)
// #endif

namespace RNM {
template <class T> inline T Min(const T &a, const T &b) {
   return a < b ? a : b;
}
template <class T> inline T Max(const T &a, const T &b) {
   return a > b ? a : b;
}
template <class T> inline T Abs(const T &a) { return a < 0 ? -a : a; }

template <class T> inline void Exchange(T &a, T &b) {
   T c = a;
   a   = b;
   b   = c;
}
template <class T> inline T Max(const T &a, const T &b, const T &c) {
   return Max(Max(a, b), c);
}
template <class T> inline T Min(const T &a, const T &b, const T &c) {
   return Min(Min(a, b), c);
}
} // namespace RNM

template <class R> class KNMK_;
template <class R> class KNM_;
template <class R> class KN_;

template <class R> class KNMK;
template <class R> class KNM;
template <class R> class KN;

template <class R> class Add_KN_;
template <class R> class Sub_KN_;
template <class R> class Mulc_KN_;
template <class R> class Add_Mulc_KN_;
template <class R> class Mul_KNM_KN_;
template <class R> class DotMul_KN_;

#ifndef ffassert
#define ffassert assert
#endif

class SubArray {
 public:
   const long n, step, start;
   explicit SubArray(long nn, long sta = 0, long s = 1)
       : n(nn), step(s), start(sta) {}
   // SubArray(const FromTo& ft) : n(ft.to-ft.from+1),step(1),start(ft.from) {}
   // SubArray(const ShapeOfArray & ); // all
   long end() const { return start + step * n; }
   long last() const { return start + step * (n - 1); }
   long len1() const { return step * (n - 1); }
};

class ShapeOfArray {
 protected:
 public:
   long n;    //   n  nb of item
   long step; //   step  nb of between 2 item
   long next; //  the   next array of same type in matrix for subarray
              // by default  no next ( in case of KN, and KNM  -next is
              // a counter of destruction  (use in frefem++)
   ShapeOfArray(const ShapeOfArray &s, long nn) : n(s.n), step(s.n), next(nn) {}
   ShapeOfArray(long nn) : n(nn), step(1), next(-1) {}
   ShapeOfArray(long nn, long s) : n(nn), step(s), next(-1) {}
   ShapeOfArray(long nn, long s, long nextt) : n(nn), step(s), next(nextt) {}
   ShapeOfArray(const ShapeOfArray &old, long stepo, long start)
       : n(old.n - start), step(old.step * stepo), next(old.next) {
      assert(n >= 0);
   }
   ShapeOfArray(const ShapeOfArray &old, const SubArray &sub)
       : n(sub.n), step(old.step * sub.step), next(old.next) {
      assert((sub.last()) * old.step <= old.last());
   }

   long end() const { return n * step; }
   long last() const { return (n - 1) * step; }
   long constant() const { return step == 0; }
   long index(long k) const {
      assert((k >= 0) && ((k < n) || !step));
      return step * k;
   }
   ShapeOfArray operator*(long stepp) const {
      return ShapeOfArray(n, step * stepp, next);
   }
   bool SameShape(const ShapeOfArray &a) const {
      return !step || !a.step || a.n == n;
   }
   long N(const ShapeOfArray &a) { return step ? n : a.n; } // size of 2 shape

   long operator[](long k) const {
      // if( k<0 || ( k<n && !step) )
      //   std::cout << "k,n,step=" << k << " " << n << " " << step <<
      //   std::endl;
      assert((k >= 0) && ((k < n) || !step));
      return step * k;
   }

   void init(long nn, long s = 1, long nextt = -1) {
      n    = nn;
      step = s;
      next = nextt;
   }
};

std::ostream &operator<<(std::ostream &f, const ShapeOfArray &s);
template <class R> std::ostream &operator<<(std::ostream &f, const KN_<R> &v);
template <class R> std::istream &operator>>(std::istream &f, KN_<R> &v);
template <class R> std::istream &operator>>(std::istream &f, KN<R> &v);

template <class R> class KN_ : public ShapeOfArray {
 protected:
   R *v;

 public:
   typedef R K; // type of
   long N() const { return n; }
   bool unset() const { return !v; }
   void set(R *vv, int nn, int st = 1, int nx = -1) {
      v    = vv;
      n    = nn;
      step = st;
      next = nx;
   }
   long size() const { return step ? n * step : n; }
   operator R *() const { return v; }

   KN_(const KN_<R> &u) : ShapeOfArray(u), v(u.v) {}
   KN_(const KN_<R> &U, const SubArray &sa)
       : ShapeOfArray(U, sa), v(U.v + U.index(sa.start)) {}

   KN_ operator()(const SubArray &sa) const {
      return KN_(*this, sa);
   } // sub array

   friend KN_<R> sub_array(const KN_<R> u, int i0, int l) {
      KN_<R> sub_u(u(SubArray(l, i0)));
      return sub_u;
   }

   R &operator[](long i) const { return v[index(i)]; }
   R &operator()(long i) const { return v[index(i)]; }
   R &operator[](int i) const { return v[index(i)]; }
   R &operator()(int i) const { return v[index(i)]; }
   R &operator[](unsigned int i) const { return v[index(i)]; }
   R &operator()(unsigned int i) const { return v[index(i)]; }

   R operator,(const KN_<R> &v) const; // dot  product

   KN_ &operator=(const KN_<R> &u);
   KN_ &operator+=(const KN_<R> &u);
   KN_ &operator-=(const KN_<R> &u);

   KN_ &operator*=(const KN_<R> &u);
   KN_ &operator/=(const KN_<R> &u);

   KN_ &operator=(R a);
   KN_ &operator+=(R a);
   KN_ &operator-=(R a);
   KN_ &operator/=(R a);
   KN_ &operator*=(R a);

   KN_ &operator=(R *a) { return operator=(KN_<R>(a, n)); }
   KN_ &operator+=(R *a) { return operator+=(KN_<R>(a, n)); }
   KN_ &operator-=(R *a) { return operator-=(KN_<R>(a, n)); }
   KN_ &operator*=(R *a) { return operator*=(KN_<R>(a, n)); }
   KN_ &operator/=(R *a) { return operator/=(KN_<R>(a, n)); }

   R min() const;
   R max() const;
   R sum() const;
   double norm() const;
   double l2() const;
   double l1() const;
   double linfty() const;
   double lp(double p) const;

   template <class T> long last(const T &) const;
   template <class T> long first(const T &) const;

   void map(R (*f)(R)); // apply the f fonction a all element of the array
   void map(
       R (*f)(const R &)); // apply the f fonction a all element of the array

   KN_ &operator=(const Add_KN_<R> &u);
   KN_ &operator+=(const Add_KN_<R> &u);
   KN_ &operator-=(const Add_KN_<R> &u);
   KN_ &operator*=(const Add_KN_<R> &u);
   KN_ &operator/=(const Add_KN_<R> &u);

   KN_ &operator=(const Sub_KN_<R> &u);
   KN_ &operator-=(const Sub_KN_<R> &u);
   KN_ &operator+=(const Sub_KN_<R> &u);
   KN_ &operator*=(const Sub_KN_<R> &u);
   KN_ &operator/=(const Sub_KN_<R> &u);

   KN_ &operator=(const Mulc_KN_<R> &u);
   KN_ &operator+=(const Mulc_KN_<R> &u);
   KN_ &operator-=(const Mulc_KN_<R> &u);
   KN_ &operator*=(const Mulc_KN_<R> &u);
   KN_ &operator/=(const Mulc_KN_<R> &u);

   KN_ &operator=(const Add_Mulc_KN_<R> &u);
   KN_ &operator+=(const Add_Mulc_KN_<R> &u);
   KN_ &operator-=(const Add_Mulc_KN_<R> &u);
   KN_ &operator*=(const Add_Mulc_KN_<R> &u);
   KN_ &operator/=(const Add_Mulc_KN_<R> &u);

   KN_ &operator=(const Mul_KNM_KN_<R> &u);
   KN_ &operator+=(const Mul_KNM_KN_<R> &u);
   KN_ &operator-=(const Mul_KNM_KN_<R> &u);
   KN_ &operator*=(const Mul_KNM_KN_<R> &u);
   KN_ &operator/=(const Mul_KNM_KN_<R> &u);

   KN_ &operator=(const DotMul_KN_<R> &u);
   KN_ &operator+=(const DotMul_KN_<R> &u);
   KN_ &operator-=(const DotMul_KN_<R> &u);
   KN_ &operator*=(const DotMul_KN_<R> &u);
   KN_ &operator/=(const DotMul_KN_<R> &u);

   friend std::ostream &operator<< <R>(std::ostream &f, const KN_<R> &v);

   KN_(R *u, const ShapeOfArray &s) : ShapeOfArray(s), v(u) {}
   KN_(R *u, long nn, long s) : ShapeOfArray(nn, s), v(u) {}
   KN_(R *u, long nn, long s, long nextt) : ShapeOfArray(nn, s, nextt), v(u) {}
   KN_(R *u, long nn) : ShapeOfArray(nn), v(u) {}

 private:
   KN_ &operator++() {
      assert(next >= 0);
      v += next;
      return *this;
   } //    ++U
   KN_ &operator--() {
      assert(next >= 0);
      v -= next;
      return *this;
   } //    --U
   KN_ operator++(int) {
      assert(next >= 0);
      KN_ old = *this;
      v       = v + next;
      return old;
   } // U++
   KN_ operator--(int) {
      assert(next >= 0);
      KN_ old = *this;
      v       = v - next;
      return old;
   } // U++

   KN_(const KN_<R> &u, long offset) : ShapeOfArray(u), v(&u[offset]) {}
   KN_(const KN_<R> &u, const ShapeOfArray &sh, long startv = 0)
       : ShapeOfArray(sh * u.step), v(&u[startv]) {}
   KN_(const KN_<R> &u, long nnext, const ShapeOfArray &sh, long startv = 0)
       : ShapeOfArray(sh.n, sh.step * u.step, nnext), v(&u[startv]) {}

   friend class KNM_<R>;
   friend class KNMK_<R>;
   friend class KN<R>;
   friend class KNM<R>;
   friend class KNMK<R>;
};

template <class R> class KNM_ : public KN_<R> {
 public:
   ShapeOfArray shapei;
   ShapeOfArray shapej;

 public:
   long IsVector1() const { return (shapei.n * shapej.n) == this->n; }
   long N() const { return shapei.n; }
   long M() const { return shapej.n; }
   long size() const { return shapei.n * shapej.n; }

   KNM_(R *u, const ShapeOfArray &s, const ShapeOfArray &si,
        const ShapeOfArray &sj)
       : KN_<R>(u, s), shapei(si), shapej(sj) {}
   KNM_(R *u, long nn, long mm)
       : KN_<R>(u, ShapeOfArray(nn * mm)), shapei(nn, 1, nn),
         shapej(mm, nn, 1) {}
   KNM_(R *u, long nn, long mm, long s)
       : KN_<R>(u, ShapeOfArray(nn * mm, s)), shapei(nn, 1, nn),
         shapej(mm, nn, 1) {}
   KNM_(KN_<R> u, long n, long m)
       : KN_<R>(u, ShapeOfArray(m * n)), shapei(n, 1, n), shapej(m, n, 1) {}

   KNM_(const KN_<R> &u, const ShapeOfArray &si, const ShapeOfArray &sj,
        long offset = 0)
       : KN_<R>(&u[offset], si.last() + sj.last() + 1, u.step), shapei(si),
         shapej(sj) {
      assert(offset >= 0 && this->n + (this->v - (R *)u) <= u.n);
   }

   KNM_(const KN_<R> &u, const ShapeOfArray &si, const ShapeOfArray &sj,
        long offset, long nnext)
       : KN_<R>(&u[offset], si.last() + sj.last() + 1, u.step, nnext),
         shapei(si), shapej(sj) {
      assert(offset >= 0 && this->n + (this->v - (R *)u) <= u.n);
   }

   KNM_(const KNM_<R> &u) : KN_<R>(u), shapei(u.shapei), shapej(u.shapej) {}

   KNM_(KNM_<R> U, const SubArray &si, const SubArray &sj)
       : KN_<R>(U, SubArray(U.ij(si.len1(), sj.len1()) + 1,
                            U.ij(si.start, sj.start))),
         shapei(U.shapei, si), shapej(U.shapej, sj) {}

   KNM_(KNM_<R> U, const SubArray &sa, const SubArray &si, const SubArray &sj)
       : KN_<R>(U, SubArray(sa)), shapei(U.shapei, si), shapej(U.shapej, sj) {}

   KNM_ operator()(const SubArray &sa, const SubArray &sb) const {
      return KNM_(*this, sa, sb);
   }

   long ij(long i, long j) const { return shapei.index(i) + shapej.index(j); }
   long indexij(long i, long j) const {
      return this->index(shapei.index(i) + shapej.index(j));
   }
   R &operator()(long i, long j) const { return this->v[indexij(i, j)]; }
   R &operator()(int i, int j) const { return this->v[indexij(i, j)]; }

   KN_<R> operator()(const char, long j) const // une colonne j  ('.',j)
   {
      return KN_<R>(&this->v[this->index(shapej.index(j))],
                    shapei * this->step);
   }
   KN_<R> operator()(long i, const char) const // une ligne i  (i,'.')
   {
      return KN_<R>(&this->v[this->index(shapei.index(i))],
                    shapej * this->step);
   }
   KN_<R> operator()(const char, int j) const // une colonne j  ('.',j)
   {
      return KN_<R>(&this->v[this->index(shapej.index(j))],
                    shapei * this->step);
   }
   KN_<R> operator()(int i, const char) const // une ligne i  (i,'.')
   {
      return KN_<R>(&this->v[this->index(shapei.index(i))],
                    shapej * this->step);
   }
   KN_<R> operator()(const char, const char) const // tous
   {
      return *this;
   }
   KNM_<R> t() const { return KNM_<R>(this->v, *this, shapej, shapei); }

   KNM_ &operator=(const KNM_<R> &u);
   KNM_ &operator=(R a);
   KNM_ &operator+=(R a);
   KNM_ &operator-=(R a);
   KNM_ &operator/=(R a);
   KNM_ &operator*=(R a);
   KNM_ &operator+=(const KNM_<R> &u);
   KNM_ &operator-=(const KNM_<R> &u);
   KNM_ &operator*=(const KNM_<R> &u);
   KNM_ &operator/=(const KNM_<R> &u);

   friend class KN_<R>;
   friend class KNMK_<R>;
   friend class KN<R>;
   friend class KNM<R>;
   friend class KNMK<R>;
};

template <class R> class KNMK_ : public KN_<R> {
   friend class KNMK<R>;

 public:
   ShapeOfArray shapei;
   ShapeOfArray shapej;
   ShapeOfArray shapek;

 public:
   long IsVector1() const {
      return (shapei.n * shapej.n * shapek.n) == this->n;
   }
   long N() const { return shapei.n; }
   long M() const { return shapej.n; }
   long K() const { return shapek.n; }
   long size() const { return shapei.n * shapej.n * shapek.n; }
   KNMK_(const ShapeOfArray &s, const ShapeOfArray &si, const ShapeOfArray &sj,
         const ShapeOfArray &sk, R *u)
       : KN_<R>(u, s), shapei(si), shapej(sj), shapek(sk) {}

   KNMK_(R *u, long n, long m, long k)
       : KN_<R>(u, ShapeOfArray(n * m * k)), shapei(n, 1, n), shapej(m, n, 1),
         shapek(k, n * m, n * m){};

   //  KNMK_(const KN_<R> & u,long n,long m,long k)
   //   : KN_<R>(ShapeOfArray(n*m*k)),shapei(n,1,n),shapekj(m,n,1),u),
   //     shapek(k,n*m,n*m){};

   KNMK_(const KNMK_<R> &U, const SubArray &si, const SubArray &sj,
         const SubArray &sk)
       : KN_<R>(U, SubArray(U.ijk(si.len1(), sj.len1(), sk.len1()) + 1,
                            U.ijk(si.start, sj.start, sk.start))),
         shapei(U.shapei, si), shapej(U.shapej, sj), shapek(U.shapek, sk) {}

   KNMK_(const KNMK_<R> &u)
       : KN_<R>(u), shapei(u.shapei), shapej(u.shapej), shapek(u.shapek) {}

   long ijk(long i, long j, long k) const {
      return shapei.index(i) + shapej.index(j) + shapek.index(k);
   }
   long indexijk(long i, long j, long k) const {
      return this->index(shapei.index(i) + shapej.index(j) + shapek.index(k));
   }

   R &operator()(long i, long j, long k) const {
      return this->v[indexijk(i, j, k)];
   }
   R &operator()(int i, int j, int k) const {
      return this->v[indexijk(i, j, k)];
   }

   //  pas de tableau suivant
   KN_<R> operator()(const char, long j, long k) const { // le tableau (.,j,k)
      return KN_<R>(*this, -1, shapei, shapej[j] + shapek[k]);
   }
   // KN_<R>  operator()(long i,const char ,long k)  const  { // le tableau
   // (i,.,k)
   //        return KN_<R>(*this,-1,shapej,shapei[i]+shapek[k]);}
   // KN_<R>  operator()(long i,long j,const char )  const  { // le tableau
   // (i,j,.)
   //        return KN_<R>(*this,-1,shapek,shapei[i]+shapej[j]);}
   //
   KN_<R> operator()(const char, int j, int k) const { // le tableau (.,j,k)
      return KN_<R>(*this, -1, shapei, shapej[j] + shapek[k]);
   }
   // KN_<R>  operator()(int i,const char ,int k)  const  { // le tableau
   // (i,.,k)
   //        return KN_<R>(*this,-1,shapej,shapei[i]+shapek[k]);}
   // KN_<R>  operator()(int i,int j,const char )  const  { // le tableau
   // (i,j,.)
   //        return KN_<R>(*this,-1,shapek,shapei[i]+shapej[j]);}
   //
   KNM_<R> operator()(const char, const char,
                      long k) const { // le tableau (.,.,k)
      return KNM_<R>(*this, shapei, shapej, shapek[k], shapek.next);
   } // step = n*m
   // //attention les suivants ne marche pas
   // KNM_<R>  operator()(const char ,long j,const char )  const  { // le
   // tableau (.,j,.)
   //        return KNM_<R>(*this,shapei,shapek,shapej[j],-1/*shapej.next*/);}
   //        // step = n
   //
   // KNM_<R>  operator()(long i,const char ,const char )  const  { // le
   // tableau (i,.,.)
   //        return KNM_<R>(*this,shapej,shapek,shapei[i],-1/*shapei.next*/);}
   //        // step = 1
   //
   KNM_<R> operator()(const char, const char,
                      int k) const { // le tableau (.,.,k)
      return KNM_<R>(*this, shapei, shapej, shapek[k], shapek.next);
   } // step = n*m
     // //attention les suivants ne marche pas
   // KNM_<R>  operator()(const char ,int j,const char )  const  { // le tableau
   // (.,j,.)
   //        return KNM_<R>(*this,shapei,shapek,shapej[j],-1/*shapej.next*/);}
   //        // step = n
   //
   // KNM_<R>  operator()(int i,const char ,const char )  const  { // le tableau
   // (i,.,.)
   //        return KNM_<R>(*this,shapej,shapek,shapei[i],-1/*shapei.next*/);}
   //        // step = 1

   KNMK_ &operator=(const KNMK_<R> &u);
   KNMK_ &operator+=(const KNMK_<R> &u);
   KNMK_ &operator-=(const KNMK_<R> &u);
   KNMK_ &operator/=(const KNMK_<R> &u);
   KNMK_ &operator*=(const KNMK_<R> &u);
   KNMK_ &operator=(R a);
   KNMK_ &operator+=(R a);
   KNMK_ &operator-=(R a);
   KNMK_ &operator/=(R a);
   KNMK_ &operator*=(R a);

   KNMK_ operator()(SubArray si, SubArray sj, SubArray sk) const {
      return KNMK_(*this, si, sj, sk);
   }

 private:
   friend class KNM_<R>;
   friend class KN_<R>;
};

//------------------------------------

template <class R> class KN : public KN_<R> {
 public:
   typedef R K;

   KN() : KN_<R>(0, 0) {}
   KN(long nn) : KN_<R>(new R[nn], nn) {}
   KN(long nn, R *p) : KN_<R>(new R[nn], nn) {
      KN_<R>::operator=(KN_<R>(p, nn));
   }
   KN(long nn, R (*f)(long i)) : KN_<R>(new R[nn], nn) {
      for (long i = 0; i < this->n; i++)
         this->v[i] = f(i);
   }
   KN(long nn, const R &a) : KN_<R>(new R[nn], nn) { KN_<R>::operator=(a); }
   KN(long nn, long s, const R a) : KN_<R>(new R[nn], nn, s) {
      KN_<R>::operator=(a);
   }
   template <class S> KN(const KN_<S> &s) : KN_<R>(new R[s.n], s.n) {
      for (long i = 0; i < this->n; i++)
         this->v[i] = s[i];
   }
   template <class S> KN(const KN_<S> &s, R (*f)(S)) : KN_<R>(new R[s.n], s.n) {
      for (long i = 0; i < this->n; i++)
         this->v[i] = f(s[i]);
   }
   KN(const KN<R> &u) : KN_<R>(new R[u.n], u.n) { KN_<R>::operator=(u); }
   KN(bool, KN<R> &u) : KN_<R>(u) {
      u.v = 0;
      u.n = 0;
   } // remove copy for return of local KN.

   KN(std::list<R> &u) : KN_<R>(new R[u.size()], u.size()) {
      int ii = 0;
      for (auto it = u.begin(); it != u.end(); ++it)
         this->v[ii++] = *it;
   }

   ~KN() { delete[] this->v; }

   void CheckSet() {
      if (!(this->n)) {
         std::cerr << "Error RNM set array\n";
         assert(0);
         exit(1);
      }
   }
   KN &operator=(R *a) {
      CheckSet();
      return operator=(KN_<R>(a, this->n));
   }
   KN &operator+=(R *a) {
      CheckSet();
      return operator+=(KN_<R>(a, this->n));
   }
   KN &operator-=(R *a) {
      CheckSet();
      return operator-=(KN_<R>(a, this->n));
   }
   KN &operator*=(R *a) {
      CheckSet();
      return operator*=(KN_<R>(a, this->n));
   }
   KN &operator/=(R *a) {
      CheckSet();
      return operator/=(KN_<R>(a, this->n));
   }

   KN &operator=(R a) {
      if (this->unset())
         this->set(new R[1], 1, 0, 0);
      KN_<R>::operator=(a);
      return *this;
   }
   KN &operator=(const KN_<R> &a) {
      if (this->unset())
         this->set(new R[a.N()], a.N());
      KN_<R>::operator=(a);
      return *this;
   }
   KN &operator=(const KN<R> &a) {
      if (this->unset())
         this->set(new R[a.N()], a.N());
      KN_<R>::operator=(a);
      return *this;
   }
   KN &operator=(const Add_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.a.N()], u.a.N());
      KN_<R>::operator=(u);
      return *this;
   }
   KN &operator=(const Sub_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.a.N()], u.a.N());
      KN_<R>::operator=(u);
      return *this;
   }
   KN &operator=(const Mulc_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.a.N()], u.a.N());
      KN_<R>::operator=(u);
      return *this;
   }
   KN &operator=(const Add_Mulc_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.a.N()], u.a.N());
      KN_<R>::operator=(u);
      return *this;
   }
   KN &operator=(const Mul_KNM_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.b.N()], u.b.N());
      KN_<R>::operator=(u);
      return *this;
   }
   KN &operator=(const DotMul_KN_<R> &u) {
      if (this->unset())
         this->set(new R[u.a.N()], u.a.N());
      KN_<R>::operator=(u);
      return *this;
   }

   KN &operator-=(R a) {
      KN_<R>::operator-=(a);
      return *this;
   }
   KN &operator-=(const KN_<R> &a) {
      KN_<R>::operator-=(a);
      return *this;
   }
   KN &operator-=(const Add_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }
   KN &operator-=(const Sub_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }
   KN &operator-=(const DotMul_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }
   KN &operator-=(const Mulc_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }
   KN &operator-=(const Add_Mulc_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }
   KN &operator-=(const Mul_KNM_KN_<R> &u) {
      KN_<R>::operator-=(u);
      return *this;
   }

   KN &operator+=(R a) {
      KN_<R>::operator+=(a);
      return *this;
   }
   KN &operator+=(const KN_<R> &a) {
      KN_<R>::operator+=(a);
      return *this;
   }
   KN &operator+=(const Add_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }
   KN &operator+=(const Sub_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }
   KN &operator+=(const Mulc_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }
   KN &operator+=(const Add_Mulc_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }
   KN &operator+=(const Mul_KNM_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }
   KN &operator+=(const DotMul_KN_<R> &u) {
      KN_<R>::operator+=(u);
      return *this;
   }

   KN &operator/=(R a) {
      KN_<R>::operator/=(a);
      return *this;
   }
   KN &operator/=(const KN_<R> &a) {
      KN_<R>::operator/=(a);
      return *this;
   }
   KN &operator/=(const Add_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }
   KN &operator/=(const Sub_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }
   KN &operator/=(const Mulc_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }
   KN &operator/=(const Add_Mulc_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }
   KN &operator/=(const Mul_KNM_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }
   KN &operator/=(const DotMul_KN_<R> &u) {
      KN_<R>::operator/=(u);
      return *this;
   }

   KN &operator*=(R a) {
      KN_<R>::operator*=(a);
      return *this;
   }
   KN &operator*=(const KN_<R> &a) {
      KN_<R>::operator*=(a);
      return *this;
   }
   KN &operator*=(const Add_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }
   KN &operator*=(const Sub_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }
   KN &operator*=(const Mulc_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }
   KN &operator*=(const Add_Mulc_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }
   KN &operator*=(const Mul_KNM_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }
   KN &operator*=(const DotMul_KN_<R> &u) {
      KN_<R>::operator*=(u);
      return *this;
   }

   static void fill0(R *v, int n) {
      if (n && v)
         for (int i = 0; i < n; ++i)
            v[i] = R();
   }
   void init(long nn) {
      this->n    = nn;
      this->step = 1;
      this->next = -1;
      if (this->v)
         delete[] this->v;
      this->v = new R[nn];
      fill0(this->v, this->n);
   }
   void init() {
      this->n    = 0;
      this->step = 1;
      this->next = -1;
      this->v    = 0;
   }
   void init(const KN_<R> &a) {
      init(a.N());
      operator=(a);
   }
   void resize(long nn) {
      if (nn != this->n) {
         R *vo   = this->v;
         long no = std::min(this->n, nn), so = this->step;
         ShapeOfArray::init(nn);
         this->v = new R[this->n];
         // copy
         if (this->v && vo)
            for (long i = 0, j = 0; j < no; i++, j += so)
               this->v[i] = vo[j];
         delete[] vo;
      }
   } //  mars 2010
   void destroy() {
      assert((this->next) < 0);
      if (this->next++ == -1) {
         delete[] this->v;
         this->v = 0;
         this->n = 0;
      }
   } //  mars 2010
   void increment() {
      assert((this->next) < 0);
      this->next--;
   }
   void swap(KN_<R> &a) {
      assert(a.n == this->n);
      R *vo   = this->v;
      this->v = a.v;
      a.v     = vo;
   }
   void load(const std::string &filename) {
      std::ifstream f(filename.c_str(), std::ifstream::in);
      int lend;
      f >> lend;
      init(lend);
      for (int i = 0; i < lend; ++i)
         f >> this->v[i];
      f.close();
   }
};

//  Array with 2 indices
//  ---------------------
template <class R> class KNM : public KNM_<R> {
 public:
   KNM() : KNM_<R>(0, 0, 0) {}
   KNM(long nn, long mm) : KNM_<R>(new R[nn * mm], nn, mm) {}
   KNM(long nn, long mm, R *p) : KNM_<R>(new R[nn], nn, mm) {
      KNM_<R>::operator=(KNM_<R>(p, nn, mm));
   }
   KNM(const KNM<R> &u) : KNM_<R>(new R[u.size()], u.N(), u.M()) {
      KN_<R>::operator=(u);
   }
   explicit KNM(const KNM_<R> &u) : KNM_<R>(new R[u.size()], u.N(), u.M()) {
      KNM_<R>::operator=(u);
   }
   ~KNM() { delete[] this->v; }

   KNM &operator=(const KNM_<R> &u) {
      if (this->unset())
         this->init(u.N(), u.M());
      KNM_<R>::operator=(u);
      return *this;
   }

   KNM &operator=(R a) {
      if (this->unset())
         assert(0 && " KNM operator=(double)");
      KNM_<R>::operator=(a);
      return *this;
   }
   KNM &operator+=(R a) {
      if (this->unset())
         assert(0 && " KNM operator+=(double)");
      KNM_<R>::operator+=(a);
      return *this;
   }
   KNM &operator-=(R a) {
      if (this->unset())
         assert(0 && " KNM operator-=(double)");
      KNM_<R>::operator-=(a);
      return *this;
   }
   KNM &operator/=(R a) {
      if (this->unset())
         assert(0 && " KNM operator/=(double)");
      KNM_<R>::operator/=(a);
      return *this;
   }
   KNM &operator*=(R a) {
      if (this->unset())
         assert(0 && " KNM operator*=(double)");
      KNM_<R>::operator*=(a);
      return *this;
   }
   KNM &operator+=(const KNM_<R> &u) {
      if (this->unset())
         this->init(u.N(), u.M());
      KNM_<R>::operator+=(u);
      return *this;
   }
   KNM &operator-=(const KNM_<R> &u) {
      if (this->unset())
         this->init(u.N(), u.M());
      KNM_<R>::operator-=(u);
      return *this;
   }

   KNM &operator/=(const KNM_<R> &u) {
      if (this->unset())
         this->init(u.N(), u.M());
      KNM_<R>::operator/=(u);
      return *this;
   }
   KNM &operator*=(const KNM_<R> &u) {
      if (this->unset())
         this->init(u.N(), u.M());
      KNM_<R>::operator*=(u);
      return *this;
   }

   void init() { //  add mars 2010 ...
      this->n    = 0;
      this->step = 1;
      this->next = -1;
      this->v    = 0;
      this->shapei.init(0);
      this->shapej.init(0);
   }

   void init(long nn, long mm) {
      ShapeOfArray::init(nn * mm);
      this->shapei.init(nn, 1, nn);
      this->shapej.init(mm, nn, 1);
      this->v = new R[nn * mm];
   }
   void destroy() {
      assert((this->next) < 0);
      if (this->next++ == -1) {
         delete[] this->v;
         this->v = 0;
         this->n = 0;
      }
   }
   void increment() {
      assert((this->next) < 0);
      this->next--;
   }
};

//  Array with 3 indices
//  ---------------------
template <class R> class KNMK : public KNMK_<R> {
 public:
   KNMK() : KNMK_<R>(0, 0, 0, 0) {}
   KNMK(long n, long m, long k) : KNMK_<R>(new R[n * m * k], n, m, k) {}
   explicit KNMK(const KNMK_<R> &u)
       : KNMK_<R>(new R[u.size()], u.N(), u.M(), u.K()) {
      KNMK_<R>::operator=(u);
   }
   KNMK(const KNMK<R> &u) : KNMK_<R>(new R[u.size()], u.N(), u.M(), u.K()) {
      KNMK_<R>::operator=(u);
   }

   ~KNMK() { delete[] this->v; }

   KNMK &operator=(const KNMK_<R> &u) {
      KNMK_<R>::operator=(u);
      return *this;
   }
   KNMK &operator=(R a) {
      KNMK_<R>::operator=(a);
      return *this;
   }
   KNMK &operator+=(R a) {
      KNMK_<R>::operator+=(a);
      return *this;
   }
   KNMK &operator-=(R a) {
      KNMK_<R>::operator-=(a);
      return *this;
   }
   KNMK &operator/=(R a) {
      KNMK_<R>::operator/=(a);
      return *this;
   }
   KNMK &operator*=(R a) {
      KNMK_<R>::operator*=(a);
      return *this;
   }
   KNMK &operator+=(const KNMK_<R> &u) {
      KNMK_<R>::operator+=(u);
      return *this;
   }
   // ici jd
   KNMK &operator-=(const KNMK_<R> &u) {
      KNMK_<R>::operator-=(u);
      return *this;
   }
   KNMK &operator*=(const KNMK_<R> &u) {
      KNMK_<R>::operator*=(u);
      return *this;
   }
   KNMK &operator/=(const KNMK_<R> &u) {
      KNMK_<R>::operator/=(u);
      return *this;
   }

   void init(long nn, long mm, long kk) {
      ShapeOfArray::init(nn * mm * kk);
      this->shapei.init(nn, 1, nn);
      this->shapej.init(mm, nn, 1);
      this->shapek.init(kk, nn * mm, nn * mm);
      this->v = new R[nn * mm * kk];
   }
};

//  -------------  optimization ---------------------
template <class R> class Add_KN_ {
 public:
   const KN_<R> a;
   const KN_<R> b;
   Add_KN_(const KN_<R> &aa, const KN_<R> &bb) : a(aa), b(bb) {
      assert(SameShape(a, b));
   }
};

template <class R> class Sub_KN_ {
 public:
   const KN_<R> a;
   const KN_<R> b;
   Sub_KN_(const KN_<R> &aa, const KN_<R> &bb) : a(aa), b(bb) {
      assert(SameShape(a, b));
   }
};

template <class R> class DotMul_KN_ {
 public:
   const KN_<R> a;
   const KN_<R> b;
   DotMul_KN_(const KN_<R> &aa, const KN_<R> &bb) : a(aa), b(bb) {
      assert(SameShape(a, b));
   }
};

template <class R> class Mulc_KN_ {
 public:
   const KN_<R> a;
   R b;
   Mulc_KN_(const KN_<R> &aa, R bb) : a(aa), b(bb) {}
   Mulc_KN_(const Mulc_KN_<R> &aa, R bb) : a(aa.a), b(aa.b * bb) {}
   Mulc_KN_ operator-() const { return Mulc_KN_(a, -b); }
};

template <class R> class Add_Mulc_KN_ {
 public:
   const KN_<R> a, b;
   const R ca, cb;
   Add_Mulc_KN_(const Mulc_KN_<R> &aa, const Mulc_KN_<R> &bb)
       : a(aa.a), b(bb.a), ca(aa.b), cb(bb.b) {
      assert(SameShape(a, b));
   }
   Add_Mulc_KN_(const Mulc_KN_<R> &aa, const KN_<R> &bb, const R cbb)
       : a(aa.a), b(bb), ca(aa.b), cb(cbb) {
      assert(SameShape(a, b));
   }
   Add_Mulc_KN_(const KN_<R> &aa, const R caa, const KN_<R> &bb, const R cbb)
       : a(aa), b(bb), ca(caa), cb(cbb) {
      assert(SameShape(a, b));
   }
};

template <class R> class Mul_KNM_KN_ {
 public:
   const KNM_<R> &A;
   const KN_<R> &b;
   Mul_KNM_KN_(const KNM_<R> &aa, const KN_<R> &bb) : A(aa), b(bb) {
      assert(SameShape(A.shapej, b));
   }
};

std::ostream &operator<<(std::ostream &f, const ShapeOfArray &s);

template <class R> std::ostream &operator<<(std::ostream &f, const KN_<R> &v);
template <class R> std::ostream &operator<<(std::ostream &f, const KNM_<R> &v);
template <class R>
inline std::ostream &operator<<(std::ostream &f, const KN<R> &v) {
   return f << (const KN_<R> &)v;
}
template <class R>
inline std::ostream &operator<<(std::ostream &f, const KNM<R> &v) {
   return f << (const KNM_<R> &)v;
}

template <class R>
inline Add_KN_<R> operator+(const KN_<R> &a, const KN_<R> &b) {
   return Add_KN_<R>(a, b);
}
template <class R>
inline Sub_KN_<R> operator-(const KN_<R> &a, const KN_<R> &b) {
   return Sub_KN_<R>(a, b);
}
template <class R> inline Mulc_KN_<R> operator*(const KN_<R> &a, const R &b) {
   return Mulc_KN_<R>(a, b);
}
template <class R> inline Mulc_KN_<R> operator/(const KN_<R> &a, const R &b) {
   return Mulc_KN_<R>(a, 1 / b);
}
template <class R> inline Mulc_KN_<R> operator*(const R &b, const KN_<R> &a) {
   return Mulc_KN_<R>(a, b);
}
template <class R> inline Mulc_KN_<R> operator-(const KN_<R> &a) {
   return Mulc_KN_<R>(a, -1);
}
template <class R>
inline DotMul_KN_<R> operator*(const KN_<R> &a, const KN_<R> &b) {
   return DotMul_KN_<R>(a, b);
}

template <class R>
inline Add_Mulc_KN_<R> operator+(const Mulc_KN_<R> &a, const Mulc_KN_<R> &b) {
   return Add_Mulc_KN_<R>(a, b);
}
template <class R>
inline Add_Mulc_KN_<R> operator-(const Mulc_KN_<R> &a, const Mulc_KN_<R> &b) {
   return Add_Mulc_KN_<R>(a, b.a, -b.b);
}

template <class R>
inline Add_Mulc_KN_<R> operator+(const Mulc_KN_<R> &a, const KN_<R> &b) {
   return Add_Mulc_KN_<R>(a, b, R(1));
}
template <class R>
inline Add_Mulc_KN_<R> operator-(const Mulc_KN_<R> &a, const KN_<R> &b) {
   return Add_Mulc_KN_<R>(a, b, R(-1));
}

template <class R>
inline Add_Mulc_KN_<R> operator+(const KN_<R> &b, const Mulc_KN_<R> &a) {
   return Add_Mulc_KN_<R>(a, b, R(1));
}

// modif FH mars 2007
template <class R>
inline Add_Mulc_KN_<R> operator-(const KN_<R> &a, const Mulc_KN_<R> &b) {
   return Add_Mulc_KN_<R>(a, R(1), b.a, -b.b);
} // modif FH mars 2007

template <class R>
inline Mul_KNM_KN_<R> operator*(const KNM_<R> &A, const KN_<R> &b) {
   return Mul_KNM_KN_<R>(A, b);
}

inline bool SameShape(const ShapeOfArray &a, const ShapeOfArray &b) {
   return a.SameShape(b);
}
template <class R>
inline bool SameShape(const ShapeOfArray &a, const Add_Mulc_KN_<R> &b) {
   return SameShape(a, b.a);
}
template <class R>
inline bool SameShape(const ShapeOfArray &a, const Add_KN_<R> &b) {
   return SameShape(a, b.a);
}
template <class R>
inline bool SameShape(const ShapeOfArray &a, const Sub_KN_<R> &b) {
   return SameShape(a, b.a);
}
template <class R>
inline bool SameShape(const ShapeOfArray &a, const Mulc_KN_<R> &b) {
   return SameShape(a, b.a);
}
template <class R>
inline bool SameShape(const ShapeOfArray &a, const Mul_KNM_KN_<R> &b) {
   return a.n == b.A.N();
}
inline bool SameShape(const ShapeOfArray &, const double) { return true; }

template <class R> inline long SameAdress(const KN_<R> &a, const KN_<R> &b) {
   return &a[0] == &b[0];
}

#include "RNM_tpl.hpp"
#ifdef assert
#undef assert
#endif
#endif
