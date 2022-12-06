
/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

template <class R> KN_<R> &KN_<R>::operator oper(const Mul_KNM_KN_<R> &u) {
   assert(SameShape(u.A.shapei) && !constant());
   R *l(v);
   KN_<R> li(u.A(0, '.')); //  first line
   for (long i = 0; i < n; i++, l += step, ++li)
      *l oper(li, u.b);
   return *this;
}

template <class R> KN_<R> &KN_<R>::operator oper(const Add_KN_<R> &u) {
   assert(u.a.N() == N());
   long stepa(u.a.step), stepb(u.b.step);
   R *l(v);
   R *aa(u.a), *bb(u.b);
   for (long i = 0; i < n; i++, l += step, aa += stepa, bb += stepb)
      *l oper *aa + *bb;
   return *this;
}

template <class R> KN_<R> &KN_<R>::operator oper(const DotMul_KN_<R> &u) {
   assert(u.a.N() == N());
   long stepa(u.a.step), stepb(u.b.step);
   R *l(v);
   R *aa(u.a), *bb(u.b);
   for (long i = 0; i < n; i++, l += step, aa += stepa, bb += stepb)
      *l oper(*aa) * (*bb);
   return *this;
}

template <class R> KN_<R> &KN_<R>::operator oper(const Sub_KN_<R> &u) {
   assert(u.a.N() == N());
   long stepa(u.a.step), stepb(u.b.step);
   R *l(v);
   R *aa(u.a), *bb(u.b);
   for (long i = 0; i < n; i++, l += step, aa += stepa, bb += stepb)
      *l oper *aa - *bb;
   return *this;
}

template <class R> KN_<R> &KN_<R>::operator oper(const Mulc_KN_<R> &u) {
   assert(u.a.N() == N());
   long stepa(u.a.step);
   R *l(v);
   R *aa(u.a), bb(u.b);
   for (long i = 0; i < n; i++, l += step, aa += stepa)
      *l oper *aa *bb;
   return *this;
}

template <class R> KN_<R> &KN_<R>::operator oper(const Add_Mulc_KN_<R> &u) {
   assert(u.a.N() == N());
   const long stepa(u.a.step), stepb(u.b.step);
   const R ca(u.ca), cb(u.cb);
   R *l(v);
   const R *aa(u.a), *bb(u.b);
   for (long i = 0; i < n; i++, l += step, aa += stepa, bb += stepb)
      *l oper *aa *ca + *bb *cb;
   return *this;
}
