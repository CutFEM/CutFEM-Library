

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
