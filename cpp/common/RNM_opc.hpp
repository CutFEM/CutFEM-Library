

template <class R> inline KN_<R> &KN_<R>::operator oper(R a) {
   R *l(v);
   for (long i = 0; i < n; i++, l += step)
      *l oper a;
   return *this;
}

template <class R> inline KNM_<R> &KNM_<R>::operator oper(R a) {
   if (IsVector1())
      KN_<R>::operator oper(a);
   else {
      KN_<R> lj(operator()('.', 0)); //  (.,.,O)
      for (long j = 0; j < M(); ++j, ++lj)
         lj oper a;
   }
   return *this;
}

template <class R> inline KNMK_<R> &KNMK_<R>::operator oper(R a) {
   if (IsVector1())
      KN_<R>::operator oper(a);
   else {
      KNM_<R> lj(operator()('.', '.', 0)); //  (.,.,O)
      long j = K();
      while (j--) {
         lj oper a;
         ++lj;
      }
   }
   return *this;
}

template <class R> inline KN_<R> &KN_<R>::operator oper(const KN_<R> &u) {
   assert(u.n == n);
   R *l(v);
   const R *r(u);
   for (long i = 0; i < n; i++, l += step, r += u.step)
      *l oper *r;
   return *this;
}

template <class R> inline KNM_<R> &KNM_<R>::operator oper(const KNM_<R> &u) {
   assert(N() == u.N() && M() == u.M());
   if (IsVector1() && u.IsVector1() &&
       shapei.step == u.shapei.step) // modif 2011 (thank to Oka)
      KN_<R>::operator oper(u);      // modif FH jan 2004
   else {
      KN_<R> lj(operator()('.', 0)); //  (.,O)
      KN_<R> uj(u('.', 0));
      long j = M();
      while (j--) {
         lj oper uj;
         ++lj;
         ++uj;
      }
   }
   return *this;
}

template <class R> inline KNMK_<R> &KNMK_<R>::operator oper(const KNMK_<R> &u) {
   assert(N() == u.N() && M() == u.M() && K() == u.K());

   if (IsVector1() && u.IsVector1() && u.N() == N() &&
       shapei.step == u.shapei.step) // modif 2011 (thank to Oka)
      KN_<R>::operator oper(u);      // modif FH 2004
   else {
      assert(K() == u.K());
      KNM_<R> lj(operator()('.', '.', 0)); //  (.,O)
      KNM_<R> uj(u('.', '.', 0));
      long j = K();
      while (j--) {
         lj oper uj;
         ++lj;
         ++uj;
      }
   }
   return *this;
}

#undef oper
