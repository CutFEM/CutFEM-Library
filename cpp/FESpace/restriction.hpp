#ifndef RESTRICTION_HPP_
#define RESTRICTION_HPP_

#include "../common/RNM.hpp"
#include "FESpace.hpp"

template <typename M>
void restriction(const M &Vh, const KN<double> &fh, const M &cutVh,
                 KN<double> &g) {

   g.init(cutVh.NbDoF());
   g = 0.0;

   for (int k = 0; k < cutVh.NbElement(); ++k) {
      const int kb = cutVh.idxElementInBackMesh(k);
      const int kg = Vh.idxElementFromBackMesh(kb);

      const typename M::FElement &FK(cutVh[k]);
      const typename M::FElement &FKg(Vh[kg]);
      assert(FK.NbDoF() == FKg.NbDoF());
      for (int ic = 0; ic < Vh.N; ++ic) {
         for (int i = FK.dfcbegin(ic); i < FK.dfcend(ic); ++i) {
            g(FK(i)) = fh(FKg(i));
         }
      }
   }
}

#endif
