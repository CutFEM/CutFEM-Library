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
