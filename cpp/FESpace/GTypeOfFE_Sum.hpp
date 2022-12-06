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

#ifndef _GTYPE_OF_FE_SUM_HPP
#define _GTYPE_OF_FE_SUM_HPP

#include "GTypeOfFE.hpp"
#include "../common/Mesh3dn.hpp"
#include "../common/Mesh2dn.hpp"
#include "../common/Mesh1dn.hpp"

template <class Mesh> class GTypeOfFESum : public GTypeOfFE<Mesh> {

 public:
   typedef GFElement<Mesh> FElement;
   typedef typename Mesh::Element Element;
   typedef typename Element::RdHat RdHat;
   typedef typename Mesh::Rd Rd;

   const int k;
   KN<const GTypeOfFE<Mesh> *> teb;
   KN<int> NN, DF,
       comp; // increasing array of N, dof and comp (to not recompute BF)

   GTypeOfFESum(const KN<GTypeOfFE<Mesh> const *> &t);

   void init(InterpolationMatrix<RdHat> &M) const;
   void FB(const What_d whatd, const Element &K, const Rd &P,
           KNMK_<R> &val) const;
   void FB(const What_d whatd, const Element &K, const Rd &P, KNMK_<R> &val,
           const RNM_ &J) const;
   virtual void get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                              KN_<double> &v) const;
   void build();
   ~GTypeOfFESum() {}
};

template <class Mesh>
GTypeOfFESum<Mesh>::GTypeOfFESum(const KN<GTypeOfFE<Mesh> const *> &t)
    : GTypeOfFE<Mesh>(t), k(t.N()), teb(t), NN(k + 1), DF(k + 1), comp(k) {
   build();

   this->polynomialOrder = 0;
   for (int i = 0; i < k; ++i) {
      this->polynomialOrder =
          std::max(this->polynomialOrder, t[i]->polynomialOrder);
   }
   this->basisFctType = t[0]->basisFctType;
   for (int i = 1; i < k; ++i) {
      if (t[i]->basisFctType != this->basisFctType)
         this->basisFctType = BasisFctType::UNDEFINED;
   }
}

template <class Mesh> void GTypeOfFESum<Mesh>::build() {
   std::map<const GTypeOfFE<Mesh> *, int> m;
   int i = k, j;
   while (i--)       // backward to get comp[i] <=i
      m[teb[i]] = i; // saving adress of FE

   i = k;
   while (i--)
      comp[i] = m[teb[i]]; //  take the value of smallest comp to not
   // recompute BF comp[i] <=i
   // reserve space intervals
   int n = 0, N = 0;
   for (j = 0; j < k; j++) {
      NN[j] = N;
      N += teb[j]->N;
   }
   NN[k] = N;
   //  reserve DF intervals
   n     = 0;
   for (j = 0; j < k; j++) {
      DF[j] = n;
      n += teb[j]->nbDoF;
   }
   DF[k] = n;

   // access ALL the subFE
   int ii = 0;
   for (int i = 0; i < k; ++i) {
      for (int j = 0; j < teb[i]->nbOfFE; ++j) {
         this->Sub_ToFE[ii++] = teb[i]->Sub_ToFE[j];
      }
   }

   // get the df of each componant
   int c = 0, c0 = 0;
   for (int i = 0; i < this->nbOfFE; i++) {
      int N     = this->Sub_ToFE[i]->N;
      int ndofi = this->Sub_ToFE[i]->nbDoF;
      for (int j = 0; j < N; ++j) {
         this->begin_dfcomp[c] = c0 + this->Sub_ToFE[i]->begin_dfcomp[j];
         this->end_dfcomp[c]   = c0 + this->Sub_ToFE[i]->end_dfcomp[j];
         c++;
      }
      c0 += ndofi;
   }

   // construction de l'interpolation .
   int npi = 0;
   int nki = 0;
   for (int i = 0; i < this->nbOfFE; i++) {
      npi += this->Sub_ToFE[i]->NbPtforInterpolation;
      nki += this->Sub_ToFE[i]->NbcoefforInterpolation;
   }

   KN<int> numPt(npi);
   {
      std::map<RdHat, int, lessRd> mpt;
      int npp = 0, kkk = 0;
      KN<RdHat> Ptt(npi);
      for (int i = 0; i < this->nbOfFE; i++) {
         const GTypeOfFE<Mesh> &ti = *this->Sub_ToFE[i];

         for (int p = 0; p < ti.NbPtforInterpolation; ++p, ++kkk) {
            Ptt[kkk] = ti.Pt_Pi_h[p];

            if (mpt.find(Ptt[kkk]) == mpt.end())
               mpt[Ptt[kkk]] = npp++;
            numPt[kkk] = mpt[Ptt[kkk]];
         }
      }

      this->NbPtforInterpolation = npp;
      this->Pt_Pi_h.init(npp);
      for (int i = 0; i < npp; ++i)
         this->Pt_Pi_h[numPt[i]] = Ptt[i];

      // get the coefficient
      this->ipj_Pi_h.init(nki);
      int lll = 0;
      int ii0 = 0, jj0 = 0, pp0 = 0;
      for (int i = 0; i < this->nbOfFE; i++) {
         const GTypeOfFE<Mesh> &ti = *this->Sub_ToFE[i];
         const KN<IPJ> p(ti.ipj_Pi_h);
         for (int j = 0; j < ti.NbcoefforInterpolation; ++j, ++lll) {
            this->ipj_Pi_h[lll].i = p[j].i + ii0;
            this->ipj_Pi_h[lll].p = numPt[p[j].p + pp0];
            this->ipj_Pi_h[lll].j = p[j].j + jj0;
         }
         ii0 += ti.nbDoF;
         jj0 += ti.N;
         pp0 += ti.NbPtforInterpolation;
      }

      // for(int i=0;i<this->ipj_Pi_h.N();++i) {
      //   std::cout << this->ipj_Pi_h[i].i << "\t"
      // 		<< this->ipj_Pi_h[i].p << "\t"
      // 		<< this->ipj_Pi_h[i].j << "\n";
      // }
      // getchar();
   }
}

template <class Mesh>
void GTypeOfFESum<Mesh>::FB(const What_d whatd, const Element &K, const Rd &P,
                            KNMK_<R> &val) const {
   val = 0.0;
   SubArray t(val.K());

   for (int i = 0; i < k; i++) {
      int j   = comp[i];
      int ni  = NN[i];
      int di  = DF[i];
      int i1  = i + 1;
      int nii = NN[i1];
      int dii = DF[i1];

      assert(ni < nii && di < dii);
      RNMK_ v(val(SubArray(dii - di, di), SubArray(nii - ni, ni), t));
      if (i <= j) {
         teb[i]->FB(whatd, K, P, v);
      } else {
         v = val(SubArray(DF[j + 1] - DF[j], DF[j]),
                 SubArray(NN[j + 1] - NN[j], NN[j]), t);
      }
   }
}

static KNM<R> mappDx(const KNM_<R> &A, const KNM_<R> &w) {
   int op[3] = {op_dx, op_dy, op_dz};
   KNM<R> x(w.N(), w.M());
   x = w; // phi(x) does not change
   for (int e = 0; e < w.N(); ++e) {
      for (int i = 0; i < A.N(); ++i) {
         x(e, op[i]) = 0;
         for (int j = 0; j < A.M(); ++j) {
            x(e, op[i]) += A(i, j) * w(e, op[j]);
         }
      }
   }
   return x;
}
template <class Mesh>
void GTypeOfFESum<Mesh>::FB(const What_d whatd, const Element &K, const Rd &P,
                            KNMK_<R> &val, const RNM_ &J) const {
   assert(0);
   // val=0.0;
   // SubArray t(val.M());
   // for (int i=0;i<k;i++) {
   //   int j=comp[i];
   //   int ni=NN[i];
   //   int di=DF[i];
   //   int i1=i+1;
   //   int nii=NN[i1];
   //   int dii=DF[i1];
   //
   //   assert(ni<nii && di < dii);
   //   RNM_ v(val(SubArray(dii-di,di),t));
   //   if (j<=i) {
   //     teb[i]->FB(whatd,K,P,v);
   //     v = mappDx(J, v);
   //   }
   //   else
   //   v=val(SubArray(DF[j+1]-DF[j],DF[j]),t);
   // }
}

template <class Mesh>
void GTypeOfFESum<Mesh>::get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                                       KN_<double> &v) const {
   int k0 = 0;
   for (int i = 0; i < k; i++) {
      int n = teb[i]->ipj_Pi_h.N();
      KN_<R> sv(v(SubArray(n, k0)));
      teb[i]->get_Coef_Pi_h(K, sv);
      k0 += n;
   }
   assert(this->ipj_Pi_h.N() == k0);
}

// template<class Mesh>
// void GTypeOfFESum<Mesh>::init(InterpolationMatrix<RdHat> & M) const
// {
//   teb[0]->init(M);
// }

#endif
