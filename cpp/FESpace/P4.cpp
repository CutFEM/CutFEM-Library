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

#include "FESpace.hpp"

class TypeOfFE_P4Lagrange2d : public GTypeOfFE<Mesh2> {
   typedef Mesh2 Mesh;
   typedef typename Mesh::Element Element;

 public:
   static const int k   = 4;
   static const int ndf = (k + 2) * (k + 1) / 2;
   static int Data[];
   // static double alpha_Pi_h[];
   static const int nn[15][4];
   static const int aa[15][4];
   static const int ff[15];
   static const int il[15];
   static const int jl[15];
   static const int kl[15];

   TypeOfFE_P4Lagrange2d() : GTypeOfFE<Mesh2>(15, 1, Data, 15 + 6, 15, 0) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P4;
      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R2 Pt[15] = {
          R2(0 / 4., 0 / 4.), R2(4 / 4., 0 / 4.), R2(0 / 4., 4 / 4.),
          R2(3 / 4., 1 / 4.), R2(2 / 4., 2 / 4.), R2(1 / 4., 3 / 4.),
          R2(0 / 4., 3 / 4.), R2(0 / 4., 2 / 4.), R2(0 / 4., 1 / 4.),
          R2(1 / 4., 0 / 4.), R2(2 / 4., 0 / 4.), R2(3 / 4., 0 / 4.),
          R2(1 / 4., 2 / 4.), R2(2 / 4., 1 / 4.), R2(1 / 4., 1 / 4.)};

      int other[15] = {0, 1, 2, 5, 4, 3, 8, 7, 6, 11, 10, 9, 12, 13, 14};
      int kk        = 0;

      for (int i = 0; i < ndf; i++) {
         ipj_Pi_h[kk++] = IPJ(i, i, 0);
         if (other[i] != i) {
            ipj_Pi_h[kk++] = IPJ(i, other[i], 0);
         }

         Pt_Pi_h[i] = Pt[i];
      }

      assert(Pt_Pi_h.N() == ndf);
      assert(ipj_Pi_h.N() == kk);
   }

   void FB(const What_d, const Element &K, const Rd &PHat, RNMK_ &val) const;
   void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const {
      for (int i = 0; i < 15 + 6; ++i) {
         v[i] = 1;
      }

      int e0     = K.EdgeOrientation(0);
      int e1     = K.EdgeOrientation(1);
      int e2     = K.EdgeOrientation(2);
      int ooo[6] = {e0, e0, e1, e1, e2, e2};
      /*   3,4
       *   5,
       *   6,7
       *   8,9,
       *   10,
       *   11,12,
       *   13,14,
       *   15
       *   16,17
       */
      int iii[6] = {3, 6, 8, 11, 13, 16};
      int jjj[6] = {};

      for (int i = 0; i < 6; ++i) {
         jjj[i] = iii[i] + 1; // si orient = -1
      }

      for (int i = 0; i < 6; ++i) {
         if (ooo[i] == 1) {
            v[jjj[i]] = 0;
         } else {
            v[iii[i]] = 0;
         }
      }
   }
};

// on what     nu df on node node of df
int TypeOfFE_P4Lagrange2d::Data[] = {
    0, 1, 2,  3,  3,  3,  4,  4,
    4, 5, 5,  5,  6,  6,  6, // the support number  of the node of the df
    0, 0, 0,  0,  1,  2,  0,  1,
    2, 0, 1,  2,  0,  1,  2, // the number of the df on  the node
    0, 1, 2,  3,  3,  3,  4,  4,
    4, 5, 5,  5,  6,  6,  6, // the node of the df
    0, 1, 2,  3,  4,  5,  6,  7,
    8, 9, 10, 11, 12, 13, 14, // which are de df on sub FE
    1, 1, 1,  0,              // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, 15};

void TypeOfFE_P4Lagrange2d::FB(const What_d whatd, const Element &K,
                               const RdHat &PHat, RNMK_ &val) const {
   R2 A(K[0]), B(K[1]), C(K[2]);
   R l0 = 1. - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
   R L[3] = {l0 * k, l1 * k, l2 * k};

   assert(val.N() >= 10);
   assert(val.M() == 1);
   // Attention il faut renumeroter les fonction de bases
   // car dans freefem++, il y a un node par sommet, arete or element
   // et la numerotation naturelle  mais 2 noud pas arete
   // donc p est la perumation
   // echange de numerotation si les arete sont dans le mauvais sens
   int p[15] = {};

   for (int i = 0; i < 15; ++i) {
      p[i] = i;
   }

   if (K.EdgeOrientation(0) < 0) {
      std::swap(p[3], p[5]); // 3,4
   }

   if (K.EdgeOrientation(1) < 0) {
      std::swap(p[6], p[8]); // 5,6
   }

   if (K.EdgeOrientation(2) < 0) {
      std::swap(p[9], p[11]); // 7,8
   }

   val = 0;
   /*
    * //  les fonction de base du Pk Lagrange sont
    */

   if (whatd & Fop_D0) {
      RN_ f0(val('.', 0, op_id));

      for (int df = 0; df < ndf; df++) {
         int pdf = p[df];
         R f     = 1. / ff[df];

         for (int i = 0; i < k; ++i) {
            f *= L[nn[df][i]] - aa[df][i];
         }

         f0[pdf] = f;
      }
   }

   if (whatd & Fop_D1 || whatd & Fop_D2) {

      // if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||
      // whatd[op_dxy]) {
      R2 D[] = {K.H(0) * k, K.H(1) * k, K.H(2) * k};
      if (whatd & Fop_D1) {
         for (int df = 0; df < ndf; df++) {
            int pdf = p[df];
            R fx = 0., fy = 0., f = 1. / ff[df];

            for (int i = 0; i < k; ++i) {
               int n = nn[df][i];
               R Ln  = L[n] - aa[df][i];
               fx    = fx * Ln + f * D[n].x;
               fy    = fy * Ln + f * D[n].y;
               f     = f * Ln;
            }
            val(pdf, 0, op_dx) = fx;
            val(pdf, 0, op_dy) = fy;
         }
      }

      if (whatd & Fop_D2) {
         for (int df = 0; df < ndf; df++) {
            int pdf = p[df];
            R fx = 0., fy = 0., f = 1. / ff[df];
            R fxx = 0., fyy = 0., fxy = 0.;

            for (int i = 0; i < k; ++i) {
               int n = nn[df][i];
               R Ln  = L[n] - aa[df][i];
               fxx   = fxx * Ln + 2. * fx * D[n].x;
               fyy   = fyy * Ln + 2. * fy * D[n].y;
               fxy   = fxy * Ln + fx * D[n].y + fy * D[n].x;
               fx    = fx * Ln + f * D[n].x;
               fy    = fy * Ln + f * D[n].y;
               f     = f * Ln;
            }

            val(pdf, 0, op_dxx) = fxx;
            val(pdf, 0, op_dyy) = fyy;
            val(pdf, 0, op_dxy) = fxy;
         }
      }
   }
}

const int TypeOfFE_P4Lagrange2d::nn[15][4] = {
    {0, 0, 0, 0}, {1, 1, 1, 1}, {2, 2, 2, 2}, {1, 1, 1, 2}, {1, 1, 2, 2},
    {1, 2, 2, 2}, {0, 2, 2, 2}, {0, 0, 2, 2}, {0, 0, 0, 2}, {0, 0, 0, 1},
    {0, 0, 1, 1}, {0, 1, 1, 1}, {0, 1, 2, 2}, {0, 1, 1, 2}, {0, 0, 1, 2}};
const int TypeOfFE_P4Lagrange2d::aa[15][4] = {
    {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 0}, {0, 1, 0, 1},
    {0, 0, 1, 2}, {0, 0, 1, 2}, {0, 1, 0, 1}, {0, 1, 2, 0}, {0, 1, 2, 0},
    {0, 1, 0, 1}, {0, 0, 1, 2}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}};
const int TypeOfFE_P4Lagrange2d::ff[15] = {24, 24, 24, 6, 4, 6, 6, 4,
                                           6,  6,  4,  6, 2, 2, 2};
const int TypeOfFE_P4Lagrange2d::il[15] = {4, 0, 0, 0, 0, 0, 1, 2,
                                           3, 3, 2, 1, 1, 1, 2};
const int TypeOfFE_P4Lagrange2d::jl[15] = {0, 4, 0, 3, 2, 1, 0, 0,
                                           0, 1, 2, 3, 1, 2, 1};
const int TypeOfFE_P4Lagrange2d::kl[15] = {0, 0, 4, 1, 2, 3, 3, 2,
                                           1, 0, 0, 0, 2, 1, 1};

static TypeOfFE_P4Lagrange2d P4_2d;
GTypeOfFE<Mesh2> &P4Lagrange2d(P4_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P4 = P4_2d;
