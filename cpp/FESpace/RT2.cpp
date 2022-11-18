#include "RT1.cpp"

class TypeOfFE_RT2_2d : public InitTypeOfRTk_2d, public GTypeOfFE<Mesh2> {
   typedef Mesh2 Mesh;               // Define 2D mesh as Mesh
   typedef typename Mesh::Element E; // Define mesh element as E

 public:
   static double Pi_h_coef[];
   bool Ortho;

   TypeOfFE_RT2_2d(bool ortho)
       : InitTypeOfRTk_2d(2), GTypeOfFE<Mesh2>(
                                  ndf, 2, Data,
                                  2 * 3 * 3 * QFE.n +
                                      QFK.n * 4 * 3, // nb coef mat interpole
                                  3 * QFE.n + QFK.n, // nb P interpolation
                                  0),
         Ortho(ortho) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::RT2;
      GTypeOfFE<Mesh>::polynomialOrder = 3;

      int dofE  = this->k + 1;             // == 3
      int dofKs = (dofE - 1) * (dofE) / 2; //= = 3 ..(

      ffassert(dofKs == 3);
      ffassert(dofE == 3);

      int kkk = 0, i = 0;

      for (int e = 0; e < 3; ++e) {
         for (int p = 0; p < QFE.n; ++p) {
            R2 A(TriangleHat[Element::nvedge[e][0]]);
            R2 B(TriangleHat[Element::nvedge[e][1]]);

            for (int l = 0; l < dofE; ++l) {
               ipj_Pi_h[kkk++] = IPJ(dofE * e + l, i, 0);
               ipj_Pi_h[kkk++] = IPJ(dofE * e + l, i, 1);
            }

            Pt_Pi_h[i++] =
                B * (QFE[p].X()) + A * (1. - QFE[p].X()); // X=0 => A  X=1 => B;
         }
      }

      for (int p = 0; p < QFK.n; ++p) {
         int i6 = 3 * 3, i7 = i6 + 1;

         for (int l = 0; l < dofKs; ++l) {
            ipj_Pi_h[kkk++] = IPJ(i6, i, 0);
            ipj_Pi_h[kkk++] = IPJ(i6, i, 1);
            ipj_Pi_h[kkk++] = IPJ(i7, i, 0);
            ipj_Pi_h[kkk++] = IPJ(i7, i, 1);
            i6 += 2;
            i7 += 2;
         }

         Pt_Pi_h[i++] = QFK[p];
      }

      assert(kkk == this->ipj_Pi_h.N());
      assert(i == this->Pt_Pi_h.N());
   }

   void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const {
      const Element &T(K.T);
      int k          = 0;
      // magic fom:
      // inv of [[4!,3!,2!2!],[3!,2!2!,2!],[2!2!,3!,4!]]/5! =
      double c1[][3] = {
          {9, -18, 3} /* 0 */, {-18, 84, -18} /* 1 */, {3, -18, 9} /* 2 */};

      for (int i = 0; i < 3; i++) {
         R2 E(Ortho ? T.Edge(i) : -T.Edge(i).perp());
         R s = T.EdgeOrientation(i);

         for (int p = 0; p < QFE.n; ++p) {
            R l1 = QFE[p].X(), l2 = 1 - QFE[p].X();
            R l11 = l1 * l1;
            R l22 = l2 * l2;
            R l21 = l2 * l1;
            R p0  = c1[0][0] * l11 + c1[0][1] * l21 + c1[0][2] * l22; //
            R p1  = c1[1][0] * l11 + c1[1][1] * l21 + c1[1][2] * l22; //
            R p2  = c1[2][0] * l11 + c1[2][1] * l21 + c1[2][2] * l22; //
            R sa  = s * QFE[p].a;
            R cc2 = sa * p0; //
            R cc1 = sa * p1; //
            R cc0 = sa * p2; //
            if (s < 0) {
               std::swap(cc0, cc2);
            }

            v[k++] = cc0 * E.x;
            v[k++] = cc0 * E.y;
            v[k++] = cc1 * E.x;
            v[k++] = cc1 * E.y;
            v[k++] = cc2 * E.x;
            v[k++] = cc2 * E.y;
         }
      }

      R2 B[2] = {T.Edge(1), T.Edge(2)};
      if (Ortho) {
         B[0] = -B[0];
         B[1] = -B[1];
      } else {
         B[0] = B[0].perp();
         B[1] = B[1].perp();
      }

      double CK = 0.5; // dof U= [u1,u2] > |K| int_K ( B_i.U )
      R ll[3];         //, lo[3];

      for (int p = 0; p < QFK.n; ++p) {
         double w = -QFK[p].a * CK;

         ll[0] = w * (1. - QFK[p].x - QFK[p].y);
         ll[1] = w * QFK[p].x;
         ll[2] = w * QFK[p].y;

         for (int l = 0; l < 3; ++l) {
            v[k++] = ll[l] * B[0].x;
            v[k++] = ll[l] * B[0].y;
            v[k++] = ll[l] * B[1].x;
            v[k++] = ll[l] * B[1].y;
         }
      }

      assert(k == this->ipj_Pi_h.N());
   }

   void FB(const What_d whatd, const Element &K, const Rd &Phat,
           RNMK_ &val) const;

 private:
   void FB_Freefem(const What_d, const Element &K, const Rd &PHat,
                   RNMK_ &bfMat) const;
   // void FB_ID(const Element &K, const Rd &Phat, RNMK_ &bfMat) const ;
   // void FB_D1(const Element &K, const Rd &Phat, RNMK_ &bfMat) const ;
   void FB_D2(const Element &K, const Rd &Phat, RNMK_ &bfMat) const;
};

void TypeOfFE_RT2_2d::FB(const What_d whatd, const Element &K, const Rd &Phat,
                         RNMK_ &bfMat) const {

   assert(bfMat.N() >= ndf);
   assert(bfMat.M() == 2);
   bfMat = 0;

   if ((whatd & Fop_D2)) {
      FB_D2(K, Phat, bfMat);
   }
   // else if ((whatd & Fop_D1) ) { FB_D1(K, Phat, bfMat);}
   // else {FB_ID(K, Phat, bfMat);}
   else {
      FB_Freefem(whatd, K, Phat, bfMat);
   }
}
void TypeOfFE_RT2_2d::FB_Freefem(const What_d whatd, const Element &K,
                                 const R2 &Phat, RNMK_ &val) const {
   R2 X   = K(Phat);
   R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
   R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y;
   R L[3] = {l0, l1, l2};
   R eo[] = {static_cast<R>(K.EdgeOrientation(0)),
             static_cast<R>(K.EdgeOrientation(1)),
             static_cast<R>(K.EdgeOrientation(2))};

   int p[15] = {0, 1, 2,  5,  4,  3,  6, 7,
                8, 9, 10, 11, 12, 13, 14}; // Permutation for orinatation
   R2 Pm[18];                              // all the momome function ..
   double cf[][6] = {{0, -5.5, 0, -2.5, -0.5, -1.5} /* 0 */,
                     {-1.25, -1.25, 0.25, -1, 0.25, -1} /* 1 */,
                     {-5.5, 0, -0.5, -1.5, 0, -2.5} /* 2 */,
                     {0, -2.5, 0, -5.5, -1.5, -0.5} /* 3 */,
                     {0.25, -1, -1.25, -1.25, -1, 0.25} /* 4 */,
                     {-0.5, -1.5, -5.5, 0, -2.5, 0} /* 5 */,
                     {-2.5, 0, -1.5, -0.5, 0, -5.5} /* 6 */,
                     {-1, 0.25, -1, 0.25, -1.25, -1.25} /* 7 */,
                     {-1.5, -0.5, -2.5, 0, -5.5, 0} /* 8 */,
                     {30, 90, -30, 180, 30, 60} /* 9 */,
                     {90, 30, 30, 60, -30, 180} /* 10 */,
                     {30, -180, -30, -90, -60, -30} /* 11 */,
                     {60, -120, 60, -60, 120, -60} /* 12 */,
                     {-120, 60, 120, -60, 60, -60} /* 13 */,
                     {-180, 30, -60, -30, -30, -90} /* 14 */};
   int Bii[][3]   = {{0, 0, 0} /* 0 */,  {0, 1, 1} /* 1 */,  {0, 2, 2} /* 2 */,
                     {0, 1, 2} /* 3 */,  {0, 2, 0} /* 4 */,  {0, 0, 1} /* 5 */,
                     {1, 0, 0} /* 6 */,  {1, 1, 1} /* 7 */,  {1, 2, 2} /* 8 */,
                     {1, 1, 2} /* 9 */,  {1, 2, 0} /* 10 */, {1, 0, 1} /* 11 */,
                     {2, 0, 0} /* 12 */, {2, 1, 1} /* 13 */, {2, 2, 2} /* 14 */,
                     {2, 1, 2} /* 15 */, {2, 2, 0} /* 16 */, {2, 0, 1} /* 17 */};
   int fe[]       = {1, 3, 2, 6, 10, 8, 12, 17, 13};
   int k6[]       = {4, 5, 9, 11, 15, 16};
   R CKK          = K.mesure() * 2;
   R2 phi[3]      = {X - Q[0], X - Q[1], X - Q[2]}; // phi * area *2

   for (int l = 0; l < 18; ++l) {
      int i = Bii[l][0];
      int j = Bii[l][1];
      int k = Bii[l][2];
      Pm[l] = phi[i] * (L[j] * L[k] / CKK);
   }

   // static int ddd=0;

   if (eo[0] < 0) {
      std::swap(p[0], p[2]);
   }

   if (eo[1] < 0) {
      std::swap(p[3], p[5]);
   }

   if (eo[2] < 0) {
      std::swap(p[6], p[8]);
   }

   double sg[15] = {eo[0], eo[0], eo[0], eo[1], eo[1], eo[1], eo[2], eo[2],
                    eo[2], 1.,    1.,    1.,    1.,    1.,    1.};

   if (whatd & Fop_id) {
      for (int pdf = 0; pdf < 15; ++pdf) {
         int df = p[pdf];
         R2 fd(0., 0.);

         if (df < 9) {
            fd = Pm[fe[df]]; // edge function ..
         }

         for (int k = 0; k < 6; ++k) {
            fd += cf[df][k] * Pm[k6[k]];
         }

         fd *= sg[df];

         val(pdf, 0, op_id) = fd.x;
         val(pdf, 1, op_id) = fd.y;

         // std::cout << FDx[pdf].val.val << "\t" << fd.x << std::endl;
         // std::cout << FDy[pdf].val.val << "\t" << fd.y << std::endl;
      }
   }
   if ((whatd & Fop_D1) || (whatd & Fop_D2)) {
      R2 DL[3] = {K.H(0), K.H(1), K.H(2)};
      R2 Dphi1(1, 0); // D(phi.x)
      R2 Dphi2(0, 1); // D(phi.y)
      R2 DxPm[18];
      R2 DyPm[18];

      if (Ortho) {
         Dphi1 = R2(0, -1);
         Dphi2 = R2(1, 0);
      } // Correction Jan 2019 FH. 	// x,y -> (-y,x)

      for (int l = 0; l < 18; ++l) {
         int i   = Bii[l][0];
         int j   = Bii[l][1];
         int k   = Bii[l][2];
         R Ljk   = L[j] * L[k];
         R2 DLjk = L[j] * DL[k] + DL[j] * L[k];
         R2 DF1 =
             (Dphi1 * Ljk + phi[i].x * DLjk) / CKK; // BUG ici Ortho ????????
         R2 DF2  = (Dphi2 * Ljk + phi[i].y * DLjk) / CKK;
         DxPm[l] = R2(DF1.x, DF2.x);
         DyPm[l] = R2(DF1.y, DF2.y);
      }

      if (whatd & Fop_dx) {
         for (int pdf = 0; pdf < 15; ++pdf) {
            int df = p[pdf];
            R2 fd(0., 0.);
            if (df < 9) {
               fd = DxPm[fe[df]]; // edge function ..
            }
            for (int k = 0; k < 6; ++k) {
               fd += cf[df][k] * DxPm[k6[k]];
            }

            fd *= sg[df];
            val(pdf, 0, op_dx) = fd.x;
            val(pdf, 1, op_dx) = fd.y;
         }
      }
      if (whatd & Fop_dy) {
         for (int pdf = 0; pdf < 15; ++pdf) {
            int df = p[pdf];
            R2 fd(0., 0.);
            if (df < 9) {
               fd = DyPm[fe[df]]; // edge function ..
            }

            for (int k = 0; k < 6; ++k) {
               fd += cf[df][k] * DyPm[k6[k]];
            }

            fd *= sg[df];
            val(pdf, 0, op_dy) = fd.x;
            val(pdf, 1, op_dy) = fd.y;
         }
      }

      if (whatd & Fop_D2) {

         std::cout << " to do FH RT2 dxx, dyy dxy " << std::endl;
         ffassert(0);
      }
   }
}
void TypeOfFE_RT2_2d::FB_D2(const Element &K, const R2 &Phat,
                            RNMK_ &bfMat) const {
   R2 X   = K(Phat);
   R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
   R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y;
   R L[3]   = {l0, l1, l2};
   R eo[]   = {static_cast<R>(K.EdgeOrientation(0)),
               static_cast<R>(K.EdgeOrientation(1)),
               static_cast<R>(K.EdgeOrientation(2))};
   double s = 1.;

   R2 ddd[3] = {K.H(0), K.H(1), K.H(2)};
   DDiff_R2 ll1(l1, 0), ll2(l2, 1);
   ll1.d[0]     = ddd[1].x;
   ll1.d[1]     = ddd[1].y;
   ll2.d[0]     = ddd[2].x;
   ll2.d[1]     = ddd[2].y;
   ll1.val.d[0] = ddd[1].x;
   ll1.val.d[1] = ddd[1].y; // init val of dx
   ll2.val.d[0] = ddd[2].x;
   ll2.val.d[1] = ddd[2].y;

   DDiff_R2 ll0 = 1 - ll1 - ll2;
   DDiff_R2 Xx(X.x, 0), Xy(X.y, 1);
   Xx.val.d[0]      = 1; // init val of dx
   Xy.val.d[1]      = 1;
   DDiff_R2 LL[3]   = {ll0, ll1, ll2};
   DDiff_R2 PHIx[3] = {
       Xx - Q[0].x,
       Xx - Q[1].x,
       Xx - Q[2].x,
   };
   DDiff_R2 PHIy[3] = {
       Xy - Q[0].y,
       Xy - Q[1].y,
       Xy - Q[2].y,
   };
   DDiff_R2 PMx[18];
   DDiff_R2 PMy[18];

   int p[15]      = {0, 1, 2,  5,  4,  3,  6, 7,
                     8, 9, 10, 11, 12, 13, 14}; // Permutation for orinatation
   // R2 Pm[18];                                 // all the momome function ..
   double cf[][6] = {{0, -5.5, 0, -2.5, -0.5, -1.5} /* 0 */,
                     {-1.25, -1.25, 0.25, -1, 0.25, -1} /* 1 */,
                     {-5.5, 0, -0.5, -1.5, 0, -2.5} /* 2 */,
                     {0, -2.5, 0, -5.5, -1.5, -0.5} /* 3 */,
                     {0.25, -1, -1.25, -1.25, -1, 0.25} /* 4 */,
                     {-0.5, -1.5, -5.5, 0, -2.5, 0} /* 5 */,
                     {-2.5, 0, -1.5, -0.5, 0, -5.5} /* 6 */,
                     {-1, 0.25, -1, 0.25, -1.25, -1.25} /* 7 */,
                     {-1.5, -0.5, -2.5, 0, -5.5, 0} /* 8 */,
                     {30, 90, -30, 180, 30, 60} /* 9 */,
                     {90, 30, 30, 60, -30, 180} /* 10 */,
                     {30, -180, -30, -90, -60, -30} /* 11 */,
                     {60, -120, 60, -60, 120, -60} /* 12 */,
                     {-120, 60, 120, -60, 60, -60} /* 13 */,
                     {-180, 30, -60, -30, -30, -90} /* 14 */};
   int Bii[][3]   = {{0, 0, 0} /* 0 */,  {0, 1, 1} /* 1 */,  {0, 2, 2} /* 2 */,
                     {0, 1, 2} /* 3 */,  {0, 2, 0} /* 4 */,  {0, 0, 1} /* 5 */,
                     {1, 0, 0} /* 6 */,  {1, 1, 1} /* 7 */,  {1, 2, 2} /* 8 */,
                     {1, 1, 2} /* 9 */,  {1, 2, 0} /* 10 */, {1, 0, 1} /* 11 */,
                     {2, 0, 0} /* 12 */, {2, 1, 1} /* 13 */, {2, 2, 2} /* 14 */,
                     {2, 1, 2} /* 15 */, {2, 2, 0} /* 16 */, {2, 0, 1} /* 17 */};
   int fe[]       = {1, 3, 2, 6, 10, 8, 12, 17, 13};
   int k6[]       = {4, 5, 9, 11, 15, 16};
   R CKK          = K.mesure() * 2;
   // R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]};    // phi * area *2

   for (int l = 0; l < 18; ++l) {
      int i  = Bii[l][0];
      int j  = Bii[l][1];
      int k  = Bii[l][2];
      // Pm[l] = phi[i] * (L[j] * L[k] / CKK);
      PMx[l] = PHIx[i] * (LL[j] * LL[k] / CKK);
      PMy[l] = PHIy[i] * (LL[j] * LL[k] / CKK);
   }
   if (eo[0] < 0) {
      std::swap(p[0], p[2]);
   }
   if (eo[1] < 0) {
      std::swap(p[3], p[5]);
   }
   if (eo[2] < 0) {
      std::swap(p[6], p[8]);
   }

   double sg[15] = {eo[0], eo[0], eo[0], eo[1], eo[1], eo[1], eo[2], eo[2],
                    eo[2], 1.,    1.,    1.,    1.,    1.,    1.};

   for (int pdf = 0; pdf < 15; ++pdf) {
      int df = p[pdf];
      // R2 fd(0., 0.);
      DDiff_R2 FDx, FDy;
      if (df < 9) {
         // fd = Pm[fe[df]];    // edge function ..
         FDx = PMx[fe[df]];
         FDy = PMy[fe[df]];
      }

      for (int k = 0; k < 6; ++k) {
         // fd  += cf[df][k] * Pm[k6[k]];
         FDx += cf[df][k] * PMx[k6[k]];
         FDy += cf[df][k] * PMy[k6[k]];
      }

      // fd *= sg[df];
      FDx *= sg[df];
      FDy *= sg[df];

      bfMat(pdf, 0, op_id)  = FDx.val.val * s;
      bfMat(pdf, 1, op_id)  = FDy.val.val * s;
      bfMat(pdf, 0, op_dx)  = FDx.d[0].val * s;
      bfMat(pdf, 1, op_dx)  = FDy.d[0].val * s;
      bfMat(pdf, 0, op_dy)  = FDx.d[1].val * s;
      bfMat(pdf, 1, op_dy)  = FDy.d[1].val * s;
      bfMat(pdf, 0, op_dxx) = FDx.d[0].d[0] * s;
      bfMat(pdf, 1, op_dxx) = FDy.d[0].d[0] * s;
      bfMat(pdf, 0, op_dxy) = FDx.d[0].d[1] * s;
      bfMat(pdf, 1, op_dxy) = FDy.d[0].d[1] * s;
      bfMat(pdf, 0, op_dyy) = FDx.d[1].d[1] * s;
      bfMat(pdf, 1, op_dyy) = FDy.d[1].d[1] * s;
   }
}

static TypeOfFE_RT2_2d myRT2_2d(false);
GTypeOfFE<Mesh2> &RT2_2d(myRT2_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::RT2 = myRT2_2d;
