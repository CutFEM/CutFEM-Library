#include "FESpace.hpp"
#include "../util/ufunction.hpp"


class TypeOfFE_P2BRLagrange : public GTypeOfFE<Mesh2> {
  typedef  Mesh2 Mesh; // Define 2D mesh as Mesh
  typedef  typename Mesh::Element  E;

  public:
   static int Data[];

   TypeOfFE_P2BRLagrange( )
     : GTypeOfFE<Mesh2>(6 + 3 + 0,          // Number of degrees of freedom
                2,                  // Dimension N
                Data,               // Data array
                6 + 3 * (2 + 2),    // Number kPi of coefficients to build the interpolation
                9                  // number nPi of integration points to build the interpolation
       ) {
     const double gauss1 = (1. - sqrt(1. / 3.)) / 2;
     const double gauss2 = 1. - gauss1;
     const R2 Pt[] = {R2(0, 0), R2(1, 0), R2(0, 1)};

     // For the 3 vertices: 6 coefficents
     int kk = 0;

     for (int p = 0; p < 3; p++) {
       Pt_Pi_h[p] = Pt[p];
       ipj_Pi_h[kk] = IPJ(kk, p, 0);
       ++kk;
       ipj_Pi_h[kk] = IPJ(kk, p, 1);
       ++kk;
     }

     // Integration point on edge e
     int p = 3;

     for (int e = 0; e < 3; ++e) {
       R2 A = Pt[Element::nvedge[e][0]];
       R2 B = Pt[Element::nvedge[e][1]];
       Pt_Pi_h[p] = A * gauss1 + B * gauss2;
       ipj_Pi_h[kk++] = IPJ(6 + e, p, 0);    // coef = 0.5 * l_e * ne_x * sge
       ipj_Pi_h[kk++] = IPJ(6 + e, p, 1);    // coef = 0.5 * l_e * ne_y * sge
       p++;
       Pt_Pi_h[p] = A * gauss2 + B * gauss1;
       ipj_Pi_h[kk++] = IPJ(6 + e, p, 0);    // coef = 0.5 * l_e * ne_x * sge
       ipj_Pi_h[kk++] = IPJ(6 + e, p, 1);    // coef = 0.5 * l_e * ne_y * sge
       p++;
     }

     assert(Pt_Pi_h.N( ) == p);
     assert(ipj_Pi_h.N( ) == kk);
   }

   void FB(const What_d, const Element &K, const Rd &PHat, RNMK_ &val) const;
   void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_< double > &v) const;
 };

 // Data array
 int TypeOfFE_P2BRLagrange::Data[] = {
   0, 0, 1, 1, 2, 2, 3, 4, 5,
   0, 1, 0, 1, 0, 1, 0, 0, 0,
   0, 0, 1, 1, 2, 2, 3, 4, 5,
   0, 1, 2, 3, 4, 5, 6, 7, 8,
   1, 1, 0, 0,
   0,
   0,
   9};

 /*!
  * \brief Define alpha k
  * \param K const baseFElement &
  * \param v KN<double> &
  */
 // Define alpha_k
 void TypeOfFE_P2BRLagrange::get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_< double > &v) const {
   const Element &T(K.T);
   int k = 0;

   // Coefficents for the 3 vertices time the 2 components
   for (int i = 0; i < 6; i++) {
     v[k++] = 1;
   }

   // Integration on edges
   for (int i = 0; i < 3; i++) {
     R2 N(T.Edge(i).perp( ));
     N *= T.EdgeOrientation(i) * 0.5;
     v[k++] = N.x;
     v[k++] = N.y;
     v[k++] = N.x;
     v[k++] = N.y;
   }
 }

 /*!
  * \brief Shape function
  * \param whatd const bool*
  * \param const Mesh &
  * \param K const Triangle &
  * \param PHat const RdHat &
  * \param val RNMK_
  */
 // Shape function
 void TypeOfFE_P2BRLagrange::FB(const What_d whatd, const Element &K, const Rd &PHat, RNMK_ &val) const {
   int max_op = 0;
   R2 A(K[0]), B(K[1]), C(K[2]);
   R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
   // int_e_1 l0*l0 = |e_1|/3 and int_e_1 l0*l1 = |e_1|/6
   // to get the flux = 1
   R2 E[3] = {K.Edge(0), K.Edge(1), K.Edge(2)};
   double l2E[3] = {(E[0], E[0]), (E[1], E[1]), (E[2], E[2])};
   // double lE[3] = {sqrt(l2E[0]), sqrt(l2E[1]), sqrt(l2E[2])};
   double sgE[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
   R2 cN[3] = {E[0].perp( ) * (6. * sgE[0] / l2E[0]), E[1].perp( ) * (6. * sgE[1] / l2E[1]),
               E[2].perp( ) * (6. * sgE[2] / l2E[2])};


   assert(val.N( ) >= 9);
   assert(val.M( ) == 2);
   val = 0;

   if (whatd & Fop_id) {
     max_op = 1;
     RN_ f0(val('.', 0, op_id));
     RN_ f1(val('.', 1, op_id));

     f1[1] = f0[0] = l0;
     f1[3] = f0[2] = l1;
     f1[5] = f0[4] = l2;

     f0[6] = cN[0].x * l1 * l2;    // opposite to the vertex 0
     f0[7] = cN[1].x * l0 * l2;    // opposite to the vertex 1
     f0[8] = cN[2].x * l1 * l0;    // opposite to the vertex 2

     f1[6] = cN[0].y * l1 * l2;    // opposite to the vertex 0
     f1[7] = cN[1].y * l0 * l2;    // opposite to the vertex 1
     f1[8] = cN[2].y * l1 * l0;    // opposite to the vertex 2
   }

   if ((whatd & Fop_D1) || (whatd & Fop_D2) ) {
     R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));

     if (whatd & Fop_dx) {
       max_op = 4;
       RN_ f0x(val('.', 0, op_dx));
       RN_ f1x(val('.', 1, op_dx));

       f1x[1] = f0x[0] = Dl0.x;
       f1x[3] = f0x[2] = Dl1.x;
       f1x[5] = f0x[4] = Dl2.x;

       f0x[6] = cN[0].x * (Dl1.x * l2 + Dl2.x * l1);
       f0x[7] = cN[1].x * (Dl2.x * l0 + Dl0.x * l2);
       f0x[8] = cN[2].x * (Dl0.x * l1 + Dl1.x * l0);

       f1x[6] = cN[0].y * (Dl1.x * l2 + Dl2.x * l1);
       f1x[7] = cN[1].y * (Dl2.x * l0 + Dl0.x * l2);
       f1x[8] = cN[2].y * (Dl0.x * l1 + Dl1.x * l0);
     }

     if (whatd & Fop_dy) {
       RN_ f0y(val('.', 0, op_dy));
       RN_ f1y(val('.', 1, op_dy));

       f1y[1] = f0y[0] = Dl0.y;
       f1y[3] = f0y[2] = Dl1.y;
       f1y[5] = f0y[4] = Dl2.y;

       f0y[6] = cN[0].x * (Dl1.y * l2 + Dl2.y * l1);
       f0y[7] = cN[1].x * (Dl2.y * l0 + Dl0.y * l2);
       f0y[8] = cN[2].x * (Dl0.y * l1 + Dl1.y * l0);

       f1y[6] = cN[0].y * (Dl1.y * l2 + Dl2.y * l1);
       f1y[7] = cN[1].y * (Dl2.y * l0 + Dl0.y * l2);
       f1y[8] = cN[2].y * (Dl0.y * l1 + Dl1.y * l0);
     }

     if (whatd & Fop_dxx) {
       max_op = 10;
       RN_ f0xx(val('.', 0, op_dxx));
       RN_ f1xx(val('.', 1, op_dxx));

       f0xx[6] = 2 * cN[0].x * Dl1.x * Dl2.x;
       f0xx[7] = 2 * cN[1].x * Dl0.x * Dl2.x;
       f0xx[8] = 2 * cN[2].x * Dl0.x * Dl1.x;
       f1xx[6] = 2 * cN[0].y * Dl1.x * Dl2.x;
       f1xx[7] = 2 * cN[1].y * Dl0.x * Dl2.x;
       f1xx[8] = 2 * cN[2].y * Dl0.x * Dl1.x;
     }

     if (whatd & Fop_dyy) {
       RN_ f0yy(val('.', 0, op_dyy));
       RN_ f1yy(val('.', 1, op_dyy));

       f0yy[6] = 2 * cN[0].x * Dl1.y * Dl2.y;
       f0yy[7] = 2 * cN[1].x * Dl0.y * Dl2.y;
       f0yy[8] = 2 * cN[2].x * Dl0.y * Dl1.y;
       f1yy[6] = 2 * cN[0].y * Dl1.y * Dl2.y;
       f1yy[7] = 2 * cN[1].y * Dl0.y * Dl2.y;
       f1yy[8] = 2 * cN[2].y * Dl0.y * Dl1.y;
     }

     if (whatd & Fop_dxy) {
       assert(val.K( ) > op_dxy);
       RN_ f0xy(val('.', 0, op_dxy));
       RN_ f1xy(val('.', 1, op_dxy));

       f0xy[6] = cN[0].x * (Dl1.x * Dl2.y + Dl1.y * Dl2.x);
       f0xy[7] = cN[1].x * (Dl0.x * Dl2.y + Dl0.y * Dl2.x);
       f0xy[8] = cN[2].x * (Dl0.x * Dl1.y + Dl0.y * Dl1.x);
       f1xy[6] = cN[0].y * (Dl1.x * Dl2.y + Dl1.y * Dl2.x);
       f1xy[7] = cN[1].y * (Dl0.x * Dl2.y + Dl0.y * Dl2.x);
       f1xy[8] = cN[2].y * (Dl0.x * Dl1.y + Dl0.y * Dl1.x);
     }
   }

   // Now, remove the flux part on 6 first DOF
   // w_i = w_i - a_i w_{k_i} - b_i w_{l_i}
   {
     int k[6] = {6 + 1, 6 + 1, 6 + 2, 6 + 2, 6 + 0, 6 + 0};
     int l[6] = {6 + 2, 6 + 2, 6 + 0, 6 + 0, 6 + 1, 6 + 1};
     R2 eN[3] = {E[0].perp( ) * (0.5 * sgE[0]), E[1].perp( ) * (0.5 * sgE[1]),
                 E[2].perp( ) * (0.5 * sgE[2])};
     double a[6] = {eN[1].x, eN[1].y, eN[2].x, eN[2].y, eN[0].x, eN[0].y};
     double b[6] = {eN[2].x, eN[2].y, eN[0].x, eN[0].y, eN[1].x, eN[1].y};
     int nop = 0;

     int vop[max_op] = {};


     for (int j = 0; j < max_op; j++) {
         vop[nop++] = j;
     }

     for (int i = 0; i < 6; ++i) {
       for (int jj = 0; jj < nop; ++jj) {
         int j = vop[jj];
         val(i, 0, j) -= a[i] * val(k[i], 0, j) + b[i] * val(l[i], 0, j);
         val(i, 1, j) -= a[i] * val(k[i], 1, j) + b[i] * val(l[i], 1, j);
       }
     }
   }
 }

 static TypeOfFE_P2BRLagrange  myP2BR_2d;
 GTypeOfFE<Mesh2> & P2BR_2d(myP2BR_2d);
 template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P2BR=myP2BR_2d;
