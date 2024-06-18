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

// P1
class TypeOfFE_P3dcLagrange2d : public GTypeOfFE<Mesh2> {

   typedef Mesh2 Mesh;
   typedef typename Mesh::Element E;
   // static const int nbNodeOnItem[4];
 public:
   static const int k   = 3;
   static const int ndf = (k + 2) * (k + 1) / 2;
   static int Data[];
   static double alpha_Pi_h[];
   static const int nn[10][3];
   static const int aa[10][3];
   static const int ff[10];
   static const int il[10];
   static const int jl[10];
   static const int kl[10];
   // dof, dim Im, Data, coefInterp, nbPt, coeff
   TypeOfFE_P3dcLagrange2d()
       : GTypeOfFE<Mesh2>(10, 1, Data, 10, 10, alpha_Pi_h) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P3dc;
      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R2 Pt[10] = {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.),
                                R2(0 / 3., 3 / 3.), R2(2 / 3., 1 / 3.),
                                R2(1 / 3., 2 / 3.), R2(0 / 3., 2 / 3.),
                                R2(0 / 3., 1 / 3.), R2(1 / 3., 0 / 3.),
                                R2(2 / 3., 0 / 3.), R2(1 / 3., 1 / 3.)};

      for (int i = 0; i < ndf; ++i) {
         Pt_Pi_h[i]  = Pt[i];
         ipj_Pi_h[i] = IPJ(i, i, 0);
      }
   }

   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

// const int TypeOfFE_P3dcLagrange2d::nbNodeOnItem[4] = {1,2,1,0};
int TypeOfFE_P3dcLagrange2d::Data[] = {
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, // the support number  of the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // the number of the df on  the node
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // which are de df on sub FE
    0, 0, 1, 0,                   // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, // begin_dfcomp
    10 // end_dfcomp
};

double TypeOfFE_P3dcLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1.,
                                                1., 1., 1., 1., 1.};

// const R2 TypeOfFE_P3dcLagrange2d::G(1. / 3., 1. / 3.);
// const R TypeOfFE_P3dcLagrange2d::cshrink = 1 - 1e-2;
// const R TypeOfFE_P3dcLagrange2d::cshrink1 = 1. /
// TypeOfFE_P3dcLagrange2d::cshrink;

void TypeOfFE_P3dcLagrange2d::FB(const What_d whatd, const Element &K,
                                 const Rd &PHat, RNMK_ &val) const {
   R2 A(K[0]), B(K[1]), C(K[2]);
   R l0 = 1. - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
   R L[3] = {l0 * k, l1 * k, l2 * k};
   assert(val.N() >= 10);
   assert(val.M() == 1);

   int p[10] = {};

   for (int i = 0; i < 10; ++i) {
      p[i] = i;
   }
   val = 0;
   RN_ f0(val('.', 0, op_id));

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

const int TypeOfFE_P3dcLagrange2d::nn[10][3] = {
    {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {1, 1, 2}, {1, 2, 2},
    {0, 2, 2}, {0, 0, 2}, {0, 0, 1}, {0, 1, 1}, {0, 1, 2}};
const int TypeOfFE_P3dcLagrange2d::aa[10][3] = {
    {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 0}, {0, 0, 1},
    {0, 0, 1}, {0, 1, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
const int TypeOfFE_P3dcLagrange2d::ff[10] = {6, 6, 6, 2, 2, 2, 2, 2, 2, 1};
const int TypeOfFE_P3dcLagrange2d::il[10] = {3, 0, 0, 0, 0, 1, 2, 2, 1, 1};
const int TypeOfFE_P3dcLagrange2d::jl[10] = {0, 3, 0, 2, 1, 0, 0, 1, 2, 1};
const int TypeOfFE_P3dcLagrange2d::kl[10] = {0, 0, 3, 1, 2, 2, 1, 0, 0, 1};

static TypeOfFE_P3dcLagrange2d P3dc_2d;
GTypeOfFE<Mesh2> &P3dcLagrange2d(P3dc_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P3dc = P3dc_2d;



// P3 Lagrange on quadrilaterals in 2D
class TypeOfFE_P3QdcLagrange2d : public GTypeOfFE<MeshQuad2> {
    typedef MeshQuad2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 3;
    static const int ndf = (k + 1) * (k + 1);
    static int Data[];
    //    static double alpha_Pi_h[];

    //! Unsure if the initialization list below is correct
    // TypeOfFE_P3QLagrange2d() : GTypeOfFE<MeshQuad2>(16, 1, Data, 20, 16, alpha_Pi_h) {
    TypeOfFE_P3QdcLagrange2d() : GTypeOfFE<MeshQuad2>(16, 1, Data, 24, 16, 0) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P3;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R2 Pt[16] = {R2(0., 0.),         R2(1., 0.),         R2(1., 1.),         R2(0., 1.),
                                  R2(1. / 3, 0.),     R2(2. / 3, 0.),     R2(1., 1. / 3),     R2(1., 2. / 3),
                                  R2(2. / 3, 1.),     R2(1. / 3, 1.),     R2(0., 2. / 3),     R2(0., 1. / 3),
                                  R2(1. / 3, 1. / 3), R2(2. / 3, 1. / 3), R2(2. / 3, 2. / 3), R2(1. / 3, 2. / 3)};

        int other[16] = {-1, -1, -1, -1, 5, 4, 7, 6, 9, 8, 11, 10, -1, -1, -1, -1};
        int kk        = 0;

        for (int i = 0; i < ndf; i++) {
            ipj_Pi_h[kk++] = IPJ(i, i, 0);
            if (other[i] >= 0) {
                ipj_Pi_h[kk++] = IPJ(i, other[i], 0);
            }

            Pt_Pi_h[i] = Pt[i];
        }

        // for (int i = 0; i < ndf; ++i) {
        //     Pt_Pi_h[i]  = Pt[i];
        //     ipj_Pi_h[i] = IPJ(i, i, 0);
        // }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const override;

    void get_Coef_Pi_h(const GbaseFElement<MeshQuad2> &K, KN_<double> &v) const override {
        for (int i = 0; i < 24; ++i) {
            v[i] = 1;
        }

        int e0     = K.EdgeOrientation(0);
        int e1     = K.EdgeOrientation(1);
        int e2     = K.EdgeOrientation(2);
        int e3     = K.EdgeOrientation(3);
        int ooo[8] = {e0, e0, e1, e1, e2, e2, e3, e3};
        int iii[8] = {};
        int jjj[8] = {};

        for (int i = 0; i < 8; ++i) {
            iii[i] = 4 + 2 * i; // si  orient = 1
            jjj[i] = 5 + 2 * i; // si orient = -1
        }

        for (int i = 0; i < 8; ++i) {
            if (ooo[i] == 1) {
                v[jjj[i]] = 0;
            } else {
                v[iii[i]] = 0;
            }
        }
    }
};

const int TypeOfFE_P3QdcLagrange2d::nbNodeOnItem[4] = {1, 1, 1, 0};

int TypeOfFE_P3QdcLagrange2d::Data[] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,       // the support number  of the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, // the number of the df on  the node
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, // which are de df on sub FE
    0, 0, 1, 0,                                           // nb node on what
    0,                                                    // for each compontant $j=0,N-1$ it give the sub FE associated
    0,                                                    // begin_dfcomp
    16                                                    // end_dfcomp
};

// double TypeOfFE_P3QLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

void TypeOfFE_P3QdcLagrange2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
        int p[16] = {};

    for (int i = 0; i < 16; ++i) {
        p[i] = i;
    }

    if (K.EdgeOrientation(0) < 0) {
        std::swap(p[4], p[5]); // 3,4
    }

    if (K.EdgeOrientation(1) < 0) {
        std::swap(p[6], p[7]); // 5,6
    }

    if (K.EdgeOrientation(2) < 0) {
        std::swap(p[8], p[9]); // 7,8
    }


    if (K.EdgeOrientation(3) < 0) {
        std::swap(p[10], p[11]); // 7,8
    }
    
    R x = P.x, y = P.y;

    std::array<double, 4> psi_x{-4.5 * x * x * x + 9 * x * x - 5.5 * x + 1, 0.5 * x * (9 * x * x - 9 * x + 2),
                                4.5 * x * (3 * x * x - 5 * x + 2), 4.5 * (-3 * x * x * x + 4 * x * x - x)};
    std::array<double, 4> psi_y{-4.5 * y * y * y + 9 * y * y - 5.5 * y + 1, 0.5 * y * (9 * y * y - 9 * y + 2),
                                4.5 * y * (3 * y * y - 5 * y + 2), 4.5 * (-3 * y * y * y + 4 * y * y - y)};

    R2 edge_x(K[0], K[1]); // [hx, 0]
    R2 edge_y(K[0], K[3]); // [0, hy]

    assert(val.N() >= E::nv + E::ne);
    assert(val.M() == 1);

    // Assert that points are on the reference element
    // assert((0. - 1e-12 <= x) && (x <= 1. + 1e-12));
    // assert((0. - 1e-12 <= y) && (y <= 1. + 1e-12));

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {

        // phi_l = psi_i^x * psi_j^y

        f0[0]  = psi_x[0] * psi_y[0];
        f0[1]  = psi_x[1] * psi_y[0];
        f0[2]  = psi_x[1] * psi_y[1];
        f0[3]  = psi_x[0] * psi_y[1];
        f0[p[4]]  = psi_x.at(2) * psi_y.at(0);
        f0[p[5]]  = psi_x.at(3) * psi_y.at(0);
        f0[p[6]]  = psi_x.at(1) * psi_y.at(2);
        f0[p[7]]  = psi_x.at(1) * psi_y.at(3);
        f0[p[8]]  = psi_x.at(3) * psi_y.at(1);
        f0[p[9]]  = psi_x.at(2) * psi_y.at(1);
        f0[p[10]] = psi_x.at(0) * psi_y.at(3);
        f0[p[11]] = psi_x.at(0) * psi_y.at(2);
        f0[12] = psi_x[2] * psi_y[2];
        f0[13] = psi_x[3] * psi_y[2];
        f0[14] = psi_x[3] * psi_y[3];
        f0[15] = psi_x[2] * psi_y[3];
    }

    // if (whatd & (Fop_D1 | Fop_D2 | Fop_D3)) {
    if (whatd & (Fop_D1 | Fop_D2)) {

        // D1 (dx and dy)
        RN_ f0x(val('.', 0, op_dx)); // derivatives in x direction
        RN_ f0y(val('.', 0, op_dy)); // derivatives in y direction

        R hx = Norme2(edge_x), hy = Norme2(edge_y);

        // dpsi/dx
        std::array<double, 4> d_psi_x{(-13.5 * x * x + 18. * x - 5.5) / hx, 0.5 * (27. * x * x - 18. * x + 2.) / hx,
                                      4.5 * (9. * x * x - 10. * x + 2.) / hx, 4.5 * (-9. * x * x + 8. * x - 1.) / hx};
        // dpsi/dy
        std::array<double, 4> d_psi_y{(-13.5 * y * y + 18. * y - 5.5) / hy, 0.5 * (27. * y * y - 18. * y + 2.) / hy,
                                      4.5 * (9. * y * y - 10. * y + 2.) / hy, 4.5 * (-9. * y * y + 8. * y - 1.) / hy};

        // dx
        f0x[0]  = d_psi_x[0] * psi_y[0];
        f0x[1]  = d_psi_x[1] * psi_y[0];
        f0x[2]  = d_psi_x[1] * psi_y[1];
        f0x[3]  = d_psi_x[0] * psi_y[1];
        f0x[p[4]]  = d_psi_x[2] * psi_y[0];
        f0x[p[5]]  = d_psi_x[3] * psi_y[0];
        f0x[p[6]]  = d_psi_x[1] * psi_y[2];
        f0x[p[7]]  = d_psi_x[1] * psi_y[3];
        f0x[p[8]]  = d_psi_x[3] * psi_y[1];
        f0x[p[9]]  = d_psi_x[2] * psi_y[1];
        f0x[p[10]] = d_psi_x[0] * psi_y[3];
        f0x[p[11]] = d_psi_x[0] * psi_y[2];
        f0x[12] = d_psi_x[2] * psi_y[2];
        f0x[13] = d_psi_x[3] * psi_y[2];
        f0x[14] = d_psi_x[3] * psi_y[3];
        f0x[15] = d_psi_x[2] * psi_y[3];

        // dy
        f0y[0]  = psi_x[0] * d_psi_y[0];
        f0y[1]  = psi_x[1] * d_psi_y[0];
        f0y[2]  = psi_x[1] * d_psi_y[1];
        f0y[3]  = psi_x[0] * d_psi_y[1];
        f0y[p[4]]  = psi_x[2] * d_psi_y[0];
        f0y[p[5]]  = psi_x[3] * d_psi_y[0];
        f0y[p[6]]  = psi_x[1] * d_psi_y[2];
        f0y[p[7]]  = psi_x[1] * d_psi_y[3];
        f0y[p[8]]  = psi_x[3] * d_psi_y[1];
        f0y[p[9]]  = psi_x[2] * d_psi_y[1];
        f0y[p[10]] = psi_x[0] * d_psi_y[3];
        f0y[p[11]] = psi_x[0] * d_psi_y[2];
        f0y[12] = psi_x[2] * d_psi_y[2];
        f0y[13] = psi_x[3] * d_psi_y[2];
        f0y[14] = psi_x[3] * d_psi_y[3];
        f0y[15] = psi_x[2] * d_psi_y[3];

        // std::cout << "x = " << x << "\n";
        // std::cout << "y = " << y << "\n";
        // std::cout << "d_psi_x[0] = " << d_psi_x[0] << "\n";
        // std::cout << "d_psi_y[0] = " << d_psi_y[0] << "\n";
        // std::cout << "f0x[11] = " << f0x[14]*hx << "\n";
        // std::cout << "f0y[11] = " << f0y[14]*hy << "\n";
        // getchar();

        // ddpsi/dxx
        R dd_psi_xx[] = {(-27. * x + 18.) / (hx * hx), (27. * x - 9.) / (hx * hx), (81. * x - 45.) / (hx * hx),
                         (-81 * x + 36.) / (hx * hx)};

        // ddpsi/dyy
        R dd_psi_yy[] = {(-27. * y + 18.) / (hy * hy), (27. * y - 9.) / (hy * hy), (81. * y - 45.) / (hy * hy),
                         (-81. * y + 36.) / (hy * hy)};

        R ddd_psi_xxx[] = {-27. / (hx * hx * hx), 27. / (hx * hx * hx), 81. / (hx * hx * hx), -81. / (hx * hx * hx)};

        R ddd_psi_yyy[] = {-27. / (hy * hy * hy), 27. / (hy * hy * hy), 81. / (hy * hy * hy), -81. / (hy * hy * hy)};

        // std::cout << "x = " << x << "\n";
        // std::cout << "psi_0(x) = " << psi_x[0] << "\n";
        // std::cout << "psi_1(x) = " << psi_x[1] << "\n";
        // std::cout << "psi_2(x) = " << psi_x[2] << "\n";
        // std::cout << "psi_3(x) = " << psi_x[3] << "\n";
        // std::cout << "\n";
        // std::cout << "d_psi_0(x) = " << d_psi_x[0]*hx << "\n";
        // std::cout << "d_psi_1(x) = " << d_psi_x[1]*hx << "\n";
        // std::cout << "d_psi_2(x) = " << d_psi_x[2]*hx << "\n";
        // std::cout << "d_psi_3(x) = " << d_psi_x[3]*hx << "\n";
        // std::cout << "\n";
        // std::cout << "dd_psi_0(x) = " << dd_psi_xx[0]*hx*hx << "\n";
        // std::cout << "dd_psi_1(x) = " << dd_psi_xx[1]*hx*hx << "\n";
        // std::cout << "dd_psi_2(x) = " << dd_psi_xx[2]*hx*hx << "\n";
        // std::cout << "dd_psi_3(x) = " << dd_psi_xx[3]*hx*hx << "\n";
        // std::cout << "\n";
        // std::cout << "ddd_psi_0(x) = " << ddd_psi_xxx[0]*hx*hx*hx << "\n";
        // std::cout << "ddd_psi_1(x) = " << ddd_psi_xxx[1]*hx*hx*hx << "\n";
        // std::cout << "ddd_psi_2(x) = " << ddd_psi_xxx[2]*hx*hx*hx << "\n";
        // std::cout << "ddd_psi_3(x) = " << ddd_psi_xxx[3]*hx*hx*hx << "\n";

        // getchar();

        // dxx, dyy and dxy
        if (whatd & Fop_D2) {

            RN_ f0xx(val('.', 0, op_dxx));
            RN_ f0yy(val('.', 0, op_dyy));
            RN_ f0xy(val('.', 0, op_dxy));

            // dxx
            f0xx[0]  = dd_psi_xx[0] * psi_y[0];
            f0xx[1]  = dd_psi_xx[1] * psi_y[0];
            f0xx[2]  = dd_psi_xx[1] * psi_y[1];
            f0xx[3]  = dd_psi_xx[0] * psi_y[1];
            f0xx[p[4]]  = dd_psi_xx[2] * psi_y[0];
            f0xx[p[5]]  = dd_psi_xx[3] * psi_y[0];
            f0xx[p[6]]  = dd_psi_xx[1] * psi_y[2];
            f0xx[p[7]]  = dd_psi_xx[1] * psi_y[3];
            f0xx[p[8]]  = dd_psi_xx[3] * psi_y[1];
            f0xx[p[9]]  = dd_psi_xx[2] * psi_y[1];
            f0xx[p[10]] = dd_psi_xx[0] * psi_y[3];
            f0xx[p[11]] = dd_psi_xx[0] * psi_y[2];
            f0xx[12] = dd_psi_xx[2] * psi_y[2];
            f0xx[13] = dd_psi_xx[3] * psi_y[2];
            f0xx[14] = dd_psi_xx[3] * psi_y[3];
            f0xx[15] = dd_psi_xx[2] * psi_y[3];

            // dyy
            f0yy[0]  = psi_x[0] * dd_psi_yy[0];
            f0yy[1]  = psi_x[1] * dd_psi_yy[0];
            f0yy[2]  = psi_x[1] * dd_psi_yy[1];
            f0yy[3]  = psi_x[0] * dd_psi_yy[1];
            f0yy[p[4]]  = psi_x[2] * dd_psi_yy[0];
            f0yy[p[5]]  = psi_x[3] * dd_psi_yy[0];
            f0yy[p[6]]  = psi_x[1] * dd_psi_yy[2];
            f0yy[p[7]]  = psi_x[1] * dd_psi_yy[3];
            f0yy[p[8]]  = psi_x[3] * dd_psi_yy[1];
            f0yy[p[9]]  = psi_x[2] * dd_psi_yy[1];
            f0yy[p[10]] = psi_x[0] * dd_psi_yy[3];
            f0yy[p[11]] = psi_x[0] * dd_psi_yy[2];
            f0yy[12] = psi_x[2] * dd_psi_yy[2];
            f0yy[13] = psi_x[3] * dd_psi_yy[2];
            f0yy[14] = psi_x[3] * dd_psi_yy[3];
            f0yy[15] = psi_x[2] * dd_psi_yy[3];

            // dxy
            f0xy[0]  = d_psi_x[0] * d_psi_y[0];
            f0xy[1]  = d_psi_x[1] * d_psi_y[0];
            f0xy[2]  = d_psi_x[1] * d_psi_y[1];
            f0xy[3]  = d_psi_x[0] * d_psi_y[1];
            f0xy[p[4]]  = d_psi_x[2] * d_psi_y[0];
            f0xy[p[5]]  = d_psi_x[3] * d_psi_y[0];
            f0xy[p[6]]  = d_psi_x[1] * d_psi_y[2];
            f0xy[p[7]]  = d_psi_x[1] * d_psi_y[3];
            f0xy[p[8]]  = d_psi_x[3] * d_psi_y[1];
            f0xy[p[9]]  = d_psi_x[2] * d_psi_y[1];
            f0xy[p[10]] = d_psi_x[0] * d_psi_y[3];
            f0xy[p[11]] = d_psi_x[0] * d_psi_y[2];
            f0xy[12] = d_psi_x[2] * d_psi_y[2];
            f0xy[13] = d_psi_x[3] * d_psi_y[2];
            f0xy[14] = d_psi_x[3] * d_psi_y[3];
            f0xy[15] = d_psi_x[2] * d_psi_y[3];
        }

        // // dxxx, dxxy, dxyy, dyyy
        // if (whatd & Fop_D3) {

        //     RN_ f0xxx(val('.', 0, op_dxxx));
        //     RN_ f0xxy(val('.', 0, op_dxxy));
        //     RN_ f0xyy(val('.', 0, op_dxyy));
        //     RN_ f0yyy(val('.', 0, op_dyyy));

        //     // dxxx
        //     f0xxx[0]  = ddd_psi_xxx[0] * psi_y[0];
        //     f0xxx[1]  = ddd_psi_xxx[1] * psi_y[0];
        //     f0xxx[2]  = ddd_psi_xxx[1] * psi_y[1];
        //     f0xxx[3]  = ddd_psi_xxx[0] * psi_y[1];
        //     f0xxx[p[4]]  = ddd_psi_xxx[2] * psi_y[0];
        //     f0xxx[p[5]]  = ddd_psi_xxx[3] * psi_y[0];
        //     f0xxx[p[6]]  = ddd_psi_xxx[1] * psi_y[2];
        //     f0xxx[p[7]]  = ddd_psi_xxx[1] * psi_y[3];
        //     f0xxx[p[8]]  = ddd_psi_xxx[3] * psi_y[1];
        //     f0xxx[p[9]]  = ddd_psi_xxx[2] * psi_y[1];
        //     f0xxx[p[10]] = ddd_psi_xxx[0] * psi_y[3];
        //     f0xxx[p[11]] = ddd_psi_xxx[0] * psi_y[2];
        //     f0xxx[12] = ddd_psi_xxx[2] * psi_y[2];
        //     f0xxx[13] = ddd_psi_xxx[3] * psi_y[2];
        //     f0xxx[14] = ddd_psi_xxx[3] * psi_y[3];
        //     f0xxx[15] = ddd_psi_xxx[2] * psi_y[3];

        //     // dxxy
        //     f0xxy[0]  = dd_psi_xx[0] * d_psi_y[0];
        //     f0xxy[1]  = dd_psi_xx[1] * d_psi_y[0];
        //     f0xxy[2]  = dd_psi_xx[1] * d_psi_y[1];
        //     f0xxy[3]  = dd_psi_xx[0] * d_psi_y[1];
        //     f0xxy[p[4]]  = dd_psi_xx[2] * d_psi_y[0];
        //     f0xxy[p[5]]  = dd_psi_xx[3] * d_psi_y[0];
        //     f0xxy[p[6]]  = dd_psi_xx[1] * d_psi_y[2];
        //     f0xxy[p[7]]  = dd_psi_xx[1] * d_psi_y[3];
        //     f0xxy[p[8]]  = dd_psi_xx[3] * d_psi_y[1];
        //     f0xxy[p[9]]  = dd_psi_xx[2] * d_psi_y[1];
        //     f0xxy[p[10]] = dd_psi_xx[0] * d_psi_y[3];
        //     f0xxy[p[11]] = dd_psi_xx[0] * d_psi_y[2];
        //     f0xxy[12] = dd_psi_xx[2] * d_psi_y[2];
        //     f0xxy[13] = dd_psi_xx[3] * d_psi_y[2];
        //     f0xxy[14] = dd_psi_xx[3] * d_psi_y[3];
        //     f0xxy[15] = dd_psi_xx[2] * d_psi_y[3];

        //     // dxyy
        //     f0xyy[0]  = d_psi_x[0] * dd_psi_yy[0];
        //     f0xyy[1]  = d_psi_x[1] * dd_psi_yy[0];
        //     f0xyy[2]  = d_psi_x[1] * dd_psi_yy[1];
        //     f0xyy[3]  = d_psi_x[0] * dd_psi_yy[1];
        //     f0xyy[p[4]]  = d_psi_x[2] * dd_psi_yy[0];
        //     f0xyy[p[5]]  = d_psi_x[3] * dd_psi_yy[0];
        //     f0xyy[p[6]]  = d_psi_x[1] * dd_psi_yy[2];
        //     f0xyy[p[7]]  = d_psi_x[1] * dd_psi_yy[3];
        //     f0xyy[p[8]]  = d_psi_x[3] * dd_psi_yy[1];
        //     f0xyy[p[9]]  = d_psi_x[2] * dd_psi_yy[1];
        //     f0xyy[p[10]] = d_psi_x[0] * dd_psi_yy[3];
        //     f0xyy[p[11]] = d_psi_x[0] * dd_psi_yy[2];
        //     f0xyy[12] = d_psi_x[2] * dd_psi_yy[2];
        //     f0xyy[13] = d_psi_x[3] * dd_psi_yy[2];
        //     f0xyy[14] = d_psi_x[3] * dd_psi_yy[3];
        //     f0xyy[15] = d_psi_x[2] * dd_psi_yy[3];

        //     // dyyy
        //     f0yyy[0]  = psi_x[0] * ddd_psi_yyy[0];
        //     f0yyy[1]  = psi_x[1] * ddd_psi_yyy[0];
        //     f0yyy[2]  = psi_x[1] * ddd_psi_yyy[1];
        //     f0yyy[3]  = psi_x[0] * ddd_psi_yyy[1];
        //     f0yyy[p[4]]  = psi_x[2] * ddd_psi_yyy[0];
        //     f0yyy[p[5]]  = psi_x[3] * ddd_psi_yyy[0];
        //     f0yyy[p[6]]  = psi_x[1] * ddd_psi_yyy[2];
        //     f0yyy[p[7]]  = psi_x[1] * ddd_psi_yyy[3];
        //     f0yyy[p[8]]  = psi_x[3] * ddd_psi_yyy[1];
        //     f0yyy[p[9]]  = psi_x[2] * ddd_psi_yyy[1];
        //     f0yyy[p[10]] = psi_x[0] * ddd_psi_yyy[3];
        //     f0yyy[p[11]] = psi_x[0] * ddd_psi_yyy[2];
        //     f0yyy[12] = psi_x[2] * ddd_psi_yyy[2];
        //     f0yyy[13] = psi_x[3] * ddd_psi_yyy[2];
        //     f0yyy[14] = psi_x[3] * ddd_psi_yyy[3];
        //     f0yyy[15] = psi_x[2] * ddd_psi_yyy[3];
        // }

    }
}

static TypeOfFE_P3QdcLagrange2d P3Qdc_2d;
GTypeOfFE<MeshQuad2> &P3QdcLagrange2d(P3Qdc_2d);
template <> GTypeOfFE<MeshQuad2> &DataFE<MeshQuad2>::P3dc = P3Qdc_2d;
