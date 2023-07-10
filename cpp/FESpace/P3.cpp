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

// P3 Polynomial basis {1, x, x^2, x^3}
class TypeOfFE_P3Polynomial1d : public GTypeOfFE<Mesh1> {

    typedef Mesh1 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 2;
    static const int ndf = (k + 1);
    static int Data[];
    static double alpha_Pi_h[];

    // dof, dim Im, Data, coefInterp, nbPt, coeff
    TypeOfFE_P3Polynomial1d() : GTypeOfFE<Mesh1>(4, 1, Data, 13, 4, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P3poly;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R1 Pt[4] = {R1(0.), R1(1. / 3), R1(2. / 3), R1(1.)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i] = Pt[i];
        }

        ipj_Pi_h[0] = IPJ(0, 0, 0);
        int j       = 1;
        for (int i = 0; i < 4; ++i) {
            ipj_Pi_h[j++] = IPJ(1, i, 0);
        }
        for (int i = 0; i < 4; ++i) {
            ipj_Pi_h[j++] = IPJ(2, i, 0);
        }
        for (int i = 0; i < 4; ++i) {
            ipj_Pi_h[j++] = IPJ(3, i, 0);
        }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P3Polynomial1d::nbNodeOnItem[4] = {1, 1, 0, 0};
int TypeOfFE_P3Polynomial1d::Data[]                = {
    0, 1, 2, 2, // the support number  of the node of the df
    0, 0, 0, 1, // the number of the df on  the node
    0, 1, 2, 2, // the node of the df
    0, 1, 2, 3, // which are de df on sub FE
    1, 1, 0, 0, // nb node on what
    0,          // for each compontant $j=0,N-1$ it give the sub FE associated
    0,          // begin_dfcomp
    4           // end_dfcomp
};

double TypeOfFE_P3Polynomial1d::alpha_Pi_h[] = {1., -5.5, 9. - 4.5, 1., 9., -22.5, 18., -4.5, -4.5, 13.5, -13.5, 4.5};

void TypeOfFE_P3Polynomial1d::FB(const What_d whatd, const Element &K, const R1 &P, RNMK_ &val) const {
    assert(K[0].X() < K[1].X());
    assert(0 <= P.X() && P.X() <= 1);
    R l[] = {1, P.X(), P.X() * P.X(), P.X() * P.X() * P.X()};

    assert(val.N() >= Element::nv);
    assert(val.M() == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
        f0[0] = l[0];
        f0[1] = l[1];
        f0[2] = l[2];
        f0[3] = l[3];
    }
    if (whatd & Fop_D1) {
        if (whatd & Fop_dx) {
            RN_ f0x(val('.', 0, op_dx));
            f0x[0] = 0;
            f0x[1] = 1. / K.measure();
            f0x[2] = 2 * l[1] * f0x[1];
            f0x[3] = 3 * l[2] * f0x[1];
        }
    }
}

static TypeOfFE_P3Polynomial1d P3_Polynomial_1d;
GTypeOfFE<Mesh1> &P3Polynomial1d(P3_Polynomial_1d);
template <> GTypeOfFE<Mesh1> &DataFE<Mesh1>::P3Poly = P3_Polynomial_1d;

// P3 Lagrange - 2D
class TypeOfFE_P3Lagrange2d : public GTypeOfFE<Mesh2> {
    typedef Mesh2 Mesh;
    typedef typename Mesh::Element Element;

  public:
    static const int k   = 3;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static const int nn[10][3];
    static const int aa[10][3];
    static const int ff[10];
    static const int il[10];
    static const int jl[10];
    static const int kl[10];

    TypeOfFE_P3Lagrange2d() : GTypeOfFE<Mesh2>(10, 1, Data, 16, 10, 0) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P3;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R2 Pt[10] = {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.), R2(0 / 3., 3 / 3.), R2(2 / 3., 1 / 3.),
                                  R2(1 / 3., 2 / 3.), R2(0 / 3., 2 / 3.), R2(0 / 3., 1 / 3.), R2(1 / 3., 0 / 3.),
                                  R2(2 / 3., 0 / 3.), R2(1 / 3., 1 / 3.)};

        int other[10] = {-1, -1, -1, 4, 3, 6, 5, 8, 7, -1};

        int kk = 0;

        for (int i = 0; i < ndf; i++) {
            ipj_Pi_h[kk++] = IPJ(i, i, 0);
            if (other[i] >= 0) {
                ipj_Pi_h[kk++] = IPJ(i, other[i], 0);
            }

            Pt_Pi_h[i] = Pt[i];
        }

        // for (int i = 0; i < 16; i++) {
        //     std::cout << "i = " << i << "\n";
        //     std::cout << "ipj_Pi_h[i] triangle = (" << ipj_Pi_h[i].i << ", " << ipj_Pi_h[i].p << ", " <<
        //     ipj_Pi_h[i].j << ")\n";
        // }

        // getchar();
    }
    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;

    void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const {
        for (int i = 0; i < 16; ++i) {
            v[i] = 1;
        }

        int e0     = K.EdgeOrientation(0);
        int e1     = K.EdgeOrientation(1);
        int e2     = K.EdgeOrientation(2);
        int ooo[6] = {e0, e0, e1, e1, e2, e2};
        int iii[6] = {};
        int jjj[6] = {};

        for (int i = 0; i < 6; ++i) {
            iii[i] = 3 + 2 * i; // si  orient = 1
            jjj[i] = 4 + 2 * i; // si orient = -1
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
int TypeOfFE_P3Lagrange2d::Data[] = {
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6, // the support number  of the node of the df
    0, 0, 0, 0, 1, 0, 1, 0, 1, 0, // the number of the df on  the node
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6, // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // which are de df on sub FE
    1, 1, 1, 0,                   // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0,
    10};

void TypeOfFE_P3Lagrange2d::FB(const What_d whatd, const Element &K, const Rd &PHat, RNMK_ &val) const {

    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1. - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};
    assert(val.N() >= 10);
    assert(val.M() == 1);

    int p[10] = {};

    for (int i = 0; i < 10; ++i) {
        p[i] = i;
    }

    if (K.EdgeOrientation(0) < 0) {
        std::swap(p[3], p[4]); // 3,4
    }

    if (K.EdgeOrientation(1) < 0) {
        std::swap(p[5], p[6]); // 5,6
    }

    if (K.EdgeOrientation(2) < 0) {
        std::swap(p[7], p[8]); // 7,8
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

const int TypeOfFE_P3Lagrange2d::nn[10][3] = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {1, 1, 2}, {1, 2, 2},
                                              {0, 2, 2}, {0, 0, 2}, {0, 0, 1}, {0, 1, 1}, {0, 1, 2}};
const int TypeOfFE_P3Lagrange2d::aa[10][3] = {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 0}, {0, 0, 1},
                                              {0, 0, 1}, {0, 1, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
const int TypeOfFE_P3Lagrange2d::ff[10]    = {6, 6, 6, 2, 2, 2, 2, 2, 2, 1};
const int TypeOfFE_P3Lagrange2d::il[10]    = {3, 0, 0, 0, 0, 1, 2, 2, 1, 1};
const int TypeOfFE_P3Lagrange2d::jl[10]    = {0, 3, 0, 2, 1, 0, 0, 1, 2, 1};
const int TypeOfFE_P3Lagrange2d::kl[10]    = {0, 0, 3, 1, 2, 2, 1, 0, 0, 1};

static TypeOfFE_P3Lagrange2d P3_2d;
GTypeOfFE<Mesh2> &P3Lagrange2d(P3_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P3 = P3_2d;

// P3 Quad – 2D
class TypeOfFE_P3QLagrange2d : public GTypeOfFE<MeshQuad2> {
    typedef MeshQuad2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 3;
    static const int ndf = (k + 1) * (k + 1);
    static int Data[];
    static double alpha_Pi_h[];

    //! Unsure if the initialization list below is correct
    TypeOfFE_P3QLagrange2d() : GTypeOfFE<MeshQuad2>(16, 1, Data, 16, 16, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P3;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R2 Pt[16] = {R2(0., 0.),         R2(1., 0.),         R2(1., 1.),         R2(0., 1.),
                                  R2(1. / 3, 0.),     R2(2. / 3, 0.),     R2(1., 1. / 3),     R2(1., 2. / 3),
                                  R2(2. / 3, 1.),     R2(1. / 3, 1.),     R2(0., 2. / 3),     R2(0., 1. / 3),
                                  R2(1. / 3, 1. / 3), R2(2. / 3, 1. / 3), R2(2. / 3, 2. / 3), R2(1. / 3, 2. / 3)};

        // int numbering[16] = {-1, -1, -1, -1, 5, 4, 7, 6, 9, 8, 11, 10, 15, 14, 13, 12};

        // int kk        = 0;
        // ipj_Pi_h[0]   = IPJ(0, 0, 0);
        // ipj_Pi_h[1] = IPJ(1, 1, 0);
        // ipj_Pi_h[2] = IPJ(2, 2, 0);
        // ipj_Pi_h[3] = IPJ(3, 3, 0);
        // ipj_Pi_h[4] = IPJ(4, 4, 0);
        // ipj_Pi_h[5] = IPJ(4, 5, 0);
        // ipj_Pi_h[6] = IPJ(5, 5, 0);
        // ipj_Pi_h[7] = IPJ(5, 4, 0);
        // ipj_Pi_h[8] = IPJ(6, 6, 0);
        // ipj_Pi_h[9] = IPJ(6, 7, 0);
        // ipj_Pi_h[10] = IPJ(7, 7, 0);
        // ipj_Pi_h[11] = IPJ(7, 6, 0);
        // ipj_Pi_h[12] = IPJ(8, 8, 0);
        // ipj_Pi_h[13] = IPJ(8, 9, 0);
        // ipj_Pi_h[14] = IPJ(9, 9, 0);
        // ipj_Pi_h[15] = IPJ(9, 8, 0);
        // ipj_Pi_h[16] = IPJ(10, 10, 0);
        // ipj_Pi_h[17] = IPJ(10, 11, 0);
        // ipj_Pi_h[18] = IPJ(11, 11, 0);
        // ipj_Pi_h[19] = IPJ(11, 10, 0);
        // ipj_Pi_h[20] = IPJ(12, 12, 0);
        // ipj_Pi_h[21]  = IPJ(12, 13, 0);
        // ipj_Pi_h[22]  = IPJ(12, 14, 0);
        // ipj_Pi_h[23]  = IPJ(12, 15, 0);
        // ipj_Pi_h[24]  = IPJ(13, 15, 0);
        // ipj_Pi_h[25]  = IPJ(13, 14, 0);
        // ipj_Pi_h[26]  = IPJ(13, 13, 0);
        // ipj_Pi_h[27]  = IPJ(13, 12, 0);
        // ipj_Pi_h[28]  = IPJ(14, 12, 0);
        // ipj_Pi_h[29]  = IPJ(14, 13, 0);
        // ipj_Pi_h[30]  = IPJ(14, 14, 0);
        // ipj_Pi_h[31]  = IPJ(14, 15, 0);
        // ipj_Pi_h[32]  = IPJ(15, 15, 0);
        // ipj_Pi_h[33]  = IPJ(15, 14, 0);
        // ipj_Pi_h[34]  = IPJ(15, 13, 0);
        // ipj_Pi_h[35]  = IPJ(15, 12, 0);

        for (int i = 0; i < ndf; ++i) {

            // ipj_Pi_h contains
            // 1st index: number of the DOF (0,1,2,3,..., 15)
            // 2nd index: the number of the DOF that has the basis function corresponding to 1st index evaluate to a
            // non-zero number Since for P3 quad, L_i(P_j) = delta_ij, the 1st and 2nd indices should be the same
            // However, shouldn't this be the case for P3 Lagrange on triangles too? The implementation above doesn't
            // follow this logic.
            ipj_Pi_h[i] = IPJ(i, i, 0);

            // ipj_Pi_h[kk++] = IPJ(i, i, 0);
            // if (numbering[i] >= 0) {
            //     ipj_Pi_h[kk++] = IPJ(i, numbering[i], 0);
            //     //ipj_Pi_h[kk] = IPJ(i, numbering[i+1], 0);
            // }

            Pt_Pi_h[i] = Pt[i]; // DOF points
        }
        // for (int i = 0; i < 28; i++) {
        //     std::cout << "i = " << i << "\n";
        //     std::cout << "ipj_Pi_h[i] quad = (" << ipj_Pi_h[i].i << ", " << ipj_Pi_h[i].p << ", " << ipj_Pi_h[i].j <<
        //     ")\n";
        // }

        // getchar();
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P3QLagrange2d::nbNodeOnItem[4] = {1, 1, 0, 0};

int TypeOfFE_P3QLagrange2d::Data[] = {
    0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8,       // the support number  of the node of the df
    0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3,       // the number of the df on  the node
    0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 8,       // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, // which are de df on sub FE
    1, 1, 1, 0,                                           // nb node on what
    0,                                                    // for each compontant $j=0,N-1$ it give the sub FE associated
    0,                                                    // begin_dfcomp
    16                                                    // end_dfcomp
};

double TypeOfFE_P3QLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

void TypeOfFE_P3QLagrange2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
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
    assert((0. - 1e-12 <= x) && (x <= 1. + 1e-12));
    assert((0. - 1e-12 <= y) && (y <= 1. + 1e-12));

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {

        // phi_l = psi_i^x * psi_j^y

        f0[0]  = psi_x.at(0) * psi_y.at(0);
        f0[1]  = psi_x.at(1) * psi_y.at(0);
        f0[2]  = psi_x.at(1) * psi_y.at(1);
        f0[3]  = psi_x.at(0) * psi_y.at(1);
        f0[4]  = psi_x.at(2) * psi_y.at(0);
        f0[5]  = psi_x.at(3) * psi_y.at(0);
        f0[6]  = psi_x.at(1) * psi_y.at(2);
        f0[7]  = psi_x.at(1) * psi_y.at(3);
        f0[8]  = psi_x.at(3) * psi_y.at(1);
        f0[9]  = psi_x.at(2) * psi_y.at(1);
        f0[10] = psi_x.at(0) * psi_y.at(3);
        f0[11] = psi_x.at(0) * psi_y.at(2);
        f0[12] = psi_x.at(2) * psi_y.at(2);
        f0[13] = psi_x.at(3) * psi_y.at(2);
        f0[14] = psi_x.at(3) * psi_y.at(3);
        f0[15] = psi_x.at(2) * psi_y.at(3);
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
        f0x[0]  = d_psi_x.at(0) * psi_y.at(0);
        f0x[1]  = d_psi_x.at(1) * psi_y.at(0);
        f0x[2]  = d_psi_x.at(1) * psi_y.at(1);
        f0x[3]  = d_psi_x.at(0) * psi_y.at(1);
        f0x[4]  = d_psi_x.at(2) * psi_y.at(0);
        f0x[5]  = d_psi_x.at(3) * psi_y.at(0);
        f0x[6]  = d_psi_x.at(1) * psi_y.at(2);
        f0x[7]  = d_psi_x.at(1) * psi_y.at(3);
        f0x[8]  = d_psi_x.at(3) * psi_y.at(1);
        f0x[9]  = d_psi_x.at(2) * psi_y.at(1);
        f0x[10] = d_psi_x.at(0) * psi_y.at(3);
        f0x[11] = d_psi_x.at(0) * psi_y.at(2);
        f0x[12] = d_psi_x.at(2) * psi_y.at(2);
        f0x[13] = d_psi_x.at(3) * psi_y.at(2);
        f0x[14] = d_psi_x.at(3) * psi_y.at(3);
        f0x[15] = d_psi_x.at(2) * psi_y.at(3);

        // dy
        f0y[0]  = psi_x.at(0) * d_psi_y.at(0);
        f0y[1]  = psi_x.at(1) * d_psi_y.at(0);
        f0y[2]  = psi_x.at(1) * d_psi_y.at(1);
        f0y[3]  = psi_x.at(0) * d_psi_y.at(1);
        f0y[4]  = psi_x.at(2) * d_psi_y.at(0);
        f0y[5]  = psi_x.at(3) * d_psi_y.at(0);
        f0y[6]  = psi_x.at(1) * d_psi_y.at(2);
        f0y[7]  = psi_x.at(1) * d_psi_y.at(3);
        f0y[8]  = psi_x.at(3) * d_psi_y.at(1);
        f0y[9]  = psi_x.at(2) * d_psi_y.at(1);
        f0y[10] = psi_x.at(0) * d_psi_y.at(3);
        f0y[11] = psi_x.at(0) * d_psi_y.at(2);
        f0y[12] = psi_x.at(2) * d_psi_y.at(2);
        f0y[13] = psi_x.at(3) * d_psi_y.at(2);
        f0y[14] = psi_x.at(3) * d_psi_y.at(3);
        f0y[15] = psi_x.at(2) * d_psi_y.at(3);

        // std::cout << "x = " << x << "\n";
        // std::cout << "y = " << y << "\n";
        // std::cout << "d_psi_x[0] = " << d_psi_x[0] << "\n";
        // std::cout << "d_psi_y[0] = " << d_psi_y[0] << "\n";
        // std::cout << "f0x[11] = " << f0x[14]*hx << "\n";
        // std::cout << "f0y[11] = " << f0y[14]*hy << "\n";
        // getchar();

        // // ddpsi/dxx
        // R dd_psi_xx[] = {(-27. * x + 18.) / (hx * hx), (27. * x - 9.) / (hx * hx), (81. * x - 45.) / (hx * hx),
        //                  (-81 * x + 36.) / (hx * hx)};

        // // ddpsi/dyy
        // R dd_psi_yy[] = {(-27. * y + 18.) / (hy * hy), (27. * y - 9.) / (hy * hy), (81. * y - 45.) / (hy * hy),
        //                  (-81. * y + 36.) / (hy * hy)};

        // R ddd_psi_xxx[] = {-27. / (hx * hx * hx), 27. / (hx * hx * hx), 81. / (hx * hx * hx), -81. / (hx * hx * hx)};

        // R ddd_psi_yyy[] = {-27. / (hy * hy * hy), 27. / (hy * hy * hy), 81. / (hy * hy * hy), -81. / (hy * hy * hy)};

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

        // // dxx, dyy and dxy
        // if (whatd & Fop_D2) {

        //     RN_ f0xx(val('.', 0, op_dxx));
        //     RN_ f0yy(val('.', 0, op_dyy));
        //     RN_ f0xy(val('.', 0, op_dxy));

        //     // dxx
        //     f0xx[0]  = dd_psi_xx[0] * psi_y[0];
        //     f0xx[1]  = dd_psi_xx[1] * psi_y[0];
        //     f0xx[2]  = dd_psi_xx[1] * psi_y[1];
        //     f0xx[3]  = dd_psi_xx[0] * psi_y[1];
        //     f0xx[4]  = dd_psi_xx[2] * psi_y[0];
        //     f0xx[5]  = dd_psi_xx[3] * psi_y[0];
        //     f0xx[6]  = dd_psi_xx[1] * psi_y[2];
        //     f0xx[7]  = dd_psi_xx[1] * psi_y[3];
        //     f0xx[8]  = dd_psi_xx[3] * psi_y[1];
        //     f0xx[9]  = dd_psi_xx[2] * psi_y[1];
        //     f0xx[10] = dd_psi_xx[0] * psi_y[3];
        //     f0xx[11] = dd_psi_xx[0] * psi_y[2];
        //     f0xx[12] = dd_psi_xx[2] * psi_y[2];
        //     f0xx[13] = dd_psi_xx[3] * psi_y[2];
        //     f0xx[14] = dd_psi_xx[3] * psi_y[3];
        //     f0xx[15] = dd_psi_xx[2] * psi_y[3];

        //     // dyy
        //     f0yy[0]  = psi_x[0] * dd_psi_yy[0];
        //     f0yy[1]  = psi_x[1] * dd_psi_yy[0];
        //     f0yy[2]  = psi_x[1] * dd_psi_yy[1];
        //     f0yy[3]  = psi_x[0] * dd_psi_yy[1];
        //     f0yy[4]  = psi_x[2] * dd_psi_yy[0];
        //     f0yy[5]  = psi_x[3] * dd_psi_yy[0];
        //     f0yy[6]  = psi_x[1] * dd_psi_yy[2];
        //     f0yy[7]  = psi_x[1] * dd_psi_yy[3];
        //     f0yy[8]  = psi_x[3] * dd_psi_yy[1];
        //     f0yy[9]  = psi_x[2] * dd_psi_yy[1];
        //     f0yy[10] = psi_x[0] * dd_psi_yy[3];
        //     f0yy[11] = psi_x[0] * dd_psi_yy[2];
        //     f0yy[12] = psi_x[2] * dd_psi_yy[2];
        //     f0yy[13] = psi_x[3] * dd_psi_yy[2];
        //     f0yy[14] = psi_x[3] * dd_psi_yy[3];
        //     f0yy[15] = psi_x[2] * dd_psi_yy[3];

        //     // dxy
        //     f0xy[0]  = d_psi_x[0] * d_psi_y[0];
        //     f0xy[1]  = d_psi_x[1] * d_psi_y[0];
        //     f0xy[2]  = d_psi_x[1] * d_psi_y[1];
        //     f0xy[3]  = d_psi_x[0] * d_psi_y[1];
        //     f0xy[4]  = d_psi_x[2] * d_psi_y[0];
        //     f0xy[5]  = d_psi_x[3] * d_psi_y[0];
        //     f0xy[6]  = d_psi_x[1] * d_psi_y[2];
        //     f0xy[7]  = d_psi_x[1] * d_psi_y[3];
        //     f0xy[8]  = d_psi_x[3] * d_psi_y[1];
        //     f0xy[9]  = d_psi_x[2] * d_psi_y[1];
        //     f0xy[10] = d_psi_x[0] * d_psi_y[3];
        //     f0xy[11] = d_psi_x[0] * d_psi_y[2];
        //     f0xy[12] = d_psi_x[2] * d_psi_y[2];
        //     f0xy[13] = d_psi_x[3] * d_psi_y[2];
        //     f0xy[14] = d_psi_x[3] * d_psi_y[3];
        //     f0xy[15] = d_psi_x[2] * d_psi_y[3];
        // }

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
        //     f0xxx[4]  = ddd_psi_xxx[2] * psi_y[0];
        //     f0xxx[5]  = ddd_psi_xxx[3] * psi_y[0];
        //     f0xxx[6]  = ddd_psi_xxx[1] * psi_y[2];
        //     f0xxx[7]  = ddd_psi_xxx[1] * psi_y[3];
        //     f0xxx[8]  = ddd_psi_xxx[3] * psi_y[1];
        //     f0xxx[9]  = ddd_psi_xxx[2] * psi_y[1];
        //     f0xxx[10] = ddd_psi_xxx[0] * psi_y[3];
        //     f0xxx[11] = ddd_psi_xxx[0] * psi_y[2];
        //     f0xxx[12] = ddd_psi_xxx[2] * psi_y[2];
        //     f0xxx[13] = ddd_psi_xxx[3] * psi_y[2];
        //     f0xxx[14] = ddd_psi_xxx[3] * psi_y[3];
        //     f0xxx[15] = ddd_psi_xxx[2] * psi_y[3];

        //     // dxxy
        //     f0xxy[0]  = dd_psi_xx[0] * d_psi_y[0];
        //     f0xxy[1]  = dd_psi_xx[1] * d_psi_y[0];
        //     f0xxy[2]  = dd_psi_xx[1] * d_psi_y[1];
        //     f0xxy[3]  = dd_psi_xx[0] * d_psi_y[1];
        //     f0xxy[4]  = dd_psi_xx[2] * d_psi_y[0];
        //     f0xxy[5]  = dd_psi_xx[3] * d_psi_y[0];
        //     f0xxy[6]  = dd_psi_xx[1] * d_psi_y[2];
        //     f0xxy[7]  = dd_psi_xx[1] * d_psi_y[3];
        //     f0xxy[8]  = dd_psi_xx[3] * d_psi_y[1];
        //     f0xxy[9]  = dd_psi_xx[2] * d_psi_y[1];
        //     f0xxy[10] = dd_psi_xx[0] * d_psi_y[3];
        //     f0xxy[11] = dd_psi_xx[0] * d_psi_y[2];
        //     f0xxy[12] = dd_psi_xx[2] * d_psi_y[2];
        //     f0xxy[13] = dd_psi_xx[3] * d_psi_y[2];
        //     f0xxy[14] = dd_psi_xx[3] * d_psi_y[3];
        //     f0xxy[15] = dd_psi_xx[2] * d_psi_y[3];

        //     // dxyy
        //     f0xyy[0]  = d_psi_x[0] * dd_psi_yy[0];
        //     f0xyy[1]  = d_psi_x[1] * dd_psi_yy[0];
        //     f0xyy[2]  = d_psi_x[1] * dd_psi_yy[1];
        //     f0xyy[3]  = d_psi_x[0] * dd_psi_yy[1];
        //     f0xyy[4]  = d_psi_x[2] * dd_psi_yy[0];
        //     f0xyy[5]  = d_psi_x[3] * dd_psi_yy[0];
        //     f0xyy[6]  = d_psi_x[1] * dd_psi_yy[2];
        //     f0xyy[7]  = d_psi_x[1] * dd_psi_yy[3];
        //     f0xyy[8]  = d_psi_x[3] * dd_psi_yy[1];
        //     f0xyy[9]  = d_psi_x[2] * dd_psi_yy[1];
        //     f0xyy[10] = d_psi_x[0] * dd_psi_yy[3];
        //     f0xyy[11] = d_psi_x[0] * dd_psi_yy[2];
        //     f0xyy[12] = d_psi_x[2] * dd_psi_yy[2];
        //     f0xyy[13] = d_psi_x[3] * dd_psi_yy[2];
        //     f0xyy[14] = d_psi_x[3] * dd_psi_yy[3];
        //     f0xyy[15] = d_psi_x[2] * dd_psi_yy[3];

        //     // dyyy
        //     f0yyy[0]  = psi_x[0] * ddd_psi_yyy[0];
        //     f0yyy[1]  = psi_x[1] * ddd_psi_yyy[0];
        //     f0yyy[2]  = psi_x[1] * ddd_psi_yyy[1];
        //     f0yyy[3]  = psi_x[0] * ddd_psi_yyy[1];
        //     f0yyy[4]  = psi_x[2] * ddd_psi_yyy[0];
        //     f0yyy[5]  = psi_x[3] * ddd_psi_yyy[0];
        //     f0yyy[6]  = psi_x[1] * ddd_psi_yyy[2];
        //     f0yyy[7]  = psi_x[1] * ddd_psi_yyy[3];
        //     f0yyy[8]  = psi_x[3] * ddd_psi_yyy[1];
        //     f0yyy[9]  = psi_x[2] * ddd_psi_yyy[1];
        //     f0yyy[10] = psi_x[0] * ddd_psi_yyy[3];
        //     f0yyy[11] = psi_x[0] * ddd_psi_yyy[2];
        //     f0yyy[12] = psi_x[2] * ddd_psi_yyy[2];
        //     f0yyy[13] = psi_x[3] * ddd_psi_yyy[2];
        //     f0yyy[14] = psi_x[3] * ddd_psi_yyy[3];
        //     f0yyy[15] = psi_x[2] * ddd_psi_yyy[3];
        // }
    }
}

static TypeOfFE_P3QLagrange2d P3Q_2d;
GTypeOfFE<MeshQuad2> &P3QLagrange2d(P3Q_2d);
template <> GTypeOfFE<MeshQuad2> &DataFE<MeshQuad2>::P3 = P3Q_2d;
