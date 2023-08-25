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

// P1 Polynomial basis {1, x}
class TypeOfFE_P2Polynomial1d : public GTypeOfFE<Mesh1> {

    typedef Mesh1 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 2;
    static const int ndf = (k + 1);
    static int Data[];
    static double alpha_Pi_h[];

    // dof, dim Im, Data, coefInterp, nbPt, coeff
    TypeOfFE_P2Polynomial1d() : GTypeOfFE<Mesh1>(3, 1, Data, 7, 3, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P2poly;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R1 Pt[3] = {R1(0.), R1(1. / 2), R1(1.)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i] = Pt[i];
        }

        ipj_Pi_h[0] = IPJ(0, 0, 0);
        for (int i = 0; i < 3; ++i) {
            ipj_Pi_h[i + 1] = IPJ(1, i, 0);
        }
        for (int i = 0; i < 3; ++i) {
            ipj_Pi_h[i + 4] = IPJ(2, i, 0);
        }

    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P2Polynomial1d::nbNodeOnItem[4] = {1, 1, 0, 0};
int TypeOfFE_P2Polynomial1d::Data[]                = {
    0, 1, 2,    // the support number  of the node of the df
    0, 0, 0,    // the number of the df on  the node
    0, 1, 2,    // the node of the df
    0, 1, 2,    // which are de df on sub FE
    1, 1, 0, 0, // nb node on what
    0,          // for each compontant $j=0,N-1$ it give the sub FE associated
    0,          // begin_dfcomp
    3           // end_dfcomp
};

double TypeOfFE_P2Polynomial1d::alpha_Pi_h[] = {1., -3., 4., -1., 2., -4., 2.};

void TypeOfFE_P2Polynomial1d::FB(const What_d whatd, const Element &K, const R1 &P, RNMK_ &val) const {
    assert(K[0].X() < K[1].X());
    assert(0 <= P.X() && P.X() <= 1);
    R l[] = {1, P.X(), P.X() * P.X()};

    assert(val.N() >= Element::nv);
    assert(val.M() == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
        f0[0] = l[0];
        f0[1] = l[1];
        f0[2] = l[2];
    }
    if (whatd & Fop_D1) {
        if (whatd & Fop_dx) {
            RN_ f0x(val('.', 0, op_dx));
            f0x[0] = 0;
            f0x[1] = 1. / K.measure();
            f0x[2] = 2 * l[1] * f0x[1];
        }
    }
}

static TypeOfFE_P2Polynomial1d P2_Polynomial_1d;
GTypeOfFE<Mesh1> &P2Polynomial1d(P2_Polynomial_1d);
template <> GTypeOfFE<Mesh1> &DataFE<Mesh1>::P2Poly = P2_Polynomial_1d;

// P2 - 2D
class TypeOfFE_P2Lagrange2d : public GTypeOfFE<Mesh2> {

    typedef Mesh2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 2;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double alpha_Pi_h[];

    TypeOfFE_P2Lagrange2d() : GTypeOfFE<Mesh2>(6, 1, Data, 6, 6, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P2;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R2 Pt[6] = {R2(0., 0.), R2(1., 0.), R2(0., 1.), R2(1. / 2, 1. / 2), R2(0., 1. / 2), R2(1. / 2, 0)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i]  = Pt[i];
            ipj_Pi_h[i] = IPJ(i, i, 0);
        }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P2Lagrange2d::nbNodeOnItem[4] = {1, 1, 0, 0};

int TypeOfFE_P2Lagrange2d::Data[] = {
    0, 1, 2, 3, 4, 5, // the support number  of the node of the df
    0, 0, 0, 0, 0, 0, // the number of the df on  the node
    0, 1, 2, 3, 4, 5, // the node of the df
    0, 1, 2, 3, 4, 5, // which are de df on sub FE
    1, 1, 0, 0,       // nb node on what
    0,                // for each compontant $j=0,N-1$ it give the sub FE associated
    0,                // begin_dfcomp
    6                 // end_dfcomp
};

double TypeOfFE_P2Lagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1.};

void TypeOfFE_P2Lagrange2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
    R l[] = {1. - P.sum(), P.x, P.y};

    assert(val.N() >= E::nv + E::ne);
    assert(val.M() == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
        int k = 0;
        for (int i = 0; i < E::nv; ++i)
            f0[k++] = l[i] * (2 * l[i] - 1.);
        for (int i = 0; i < E::ne; ++i)
            f0[k++] = 4. * l[E::nvedge[i][0]] * l[E::nvedge[i][1]];
    }

    if (whatd & (Fop_D1 | Fop_D2)) {
        R2 Dl[3];
        R l4[3] = {(4 * l[0] - 1), (4 * l[1] - 1), (4 * l[2] - 1)};

        K.Gradlambda(Dl);
        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        int k = 0;
        for (int i = 0; i < E::nv; ++i, ++k) {
            f0x[k] = Dl[i].x * l4[i];
            f0y[k] = Dl[i].y * l4[i];
        }
        for (int i = 0; i < E::ne; ++i, ++k) {
            int i0 = E::nvedge[i][0], i1 = E::nvedge[i][1];
            f0x[k] = 4 * (Dl[i1].x * l[i0] + Dl[i0].x * l[i1]);
            f0y[k] = 4 * (Dl[i1].y * l[i0] + Dl[i0].y * l[i1]);
        }
        assert(k == 6);

        if (whatd & Fop_D2) {
            RN_ f0xx(val('.', 0, op_dxx));
            RN_ f0yy(val('.', 0, op_dyy));
            RN_ f0xy(val('.', 0, op_dxy));

            k = 0;
            for (int i = 0; i < E::nv; ++i, ++k) {
                f0xx[k] = 4. * Dl[i].x * Dl[i].x;
                f0yy[k] = 4. * Dl[i].y * Dl[i].y;
                f0xy[k] = 4. * Dl[i].x * Dl[i].y;
            }
            for (int i = 0; i < E::ne; ++i, ++k) {
                int i0 = E::nvedge[i][0], i1 = E::nvedge[i][1];
                f0xx[k] = 8. * Dl[i0].x * Dl[i1].x;
                f0yy[k] = 8. * Dl[i0].y * Dl[i1].y;
                f0xy[k] = 4. * (Dl[i0].x * Dl[i1].y + Dl[i1].x * Dl[i0].y);
            }
        }
    }
}

static TypeOfFE_P2Lagrange2d P2_2d;
GTypeOfFE<Mesh2> &P2Lagrange2d(P2_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P2 = P2_2d;

// P2 Quad - 2D
class TypeOfFE_P2QLagrange2d : public GTypeOfFE<MeshQuad2> {
    typedef MeshQuad2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 2;
    static const int ndf = (k + 1) * (k + 1);
    static int Data[];
    static double alpha_Pi_h[];

    TypeOfFE_P2QLagrange2d() : GTypeOfFE<MeshQuad2>(9, 1, Data, 9, 9, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P2;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R2 Pt[9] = {R2(0., 0.),  R2(1., 0.),  R2(1., 1.),  R2(0., 1.),  R2(0.5, 0.),
                                 R2(1., 0.5), R2(0.5, 1.), R2(0., 0.5), R2(0.5, 0.5)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i]  = Pt[i];
            ipj_Pi_h[i] = IPJ(i, i, 0);
        }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P2QLagrange2d::nbNodeOnItem[4] = {1, 1, 1, 0};

int TypeOfFE_P2QLagrange2d::Data[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, // the support number  of the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, // the number of the df on  the node
    0, 1, 2, 3, 4, 5, 6, 7, 8, // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, // which are de df on sub FE
    1, 1, 1, 0,                // nb node on what
    0,                         // for each compontant $j=0,N-1$ it give the sub FE associated
    0,                         // begin_dfcomp
    9                          // end_dfcomp
};

double TypeOfFE_P2QLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1., 1., 1., 1.};

void TypeOfFE_P2QLagrange2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
    R psi_x[] = {1. - 3 * P.x + 2 * P.x * P.x, P.x * (2 * P.x - 1), 4 * P.x * (1 - P.x)};
    R psi_y[] = {1. - 3 * P.y + 2 * P.y * P.y, P.y * (2 * P.y - 1), 4 * P.y * (1 - P.y)};

    R2 edge_x(K[0], K[1]); // [hx, 0]
    R2 edge_y(K[0], K[3]); // [0, hy]

    R x = P.x, y = P.y;

    assert(val.N() >= E::nv + E::ne);
    assert(val.M() == 1);

    // Assert that points are on the reference element
    assert((0. - 1e-12 <= x) && (x <= 1. + 1e-12));
    assert((0. - 1e-12 <= y) && (y <= 1. + 1e-12));

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {

        // phi_l = psi_i^x * psi_j^y

        f0[0] = psi_x[0] * psi_y[0];
        f0[1] = psi_x[1] * psi_y[0];
        f0[2] = psi_x[1] * psi_y[1];
        f0[3] = psi_x[0] * psi_y[1];
        f0[4] = psi_x[2] * psi_y[0];
        f0[5] = psi_x[1] * psi_y[2];
        f0[6] = psi_x[2] * psi_y[1];
        f0[7] = psi_x[0] * psi_y[2];
        f0[8] = psi_x[2] * psi_y[2];
    }

    if (whatd & (Fop_D1 | Fop_D2)) {
        RN_ f0x(val('.', 0, op_dx)); // derivatives in x direction
        RN_ f0y(val('.', 0, op_dy)); // derivatives in y direction
        R hx = Norme2(edge_x), hy = Norme2(edge_y);
        R d_psi_x[] = {(-3 + 4 * P.x) / hx, (4 * P.x - 1) / hx, (4 - 8 * P.x) / hx};
        R d_psi_y[] = {(-3 + 4 * P.y) / hy, (4 * P.y - 1) / hy, (4 - 8 * P.y) / hy};

        // d(phi)/dx = d(psi^x)/dx * psi^y
        f0x[0] = d_psi_x[0] * psi_y[0];
        f0x[1] = d_psi_x[1] * psi_y[0];
        f0x[2] = d_psi_x[1] * psi_y[1];
        f0x[3] = d_psi_x[0] * psi_y[1];
        f0x[4] = d_psi_x[2] * psi_y[0];
        f0x[5] = d_psi_x[1] * psi_y[2];
        f0x[6] = d_psi_x[2] * psi_y[1];
        f0x[7] = d_psi_x[0] * psi_y[2];
        f0x[8] = d_psi_x[2] * psi_y[2];

        // d(phi)/dy = psi^x * d(psi^y)/dy
        f0y[0] = psi_x[0] * d_psi_y[0];
        f0y[1] = psi_x[1] * d_psi_y[0];
        f0y[2] = psi_x[1] * d_psi_y[1];
        f0y[3] = psi_x[0] * d_psi_y[1];
        f0y[4] = psi_x[2] * d_psi_y[0];
        f0y[5] = psi_x[1] * d_psi_y[2];
        f0y[6] = psi_x[2] * d_psi_y[1];
        f0y[7] = psi_x[0] * d_psi_y[2];
        f0y[8] = psi_x[2] * d_psi_y[2];


        if (whatd & Fop_D2) {
            RN_ f0xx(val('.', 0, op_dxx));
            RN_ f0yy(val('.', 0, op_dyy));
            RN_ f0xy(val('.', 0, op_dxy));

            R dd_psi_xx[] = {4. / (hx * hx), 4. / (hx * hx), -8. / (hx * hx)};
            R dd_psi_yy[] = {4. / (hy * hy), 4. / (hy * hy), -8. / (hy * hy)};

            f0xx[0] = dd_psi_xx[0] * psi_y[0];
            f0xx[1] = dd_psi_xx[1] * psi_y[0];
            f0xx[2] = dd_psi_xx[1] * psi_y[1];
            f0xx[3] = dd_psi_xx[0] * psi_y[1];
            f0xx[4] = dd_psi_xx[2] * psi_y[0];
            f0xx[5] = dd_psi_xx[1] * psi_y[2];
            f0xx[6] = dd_psi_xx[2] * psi_y[1];
            f0xx[7] = dd_psi_xx[0] * psi_y[2];
            f0xx[8] = dd_psi_xx[2] * psi_y[2];

            f0yy[0] = psi_x[0] * dd_psi_yy[0];
            f0yy[1] = psi_x[1] * dd_psi_yy[0];
            f0yy[2] = psi_x[1] * dd_psi_yy[1];
            f0yy[3] = psi_x[0] * dd_psi_yy[1];
            f0yy[4] = psi_x[2] * dd_psi_yy[0];
            f0yy[5] = psi_x[1] * dd_psi_yy[2];
            f0yy[6] = psi_x[2] * dd_psi_yy[1];
            f0yy[7] = psi_x[0] * dd_psi_yy[2];
            f0yy[8] = psi_x[2] * dd_psi_yy[2];

            f0xy[0] = d_psi_x[0] * d_psi_y[0];
            f0xy[1] = d_psi_x[1] * d_psi_y[0];
            f0xy[2] = d_psi_x[1] * d_psi_y[1];
            f0xy[3] = d_psi_x[0] * d_psi_y[1];
            f0xy[4] = d_psi_x[2] * d_psi_y[0];
            f0xy[5] = d_psi_x[1] * d_psi_y[2];
            f0xy[6] = d_psi_x[2] * d_psi_y[1];
            f0xy[7] = d_psi_x[0] * d_psi_y[2];
            f0xy[8] = d_psi_x[2] * d_psi_y[2];

        }
    }
}

static TypeOfFE_P2QLagrange2d P2Q_2d;
GTypeOfFE<MeshQuad2> &P2QLagrange2d(P2Q_2d);
template <> GTypeOfFE<MeshQuad2> &DataFE<MeshQuad2>::P2 = P2Q_2d;

// P2 - 3D
class TypeOfFE_P2Lagrange3d : public GTypeOfFE<Mesh3> {

    typedef Mesh3 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 2;
    static const int ndf = 10; //(k + 2) * (k + 1) / 2;
    static int Data[];
    static double alpha_Pi_h[];

    TypeOfFE_P2Lagrange3d() : GTypeOfFE<Mesh3>(10, 1, Data, 10, 10, alpha_Pi_h) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P2;
        GTypeOfFE<Mesh>::polynomialOrder = k;

        static const R3 Pt[10] = {R3(0., 0., 0.),        R3(1., 0., 0.),         R3(0., 1., 0.),
                                  R3(0., 0., 1.),        R3(1. / 2, 0., 0.),     R3(0., 1. / 2, 0.),
                                  R3(0., 0., 1. / 2),    R3(1. / 2, 1. / 2, 0.), R3(1. / 2, 0., 1. / 2),
                                  R3(0., 1. / 2, 1. / 2)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i]  = Pt[i];
            ipj_Pi_h[i] = IPJ(i, i, 0);
        }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P2Lagrange3d::nbNodeOnItem[4] = {1, 1, 0, 0};
int TypeOfFE_P2Lagrange3d::Data[]                = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // the support number  of the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // the number of the df on  the node
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // which are de df on sub FE
    1, 1, 0, 0,                   // nb node on what
    0,                            // for each compontant $j=0,N-1$ it give the sub FE associated
    0,                            // begin_dfcomp
    10                            // end_dfcomp
};

double TypeOfFE_P2Lagrange3d::alpha_Pi_h[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

void TypeOfFE_P2Lagrange3d::FB(const What_d whatd, const Element &K, const R3 &P, RNMK_ &val) const {
    R l[] = {1. - P.sum(), P.x, P.y, P.z};

    assert(val.N() >= E::nv + E::ne);
    assert(val.M() == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) {
        int k = 0;
        for (int i = 0; i < E::nv; ++i)
            f0[k++] = l[i] * (2 * l[i] - 1.);
        for (int i = 0; i < E::ne; ++i)
            f0[k++] = 4. * l[E::nvedge[i][0]] * l[E::nvedge[i][1]];
    }

    if (whatd & (Fop_D1 | Fop_D2)) {

        R3 Dl[4];
        R l4[4] = {(4 * l[0] - 1), (4 * l[1] - 1), (4 * l[2] - 1), (4 * l[3] - 1)};
        K.Gradlambda(Dl);

        if (whatd & (Fop_D1)) {
            RN_ f0x(val('.', 0, op_dx));
            RN_ f0y(val('.', 0, op_dy));
            RN_ f0z(val('.', 0, op_dz));

            int k = 0;
            for (int i = 0; i < E::nv; ++i, ++k) {
                f0x[k] = Dl[i].x * l4[i];
                f0y[k] = Dl[i].y * l4[i];
                f0z[k] = Dl[i].z * l4[i];
            }

            for (int i = 0; i < E::ne; ++i, ++k) {
                int i0 = E::nvedge[i][0], i1 = E::nvedge[i][1];
                f0x[k] = 4 * (Dl[i1].x * l[i0] + Dl[i0].x * l[i1]);
                f0y[k] = 4 * (Dl[i1].y * l[i0] + Dl[i0].y * l[i1]);
                f0z[k] = 4 * (Dl[i1].z * l[i0] + Dl[i0].z * l[i1]);
            }

            assert(k == 10);
        }

        if (whatd & Fop_D2) {
            RN_ f0xx(val('.', 0, op_dxx));
            RN_ f0yy(val('.', 0, op_dyy));
            RN_ f0zz(val('.', 0, op_dzz));
            RN_ f0xy(val('.', 0, op_dxy));
            RN_ f0xz(val('.', 0, op_dxz));
            RN_ f0yz(val('.', 0, op_dyz));

            int k = 0;
            for (int i = 0; i < E::nv; ++i, ++k) {
                f0xx[k] = 4. * Dl[i].x * Dl[i].x;
                f0yy[k] = 4. * Dl[i].y * Dl[i].y;
                f0zz[k] = 4. * Dl[i].z * Dl[i].z;
                f0xy[k] = 4. * Dl[i].x * Dl[i].y;
                f0xz[k] = 4. * Dl[i].x * Dl[i].z;
                f0yz[k] = 4. * Dl[i].y * Dl[i].z;
            }
            for (int i = 0; i < E::ne; ++i, ++k) {
                int i0 = E::nvedge[i][0], i1 = E::nvedge[i][1];
                f0xx[k] = 8. * Dl[i0].x * Dl[i1].x;
                f0yy[k] = 8. * Dl[i0].y * Dl[i1].y;
                f0zz[k] = 8. * Dl[i0].z * Dl[i1].z;
                f0xy[k] = 4. * (Dl[i0].x * Dl[i1].y + Dl[i1].x * Dl[i0].y);
                f0xz[k] = 4. * (Dl[i0].x * Dl[i1].z + Dl[i1].x * Dl[i0].z);
                f0yz[k] = 4. * (Dl[i0].y * Dl[i1].z + Dl[i1].y * Dl[i0].z);
            }
            assert(k == 10);
        }
    }
}

static TypeOfFE_P2Lagrange3d P2_3d;
GTypeOfFE<Mesh3> &P2Lagrange3d(P2_3d);
template <> GTypeOfFE<Mesh3> &DataFE<Mesh3>::P2 = P2_3d;
