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

class TypeOfFE_BDM1_2d : public GTypeOfFE<Mesh2> {

    typedef Mesh2 Mesh;
    const QuadratureFormular1d &QFE;

  public:
    static int Data[];
    // static double Pi_h_coef[]; // [doesnt seem to be used]
    TypeOfFE_BDM1_2d()
        : GTypeOfFE<Mesh2>(6, 2, Data,
                           2 * 2 * 3 * 2, // nb coef mat interpole
                           3 * 2          // nb P interpolation 3-> nb of edges
                           ),
          QFE(QF_GaussLegendre2) // quadrature formula with 2 points
    {

        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::BDM1;
        GTypeOfFE<Mesh>::polynomialOrder = 1;

        Triangle2 TriangleHat;
        Vertex2 verticesHat[3];
        (R2 &)verticesHat[0] = R2::KHat[0];
        (R2 &)verticesHat[1] = R2::KHat[1];
        (R2 &)verticesHat[2] = R2::KHat[2];
        int iv[3]            = {0, 1, 2};
        TriangleHat.set(verticesHat, iv, 0);

        int kkk = 0, i = 0;

        for (int e = 0; e < 3; ++e) {
            for (int p = 0; p < QFE.n; ++p) {
                R2 A(TriangleHat[Element::nvedge[e][0]]);
                R2 B(TriangleHat[Element::nvedge[e][1]]);

                ipj_Pi_h[kkk++] = IPJ(2 * e, i, 0);
                ipj_Pi_h[kkk++] = IPJ(2 * e, i, 1);
                ipj_Pi_h[kkk++] = IPJ(2 * e + 1, i, 0);
                ipj_Pi_h[kkk++] = IPJ(2 * e + 1, i, 1);

                Pt_Pi_h[i++] = B * (QFE[p].X()) + A * (1. - QFE[p].X()); // X=0 => A  X=1 => B;
            }
        }

        ffassert(kkk == this->ipj_Pi_h.N());
        ffassert(i == this->Pt_Pi_h.N());
    }

    void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K,
                       KN_<double> &v) const { // compute the coef of interpolation ...
        const Element &T = K.T;
        int k            = 0;
        double s         = 1. / sqrt(T.measure());    // gives h^-3 for Taylor-Green vortex problem
        //double s         = 1.;  // gives h^-1.15 for Taylor-Green vortex problem
        //double s = sqrt(T.measure());   // gives decreasing condition numbers

        for (int i = 0; i < 3; i++) {
            R2 E(-T.Edge(i).perp());
            R eOrientation = T.EdgeOrientation(i);

            for (int p = 0; p < QFE.n; ++p) {
                R l0 = QFE[p].X(), l1 = 1 - QFE[p].X();
                R p0   = eOrientation;   // poly othogonaux to \lambda_1
                R p1   = -3 * (l0 - l1); // poly othogonaux to \lambda_0
                R cc0  = p0 * QFE[p].a;  //
                R cc1  = p1 * QFE[p].a;  //
                v[k++] = cc0 * E.x * s;
                v[k++] = cc0 * E.y * s;
                v[k++] = cc1 * E.x * s;
                v[k++] = cc1 * E.y * s;
            }
        }

        assert(k == this->ipj_Pi_h.N());
    }

    // void FB(const What_d, const Mesh &Th, const Triangle &K, const RdHat
    // &PHat,
    //         RNMK_ &val) const;
    void FB(const What_d whatd, const Element &K, const Rd &PHat, RNMK_ &val) const;
};

// ENDOFCLASS TypeOfFE_PkEdge
int TypeOfFE_BDM1_2d::Data[] = {3, 3, 4, 4, 5, 5, // support on what
                                0, 1, 0, 1, 0, 1, // df on node
                                0, 0, 1, 1, 2, 2, // th node of df
                                // 0, 0, 0, 0, 0, 0,    // df previou FE
                                0, 1, 2, 3, 4, 5, // which df on prev
                                0, 1, 0, 0,
                                // 0, 0, 0, 0, 6, 6};
                                0, 0, 6};
void TypeOfFE_BDM1_2d::FB(const What_d whatd, const Element &K, const R2 &PHat, RNMK_ &bfMat) const {

    R2 X   = K(PHat);
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R refBaryc[3] = {l0, l1, l2};
    R2 Dl[3]      = {K.H(0), K.H(1), K.H(2)};

    assert(bfMat.N() >= 6);
    assert(bfMat.M() == 2);

    bfMat      = 0;
    R cK       = 2 * K.measure();
    double s   = sqrt(K.measure());   // gives h^-3 for Taylor-Green vortex problem
    //double s = 1.;                  // gives h^-1.15 for Taylor-Green vortex problem
    //double s   = 1./sqrt(K.measure());  // gives decreasing condition numbers
    int ortho0 = 0, ortho1 = 1;
    R s1ortho = 1;
    if (whatd & Fop_D0) {
        for (int df = 0, e = 0; e < 3; ++e) {
            int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
            R eOrientation = K.EdgeOrientation(e);
            R2 f1          = (X - Q[e]) * eOrientation / cK;
            R2 f2          = -(Dl[e1] * refBaryc[e2] + Dl[e2] * refBaryc[e1]).perp();

            // R2 f1 = 1./Dl[e].norm()*(Q[e1] - Q[e])*refBaryc[e1];
            // R2 f2 = 1./Dl[e].norm()*(Q[e2] - Q[e])*refBaryc[e2];

            bfMat(df, ortho0, op_id)   = f1.x * s;
            bfMat(df++, ortho1, op_id) = s1ortho * f1.y * s;

            bfMat(df, ortho0, op_id)   = f2.x * s;
            bfMat(df++, ortho1, op_id) = s1ortho * f2.y * s;
        }
    }

    if ((whatd & Fop_D1) || (whatd & Fop_D2)) {
        R2 Dphix(1, 0);
        R2 Dphiy(0, 1);

        if (whatd & Fop_dx) {
            for (int df = 0, e = 0; e < 3; ++e) {
                int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
                R eOrientation = K.EdgeOrientation(e);
                R2 f1          = R2(eOrientation / cK, 0.);
                R2 f2          = -(Dl[e1] * Dl[e2].x + Dl[e2] * Dl[e1].x).perp();

                bfMat(df, ortho0, op_dx)   = f1.x * s;
                bfMat(df++, ortho1, op_dx) = s1ortho * f1.y * s;

                bfMat(df, ortho0, op_dx)   = f2.x * s;
                bfMat(df++, ortho1, op_dx) = s1ortho * f2.y * s;
            }
        }

        if (whatd & Fop_dy) {
            for (int df = 0, e = 0; e < 3; ++e) {
                int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
                R eOrientation = K.EdgeOrientation(e);
                R2 f1          = R2(0., eOrientation / cK);
                R2 f2          = -(Dl[e1] * Dl[e2].y + Dl[e2] * Dl[e1].y).perp();

                bfMat(df, ortho0, op_dy)   = f1.x * s;
                bfMat(df++, ortho1, op_dy) = s1ortho * f1.y * s;

                bfMat(df, ortho0, op_dy)   = f2.x * s;
                bfMat(df++, ortho1, op_dy) = s1ortho * f2.y * s;
            }
        }
    }
}

static TypeOfFE_BDM1_2d myBDM1_2d;
GTypeOfFE<Mesh2> &BDM1_2d(myBDM1_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::BDM1 = myBDM1_2d;
