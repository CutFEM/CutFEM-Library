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
#include "../num/DA.hpp"

struct InitTypeOfRTk_2d {
    int k;    // order poly on edge
    int ndfi; // nb of internal dof
    int npe;  // nb point on edge
    int ndf;  // nb dof

    KN<R> X;      // point on edge
    KN<int> Data; // data of TypeOfFE
    const QuadratureFormular1d QFE;
    const GQuadratureFormular<R2> &QFK;

    Vertex2 verticesHat[3];
    Triangle2 TriangleHat;

    InitTypeOfRTk_2d(int KK)
        : k(KK), ndfi((k + 1) * (k)), npe(k + 1), ndf(3 * npe + ndfi), Data(4 * ndf + 7),
          QFE(-1 + 2 * npe, npe, GaussLegendre(npe), true), QFK(QuadratureFormular_T_5) {
        int ndfe = ndf - ndfi; //
        int o[5];

        o[0] = 0;
        for (int i = 1; i < 5; ++i) {
            // o[i] = o[i - 1] + ndf;
            o[i] = i * ndf;
        }

        for (int dof = 0; dof < ndf; ++dof) {
            if (dof < ndfe) {
                int e            = dof / npe;
                int n            = dof % npe;
                Data[o[0] + dof] = 3 + e;
                Data[o[1] + dof] = n;
                Data[o[2] + dof] = e;
                Data[o[3] + dof] = dof;
            } else {
                int n            = dof - ndfe;
                Data[o[0] + dof] = 6;
                Data[o[1] + dof] = n;
                Data[o[2] + dof] = 3;
                Data[o[3] + dof] = dof;
            }
        }

        Data[o[4] + 0] = 0;
        Data[o[4] + 1] = 1;
        Data[o[4] + 2] = 1;
        Data[o[4] + 3] = 0;
        Data[o[4] + 4] = 0;
        Data[o[4] + 5] = 0;
        Data[o[4] + 6] = ndf; // end_dfcomp

        (R2 &)verticesHat[0] = R2::KHat[0];
        (R2 &)verticesHat[1] = R2::KHat[1];
        (R2 &)verticesHat[2] = R2::KHat[2];
        int iv[3]            = {0, 1, 2};
        TriangleHat.set(verticesHat, iv, 0);
    }
};

class TypeOfFE_RT1_2d : public InitTypeOfRTk_2d, public GTypeOfFE<Mesh2> {
    typedef Mesh2 Mesh;               // Define 2D mesh as Mesh
    typedef typename Mesh::Element E; // Define mesh element as E

  public:
    static double Pi_h_coef[];
    // bool Ortho;

    TypeOfFE_RT1_2d()
        : InitTypeOfRTk_2d(1), GTypeOfFE<Mesh2>(ndf, 2, Data,
                                                2 * 2 * 3 * QFE.n + QFK.n * 4, // nb coef mat interpole
                                                3 * QFE.n + QFK.n              // nb P interpolation
                               ) {
        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::RT1;
        GTypeOfFE<Mesh>::polynomialOrder = 2;

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
        int i6 = 6, i7 = 7;
        for (int p = 0; p < QFK.n; ++p) {
            ipj_Pi_h[kkk++] = IPJ(i6, i, 0);
            ipj_Pi_h[kkk++] = IPJ(i6, i, 1);
            ipj_Pi_h[kkk++] = IPJ(i7, i, 0);
            ipj_Pi_h[kkk++] = IPJ(i7, i, 1);
            Pt_Pi_h[i++]    = QFK[p];
        }

        assert(kkk == this->ipj_Pi_h.N());
        assert(i == this->Pt_Pi_h.N());
    }

    void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const { // compute the coef of interpolation ...
        const Element &T     = K.T;
        int k                = 0;
        int arrEdgeOrient[3] = {T.EdgeOrientation(0), T.EdgeOrientation(1), T.EdgeOrientation(2)};
        double aa            = T.measure();
        double s[3]          = {1., 1., 1.}; //{T.lenEdge(0)/aa, T.lenEdge(1)/aa, T.lenEdge(2)/aa};

        double sb = 1.; // 1./sqrt(T.measure());
        for (int i = 0; i < 3; i++) {
            R2 E(-T.Edge(i).perp());
            R eOrientation = arrEdgeOrient[i];

            for (int p = 0; p < QFE.n; ++p) {
                R l0 = QFE[p].X(), l1 = 1 - QFE[p].X();
                R p0      = (2 * l0 - l1) * 2;            // poly othogonaux to \lambda_1
                R p1      = (2 * l1 - l0) * 2;            // poly othogonaux to \lambda_0
                R lambda1 = eOrientation * p0 * QFE[p].a; // [some quadrature function?]
                R lambda0 = eOrientation * p1 * QFE[p].a; //
                if (eOrientation < 0) {
                    std::swap(lambda1, lambda0); // exch lambda0,lambda1
                }
                v[k++] = lambda0 * E.x * s[i];
                v[k++] = lambda0 * E.y * s[i];
                v[k++] = lambda1 * E.x * s[i];
                v[k++] = lambda1 * E.y * s[i];
            }
        }

        R2 B[2] = {T.Edge(1), T.Edge(2)};
        B[0]    = B[0].perp();
        B[1]    = B[1].perp();

        double CK = 0.5; // dof U= [u1,u2] > |K| int_K ( B_i.U )

        for (int p = 0; p < QFK.n; ++p) {
            double w = QFK[p].a * CK;
            v[k++]   = w * B[0].x * sb;
            v[k++]   = w * B[0].y * sb;
            v[k++]   = w * B[1].x * sb;
            v[k++]   = w * B[1].y * sb;
        }

        assert(k == this->ipj_Pi_h.N());
    }

    void FB(const What_d, const Element &K, const Rd &PHat, RNMK_ &bfMat) const;

  private:
    void FB_Freefem(const What_d, const Element &K, const Rd &PHat, RNMK_ &bfMat) const;
    void FB_ID(const Element &K, const Rd &Phat, RNMK_ &bfMat) const;
    void FB_D1(const Element &K, const Rd &Phat, RNMK_ &bfMat) const;
    void FB_D2(const Element &K, const Rd &Phat, RNMK_ &bfMat) const;
};

void TypeOfFE_RT1_2d::FB(const What_d whatd, const Element &K, const Rd &Phat, RNMK_ &bfMat) const {

    assert(bfMat.N() >= ndf);
    assert(bfMat.M() == 2);
    bfMat = 0;

    // if ((whatd & Fop_D2) ) { FB_D2(K, Phat, bfMat);}
    // // else if ((whatd & Fop_D1) ) { FB_D1(K, Phat, bfMat);}
    // // else {FB_ID(K, Phat, bfMat);}
    // else {
    FB_Freefem(whatd, K, Phat, bfMat);
    // }
}
void TypeOfFE_RT1_2d::FB_Freefem(const What_d whatd, const Element &K, const Rd &Phat, RNMK_ &bfMat) const {
    R2 X   = K(Phat);                        // [Phat is the quadrature point X mapped to the reference
                                             // triangle (in loop of addElementMat, X is sent to Phat)]
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])}; // Triangle K node points
    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x,
      l2                = Phat.y; // [reference triangle barycentric coords]
    R refBaryc[3]       = {l0, l1, l2};
    int arrEdgeOrient[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    /*
     *
     * THE 2 DOF k=0,1  are: on edge e   f -> \int_e f \lambda_{e+k} . n_e
     * THE 2 internal dof are : f -> \int_K f e_i  where e_i is the canonical
     * basis of R^2
     *
     *
     * so the basis function are
     *
     * let call \phi_i the basic fonction of RT0 (without orientation) so the
     * normal is exterior. \phi_i (X) = ( X- Q_i ) / (2 |K|) =  \lambda_{i+1}
     * Curl( \lambda_{i+2}) - \lambda_{i+2} Curl( \lambda_{i+1})
     *
     * edge function j=0,1
     * i1= i+j+1, i2= i+2-j  remark : {i,i1,i2} <=> {i,i+1,i+2}
     * fb_i,j = \phi_i ( \lambda_{i1} - 4/3 \lambda_i) + 1/3
     * \phi_{i1}\lambda_{i1}
     *
     *
     *
     * we have 2 bubbles functions
     * fb_6 =      8    \phi_0   \lambda_{0} +   16 \phi_1   \lambda_{1}
     * fb_7 =     - 8   \phi_0   \lambda_{0} + 8 \phi_1   \lambda_{1}
     *
     *
     * such that  for i,j =0,1
     * int_K ( Bb_j  f_{6+i))  =   \delta_{ij}
     *
     *
     * so all basic d function are the sum of 3 function
     *
     * sum_{k=0}^2  c_k  phi_{p_k} lambda_{l_k}
     *
     */

    assert(bfMat.N() >= ndf);
    assert(bfMat.M() == 2);

    bfMat = 0;

    R triMeas2 = 2 * K.mesure();

    R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]}; // phi * area *2
    R2 X_dx(1, 0);
    R2 X_dy(0, 1);
    R2 phi_dx[3] = {X_dx, X_dx, X_dx};
    R2 phi_dy[3] = {X_dy, X_dy, X_dy};

    R Kmap[2][2]; // Kmap: K -> Khat
    Kmap[0][0] = 1 / triMeas2 * (Q[2].y - Q[0].y);
    Kmap[0][1] = -1 / triMeas2 * (Q[2].x - Q[0].x);
    Kmap[1][0] = -1 / triMeas2 * (Q[1].y - Q[0].y);
    Kmap[1][1] = 1 / triMeas2 * (Q[1].x - Q[0].x);

    // R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y; // [reference
    // triangle barycentric coords] R refBaryc[3] = {l0, l1, l2};
    R2 Phat_dx(Kmap[0][0], Kmap[1][0]);
    R2 Phat_dy(Kmap[0][1], Kmap[1][1]);
    R refBaryc_dx[3] = {1 - Phat_dx.x - Phat_dx.y, Phat_dx.x, Phat_dx.y};
    R refBaryc_dy[3] = {1 - Phat_dy.x - Phat_dy.y, Phat_dy.x, Phat_dy.y};
    int pI[8][3]; // store p_k
    int lI[8][3]; // store l_k
    R cI[8][3];   // store c_k
    int dof     = 0;
    // double sb = sqrt(K.mesure());
    // double aa = K.measure();
    double s[4] = {1., 1., 1., 1.}; //{K.lenEdge(0)/aa, K.lenEdge(1)/aa, K.lenEdge(2)/aa, sb};

    for (int e = 0; e < 3; ++e) { // [loops through edges]
        // int i = e;
        int ii[2]      = {(e + 1) % 3, (e + 2) % 3};
        R eOrientation = arrEdgeOrient[e] / triMeas2;
        if (eOrientation < 0) {
            std::swap(ii[0], ii[1]);
        }

        for (int j = 0; j < 2; ++j, dof++) {
            pI[dof][0] = e;
            lI[dof][0] = ii[j];
            cI[dof][0] = eOrientation;

            pI[dof][1] = e;
            lI[dof][1] = e;
            cI[dof][1] = -eOrientation * 4. / 3.;

            pI[dof][2] = ii[j];
            lI[dof][2] = ii[j];
            cI[dof][2] = eOrientation / 3.;
        }
    } // [now dof = 6, next code piece does bubbles dof=6,7]

    // FB  (x-Q_i) l_i l_j  =
    R s8 = 8 / triMeas2, s01 = s8;
    R cbb[] = {s8, 2 * s01, -s01, s8}; // { [ 8, 16], [ -8, 8] }

    // [the 2 bubbles]
    for (int j = 0; j < 2; ++j, dof++) { // [j indexes the bubble funs]
        pI[dof][0] = 0;                  // i
        lI[dof][0] = 0;
        cI[dof][0] = cbb[j];

        pI[dof][1] = 1; // i
        lI[dof][1] = 1;
        cI[dof][1] = cbb[j + 2];

        pI[dof][2] = 2;
        lI[dof][2] = 2;
        cI[dof][2] = 0; // [the third constant is zero for the bubbles]
    }

    assert(dof == 8);
    if (whatd & Fop_id) {
        for (int dof = 0; dof < 8; ++dof) {
            R2 fd(0., 0.);

            for (int k = 0; k < 3; ++k) {
                fd += (cI[dof][k] * refBaryc[lI[dof][k]]) * phi[pI[dof][k]];
            }

            bfMat(dof, 0, op_id) = fd.x * s[dof / 2];
            bfMat(dof, 1, op_id) = fd.y * s[dof / 2];
        }
    }

    if ((whatd & Fop_D1) || (whatd & Fop_D2)) {
        R2 DL[3] = {K.H(0), K.H(1), K.H(2)}; // [heights?? represents derivative of refBaryc]
        R2 e1(1, 0);                         // [initialise R2 canonical basis vecs]
        R2 e2(0, 1);

        if (whatd & Fop_dx) {
            for (int dof = 0; dof < 8; ++dof) {
                R2 fd(0., 0.); // init 0 vec

                for (int k = 0; k < 3; ++k) {
                    // std::cout << pI[dof][k] << std::endl;
                    fd += cI[dof][k] * (DL[lI[dof][k]].x * phi[pI[dof][k]] + refBaryc[lI[dof][k]] * e1);
                }

                bfMat(dof, 0, op_dx) = fd.x * s[dof / 2];
                bfMat(dof, 1, op_dx) = fd.y * s[dof / 2];
            }
        }

        if (whatd & Fop_dy) {
            for (int dof = 0; dof < 8; ++dof) {
                R2 fd(0., 0.);

                for (int k = 0; k < 3; ++k) {
                    fd += cI[dof][k] * (DL[lI[dof][k]].y * phi[pI[dof][k]] + refBaryc[lI[dof][k]] * e2);
                }

                bfMat(dof, 0, op_dy) = fd.x * s[dof / 2];
                bfMat(dof, 1, op_dy) = fd.y * s[dof / 2];
            }
        }

        if (whatd & Fop_D2) {
            if (whatd & Fop_dxx) {
                for (int dof = 0; dof < 8; ++dof) {
                    R2 fd(0., 0.); // init 0 vec

                    for (int k = 0; k < 3; ++k) { // [take expression from Fop_dx and diff wrt x]
                        // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dx[pI[dof][k]] +
                        // refBaryc_dx[lI[dof][k]] * e1); fd += cI[dof][k] *
                        // (DL[lI[dof][k]].x * phi_dx[pI[dof][k]] + DL[lI[dof][k]].x *
                        // e1);
                        fd += cI[dof][k] * (DL[lI[dof][k]].x * e1 + DL[lI[dof][k]].x * e1);
                    }

                    bfMat(dof, 0, op_dxx) = fd.x * s[dof / 2];
                    bfMat(dof, 1, op_dxx) = fd.y * s[dof / 2];
                }
            }

            if (whatd & Fop_dxy) {
                for (int dof = 0; dof < 8; ++dof) {
                    R2 fd(0., 0.); // init 0 vec

                    for (int k = 0; k < 3; ++k) { // [take expression from Fop_dx and diff wrt y]
                        // std::cout << pI[dof][k] << std::endl;
                        // fd += cI[dof][k] * (DL[lI[dof][k]].x * phi_dy[pI[dof][k]] +
                        // refBaryc_dy[lI[dof][k]] * e1); fd += cI[dof][k] *
                        // (DL[lI[dof][k]].x * phi_dy[pI[dof][k]] + DL[lI[dof][k]].y *
                        // e1);
                        fd += cI[dof][k] * (DL[lI[dof][k]].x * e2 + DL[lI[dof][k]].y * e1);
                    }

                    bfMat(dof, 0, op_dxy) = fd.x * s[dof / 2];
                    bfMat(dof, 1, op_dxy) = fd.y * s[dof / 2];
                }
            }

            if (whatd & Fop_dyy) {
                for (int dof = 0; dof < 8; ++dof) {
                    R2 fd(0., 0.);

                    for (int k = 0; k < 3; ++k) {
                        // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dy[pI[dof][k]] +
                        // refBaryc_dy[lI[dof][k]] * e2); fd += cI[dof][k] *
                        // (DL[lI[dof][k]].y * phi_dy[pI[dof][k]] + DL[lI[dof][k]].y *
                        // e2);
                        fd += cI[dof][k] * (DL[lI[dof][k]].y * e2 + DL[lI[dof][k]].y * e2);
                    }

                    bfMat(dof, 0, op_dyy) = fd.x * s[dof / 2];
                    bfMat(dof, 1, op_dyy) = fd.y * s[dof / 2];
                }
            }

            // if (whatd & Fop_dyx) { // [MIXED PARTIALS DO NOT COMMUTE]
            //   for (int dof = 0; dof < 8; ++dof) {
            //     R2 fd(0., 0.);
            //     // R2 test(0., 0.);
            //
            //     for (int k = 0; k < 3; ++k) {// [take expression from Fop_dy and
            //     diff wrt x]
            //       // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dx[pI[dof][k]] +
            //       refBaryc_dx[lI[dof][k]] * e2);
            //       // fd += cI[dof][k] * (DL[lI[dof][k]].y * phi_dx[pI[dof][k]] +
            //       DL[lI[dof][k]].x * e2); fd += cI[dof][k] * (DL[lI[dof][k]].y *
            //       e1 + DL[lI[dof][k]].x * e2);
            //     }
            //     // std::cout << (fd-test) << std::endl;
            //     bfMat(dof, 0, op_dyx) = fd.x;
            //     bfMat(dof, 1, op_dyx) = fd.y;
            //   }
            // }
        }
    }
}
// void TypeOfFE_RT1_2d::FB_ID(const Element &K, const Rd &Phat,
//                             RNMK_ &bfMat) const {

//    R2 X   = K(Phat); // [Phat is the quadrature point X mapped to the
//    reference
//                      // triangle (in loop of addElementMat, X is sent to
//                      Phat)]
//    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])}; // Triangle K node points
//    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x,
//      l2                = Phat.y; // [reference triangle barycentric coords]
//    R refBaryc[3]       = {l0, l1, l2};
//    int arrEdgeOrient[] = {K.EdgeOrientation(0), K.EdgeOrientation(1),
//                           K.EdgeOrientation(2)};

//    R triMeas2 = 2 * K.mesure();
//    R2 phi[3]  = {X - Q[0], X - Q[1], X - Q[2]}; // phi * area *2

//    int pI[8][3]; // store p_k
//    int lI[8][3]; // store l_k
//    R cI[8][3];   // store c_k
//    int dof  = 0;
//    double s = sqrt(K.mesure());

//    for (int e = 0; e < 3; ++e) {
//       // int i = e;
//       int ii[2]      = {(e + 1) % 3, (e + 2) % 3};
//       R eOrientation = arrEdgeOrient[e] / triMeas2;
//       if (eOrientation < 0) {
//          std::swap(ii[0], ii[1]);
//       }

//       for (int j = 0; j < 2; ++j, dof++) {
//          pI[dof][0] = e;
//          lI[dof][0] = ii[j];
//          cI[dof][0] = eOrientation;

//          pI[dof][1] = e;
//          lI[dof][1] = e;
//          cI[dof][1] = -eOrientation * 4. / 3.;

//          pI[dof][2] = ii[j];
//          lI[dof][2] = ii[j];
//          cI[dof][2] = eOrientation / 3.;
//       }
//    }

//    R s8 = 8 / triMeas2, s01 = s8;
//    R cbb[] = {s8, 2 * s01, -s01, s8}; // { [ 8, 16], [ -8, 8] }
//    // [the 2 bubbles]
//    for (int j = 0; j < 2; ++j, dof++) { // [j indexes the bubble funs]
//       pI[dof][0] = 0;                   // i
//       lI[dof][0] = 0;
//       cI[dof][0] = cbb[j];

//       pI[dof][1] = 1; // i
//       lI[dof][1] = 1;
//       cI[dof][1] = cbb[j + 2];

//       pI[dof][2] = 2;
//       lI[dof][2] = 2;
//       cI[dof][2] = 0; // [the third constant is zero for the bubbles]
//    }

//    assert(dof == 8);
//    for (int dof = 0; dof < 8; ++dof) {
//       R2 fd(0., 0.);

//       for (int k = 0; k < 3; ++k) {
//          fd += (cI[dof][k] * refBaryc[lI[dof][k]]) * phi[pI[dof][k]];
//       }

//       bfMat(dof, 0, op_id) = fd.x * s;
//       bfMat(dof, 1, op_id) = fd.y * s;
//    }
// };
// void TypeOfFE_RT1_2d::FB_D1(const Element &K, const Rd &Phat,
//                             RNMK_ &bfMat) const {
//    R2 X   = K(Phat); // [Phat is the quadrature point X mapped to the
//    reference
//                      // triangle (in loop of addElementMat, X is sent to
//                      Phat)]
//    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])}; // Triangle K node points
//    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x,
//      l2                = Phat.y; // [reference triangle barycentric coords]
//    R refBaryc[3]       = {l0, l1, l2};
//    int arrEdgeOrient[] = {K.EdgeOrientation(0), K.EdgeOrientation(1),
//                           K.EdgeOrientation(2)};
//    R triMeas2          = 2 * K.mesure();
//    R2 phi[3]           = {X - Q[0], X - Q[1], X - Q[2]}; // phi * area *2

//    R2 ddd[3] = {K.H(0), K.H(1), K.H(2)};
//    Diff<double, 2> ll1(l1, 0), ll2(l2, 1);
//    ll1.d[0]            = ddd[1].x;
//    ll1.d[1]            = ddd[1].y;
//    ll2.d[0]            = ddd[2].x;
//    ll2.d[1]            = ddd[2].y;
//    Diff<double, 2> ll0 = 1 - ll1 - ll2;
//    Diff<double, 2> Xx(X.x, 0), Xy(X.y, 1);
//    Diff<double, 2> LL[3]   = {ll0, ll1, ll2};
//    Diff<double, 2> PHIx[3] = {Xx - Q[0].x, Xx - Q[1].x, Xx - Q[2].x};
//    Diff<double, 2> PHIy[3] = {Xy - Q[0].y, Xy - Q[1].y, Xy - Q[2].y};
//    int pI[8][3]; // store p_k
//    int lI[8][3]; // store l_k
//    R cI[8][3];   // store c_k
//    int dof  = 0;
//    double s = sqrt(K.mesure());

//    for (int e = 0; e < 3; ++e) {
//       // int i = e;
//       int ii[2]      = {(e + 1) % 3, (e + 2) % 3};
//       R eOrientation = arrEdgeOrient[e] / triMeas2;
//       if (eOrientation < 0) {
//          std::swap(ii[0], ii[1]);
//       }

//       for (int j = 0; j < 2; ++j, dof++) {
//          pI[dof][0] = e;
//          lI[dof][0] = ii[j];
//          cI[dof][0] = eOrientation;

//          pI[dof][1] = e;
//          lI[dof][1] = e;
//          cI[dof][1] = -eOrientation * 4. / 3.;

//          pI[dof][2] = ii[j];
//          lI[dof][2] = ii[j];
//          cI[dof][2] = eOrientation / 3.;
//       }
//    }

//    R s8 = 8 / triMeas2, s01 = s8;
//    R cbb[] = {s8, 2 * s01, -s01, s8}; // { [ 8, 16], [ -8, 8] }

//    for (int j = 0; j < 2; ++j, dof++) {
//       pI[dof][0] = 0;
//       lI[dof][0] = 0;
//       cI[dof][0] = cbb[j];

//       pI[dof][1] = 1; // i
//       lI[dof][1] = 1;
//       cI[dof][1] = cbb[j + 2];

//       pI[dof][2] = 2;
//       lI[dof][2] = 2;
//       cI[dof][2] = 0;
//    }

//    assert(dof == 8);

//    for (int dof = 0; dof < 8; ++dof) {
//       R2 fd(0., 0.);
//       Diff<double, 2> FDx, FDy;

//       for (int k = 0; k < 3; ++k) {
//          fd += (cI[dof][k] * refBaryc[lI[dof][k]]) * phi[pI[dof][k]];
//          FDx += (cI[dof][k] * LL[lI[dof][k]]) * PHIx[pI[dof][k]];
//          FDy += (cI[dof][k] * LL[lI[dof][k]]) * PHIy[pI[dof][k]];
//       }

//       bfMat(dof, 0, op_id) = FDx.val * s;
//       bfMat(dof, 1, op_id) = FDy.val * s;
//       bfMat(dof, 0, op_dx) = FDx.d[0] * s;
//       bfMat(dof, 1, op_dx) = FDy.d[0] * s;
//       bfMat(dof, 0, op_dy) = FDx.d[1] * s;
//       bfMat(dof, 1, op_dy) = FDy.d[1] * s;
//    }
// };
// void TypeOfFE_RT1_2d::FB_D2(const Element &K, const Rd &Phat,
//                             RNMK_ &bfMat) const {
//    R2 X   = K(Phat); // [Phat is the quadrature point X mapped to the
//    reference
//                      // triangle (in loop of addElementMat, X is sent to
//                      Phat)]
//    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])}; // Triangle K node points
//    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x,
//      l2                = Phat.y; // [reference triangle barycentric coords]
//    R refBaryc[3]       = {l0, l1, l2};
//    int arrEdgeOrient[] = {K.EdgeOrientation(0), K.EdgeOrientation(1),
//                           K.EdgeOrientation(2)};
//    R triMeas2          = 2 * K.mesure();
//    R2 phi[3]           = {X - Q[0], X - Q[1], X - Q[2]}; // phi * area *2

//    R2 ddd[3] = {K.H(0), K.H(1), K.H(2)};
//    Diff<Diff<R, 2>, 2> ll1(l1, 0), ll2(l2, 1);
//    ll1.d[0]     = ddd[1].x;
//    ll1.d[1]     = ddd[1].y;
//    ll2.d[0]     = ddd[2].x;
//    ll2.d[1]     = ddd[2].y;
//    ll1.val.d[0] = ddd[1].x;
//    ll1.val.d[1] = ddd[1].y; // init val of dx
//    ll2.val.d[0] = ddd[2].x;
//    ll2.val.d[1] = ddd[2].y;

//    Diff<Diff<R, 2>, 2> ll0 = 1 - ll1 - ll2;
//    Diff<Diff<R, 2>, 2> Xx(X.x, 0), Xy(X.y, 1);
//    Xx.val.d[0]                 = 1; // init val of dx
//    Xy.val.d[1]                 = 1;
//    Diff<Diff<R, 2>, 2> LL[3]   = {ll0, ll1, ll2};
//    Diff<Diff<R, 2>, 2> PHIx[3] = {Xx - Q[0].x, Xx - Q[1].x, Xx - Q[2].x};
//    Diff<Diff<R, 2>, 2> PHIy[3] = {Xy - Q[0].y, Xy - Q[1].y, Xy - Q[2].y};

//    int pI[8][3]; // store p_k
//    int lI[8][3]; // store l_k
//    R cI[8][3];   // store c_k
//    int dof  = 0;
//    double s = sqrt(K.mesure());

//    for (int e = 0; e < 3; ++e) {
//       // int i = e;
//       int ii[2]      = {(e + 1) % 3, (e + 2) % 3};
//       R eOrientation = arrEdgeOrient[e] / triMeas2;
//       if (eOrientation < 0) {
//          std::swap(ii[0], ii[1]);
//       }

//       for (int j = 0; j < 2; ++j, dof++) {
//          pI[dof][0] = e;
//          lI[dof][0] = ii[j];
//          cI[dof][0] = eOrientation;

//          pI[dof][1] = e;
//          lI[dof][1] = e;
//          cI[dof][1] = -eOrientation * 4. / 3.;

//          pI[dof][2] = ii[j];
//          lI[dof][2] = ii[j];
//          cI[dof][2] = eOrientation / 3.;
//       }
//    }

//    R s8 = 8 / triMeas2, s01 = s8;
//    R cbb[] = {s8, 2 * s01, -s01, s8}; // { [ 8, 16], [ -8, 8] }

//    for (int j = 0; j < 2; ++j, dof++) {
//       pI[dof][0] = 0;
//       lI[dof][0] = 0;
//       cI[dof][0] = cbb[j];

//       pI[dof][1] = 1; // i
//       lI[dof][1] = 1;
//       cI[dof][1] = cbb[j + 2];

//       pI[dof][2] = 2;
//       lI[dof][2] = 2;
//       cI[dof][2] = 0;
//    }

//    assert(dof == 8);

//    for (int dof = 0; dof < 8; ++dof) {
//       R2 fd(0., 0.);
//       Diff<Diff<double, 2>, 2> FDx, FDy;

//       for (int k = 0; k < 3; ++k) {
//          fd += (cI[dof][k] * refBaryc[lI[dof][k]]) * phi[pI[dof][k]];
//          FDx += (cI[dof][k] * LL[lI[dof][k]]) * PHIx[pI[dof][k]];
//          FDy += (cI[dof][k] * LL[lI[dof][k]]) * PHIy[pI[dof][k]];
//       }

//       bfMat(dof, 0, op_id)  = FDx.val.val * s;
//       bfMat(dof, 1, op_id)  = FDy.val.val * s;
//       bfMat(dof, 0, op_dx)  = FDx.d[0].val * s;
//       bfMat(dof, 1, op_dx)  = FDy.d[0].val * s;
//       bfMat(dof, 0, op_dy)  = FDx.d[1].val * s;
//       bfMat(dof, 1, op_dy)  = FDy.d[1].val * s;
//       bfMat(dof, 0, op_dxx) = FDx.d[0].d[0] * s;
//       bfMat(dof, 1, op_dxx) = FDy.d[0].d[0] * s;
//       bfMat(dof, 0, op_dxy) = FDx.d[0].d[1] * s;
//       bfMat(dof, 1, op_dxy) = FDy.d[0].d[1] * s;
//       bfMat(dof, 0, op_dyy) = FDx.d[1].d[1] * s;
//       bfMat(dof, 1, op_dyy) = FDy.d[1].d[1] * s;
//    }
// };

static TypeOfFE_RT1_2d myRT1_2d;
GTypeOfFE<Mesh2> &RT1_2d(myRT1_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::RT1 = myRT1_2d;

// class TypeOfFE_RT1_3d : public GTypeOfFE< Mesh3 > {
//  public:
//   typedef Mesh3 Mesh;
//   typedef Mesh3::Element Element;
//   typedef GFElement< Mesh3 > FElement;
//   static int dfon[];
//   static const int d = Mesh::Rd::d;
//   // quadrature Formula on a face
//   static const GQuadratureFormular< R2 > QFface;
//   // quadrature Formula on an element
//   static const GQuadratureFormular< R3 > QFtetra;
//   TypeOfFE_RT1_3d( );
//   int edgeface[4][3];
//   void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const
//   RdHat &PHat,
//           RNMK_ &val) const;
//   void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
//   int ocoef, int odf,
//            int *nump) const;
// };
//
// int TypeOfFE_RT1_3d::dfon[] = {0, 0, 3, 3};    // dofs per vertice, edge,
// face, volume
//
// // Quadrature formula on a face,
// const GQuadratureFormular< R2 >
// TypeOfFE_RT1_3d::QFface(QuadratureFormular_T_5);
//
// // Quadrature formula on the tetraedron
// const GQuadratureFormular< R3 >
// TypeOfFE_RT1_3d::QFtetra(QuadratureFormular_Tet_2);
//
// TypeOfFE_RT1_3d::TypeOfFE_RT1_3d( )
//   : GTypeOfFE< Mesh3 >(dfon, d, 1, 3 * QFface.n * 3 * Element::nf + 3 *
//   QFtetra.n * 3,
//                        Element::nf * QFface.n + QFtetra.n, false, true) {
//   assert(QFface.n);
//   assert(QFtetra.n);
//   // 4 ref tetra vertices
//   R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};
//
//   // We build the interpolation pts on the faces
//   int p = 0;
//
//   for (int f = 0; f < Element::nf; ++f) {
//     for (int q = 0; q < QFface.n; ++q, ++p) {
//       double x = QFface[q].x;
//       double y = QFface[q].y;
//       this->PtInterpolation[p] = Pt[Element::nvface[f][0]] * (1. - x - y) +
//                                  Pt[Element::nvface[f][1]] * x +
//                                  Pt[Element::nvface[f][2]] * y;
//     }
//   }
//
//   // We build the interpolation bubble pts
//   for (int q = 0; q < QFtetra.n; ++q, ++p) {
//     double x = QFtetra[q].x;
//     double y = QFtetra[q].y;
//     double z = QFtetra[q].z;
//     this->PtInterpolation[p] = Pt[0] * (1. - x - y - z) + Pt[1] * x + Pt[2] *
//     y + Pt[3] * z;
//   }
//
//   int i = 0;
//   p = 0;
//
//   for (int f = 0; f < Element::nf; f++) {        // loop on the 4 face dofs
//     for (int q = 0; q < QFface.n; ++q, p++) {    // loop on the face
//     quadrature pts
//       for (int df = 0; df < 3; df++) {           // 3 dof par face
//         int dof = 3 * f + df;                    // numero  du dof  3 dof /
//         face
//
//         for (int c = 0; c < 3; c++, i++) {    // loop on the 3 components
//           this->pInterpolation[i] = p;        // pk
//           this->cInterpolation[i] = c;        // jk
//           this->dofInterpolation[i] = dof;    // ik
//           this->coefInterpolation[i] = 0.;
//         }
//       }
//     }
//   }
//
//   {
//     p = Element::nf * QFface.n;
//
//     for (int q = 0; q < QFtetra.n; ++q, ++p) {    // loop on the volume
//     quadrature pts
//       for (int v = 12; v < 15; v++) {             // loop on the 3 volume
//       dofs
//         for (int c = 0; c < 3; c++, i++) {        // loop on the 3 components
//           this->pInterpolation[i] = p;            // pk
//           this->cInterpolation[i] = c;            // jk
//           this->dofInterpolation[i] = v;          // ik
//           this->coefInterpolation[i] = 0.;
//         }
//       }
//     }
//   }
//   // verif bonne taille
//   ffassert(p == this->PtInterpolation.N( ));
// }    // end TypeOfFE_RT1_3d()
//
// // For the coefficients of interpolation alphak in (13.1)
//
// void TypeOfFE_RT1_3d::set(const Mesh &Th, const Element &K,
// InterpolationMatrix< RdHat > &M,
//                           int ocoef, int odf, int *nump) const {
//   int i = ocoef;
//
//   // *******************************************
//   // DOFs on the 4 faces --- 3 DOFs per face //
//   // *******************************************
//
//   // inv of int(F) lamda_i lambda_j = Area(K) * (d! (1+delta ij) / (d+max n)!
//   = 0.5 * 2(1+delta
//   // ij) /(2+2)!
//   double c1[][3] = {{9, -3, -3} /* 0 */, {-3, 9, -3} /* 1 */, {-3, -3, 9} /*
//   2 */}; R3 NK[4];
//
//   K.Gradlambda(NK);    // InteriorNormal of F /h = - N  3*|K|/|F|
//   double coefK = -3. * K.mesure( );
//
//   for (int ff = 0; ff < Element::nf; ff++) {    // loop on the 4 face dofs
//     const Element::Vertex *fV[3] = {&K.at(Element::nvface[ff][0]),
//     &K.at(Element::nvface[ff][1]),
//                                     &K.at(Element::nvface[ff][2])};
//     int p[] = {0, 1, 2};
//     int fp = K.facePermutation(ff);
//     if (fp & 1) {
//       std::swap(p[0], p[1]);
//     }
//
//     if (fp & 2) {
//       std::swap(p[1], p[2]);
//     }
//
//     if (fp & 4) {
//       std::swap(p[0], p[1]);
//     }
//
//     R3 N = NK[ff];
//     N *= coefK * K.faceOrient(ff);
//
//     for (int q = 0; q < QFface.n; ++q) {    // loop on the face quadrature
//     pts
//       // the 3 lambdas of P1 on the face
//
//       double lambda[3] = {1. - QFface[q].x - QFface[q].y, QFface[q].x,
//       QFface[q].y}; R sa = QFface[q].a; R cp[3] = {(c1[0][0] * lambda[0] +
//       c1[0][1] * lambda[1] + c1[0][2] * lambda[2]) * sa,
//                  (c1[1][0] * lambda[0] + c1[1][1] * lambda[1] + c1[1][2] *
//                  lambda[2]) * sa, (c1[2][0] * lambda[0] + c1[2][1] *
//                  lambda[1] + c1[2][2] * lambda[2]) * sa};
//
//       for (int idof = 0; idof < 3; idof++) {    // loop sur 3 dof de la face
//         for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
//           M.coef[i] = N[c] * cp[p[idof]];
//         }
//       }
//     }
//   }
//
//   // cout << endl;
//   // **************************************************
//   // DOFs on the tetraedra --- 3 DOFs in the volume //
//   // **************************************************
//   // Base Piola compatible  B_i =N_i* |F_i|/6  for 3 face 1,2,3
//   double CK = -K.mesure( );    // dof U= [u1,u2,u3] > |K| int_K ( B_i.U )
//
//   for (int p = 0; p < QFtetra.n; ++p) {
//     double w = QFtetra[p].a * CK;
//
//     for (int l = 0; l < 3; l++) {
//       M.coef[i++] = w * NK[l + 1].x;
//       M.coef[i++] = w * NK[l + 1].y;
//       M.coef[i++] = w * NK[l + 1].z;
//     }
//   }
//
// }    // end set function
//
// // here the basis functions and theirs derivates
// void TypeOfFE_RT1_3d::FB(const What_d whatd, const Mesh &Th, const
// Mesh3::Element &K,
//                          const RdHat &PHat, RNMK_ &val) const {
//   assert(val.N( ) >= 15);
//   assert(val.M( ) == 3);
//
//   val = 0;
//
//   // basis functions to RT03d / multiply by sign to have a exterior normal
//   and divide by the
//   // mesure of K phi = signe * (x - qi)/ (volume*d)
//   R cc = d * K.mesure( );
//   R lambda[] = {1. - PHat.sum( ), PHat.x, PHat.y, PHat.z};
//   R3 X = K(PHat);
//   R3 phi[4] = {X - K[0], X - K[1], X - K[2], X - K[3]};    // phi * area *6
//
//   // fo contain just the sign about permutation ----- 1perm=-1 / 2perm=1 /
//   3perm=-1 double fo[4] = {(double)K.faceOrient(0), (double)K.faceOrient(1),
//   (double)K.faceOrient(2),
//                   (double)K.faceOrient(3)};
//   int p[15] = {2, 1, 0,  3,  4,  5,  8, 7,
//                6, 9, 10, 11, 12, 13, 14};    // Permutation for orientation
//                to dof
//   R3 Pm[16];                                 // all the momome function ..
//
//   for (int ff = 0, k = 0; ff < Element::nf; ff++, k += 3) {
//     // orientation de la face a envert
//     int fp = K.facePermutation(ff);
//     if (fp & 1) {
//       std::swap(p[k], p[k + 1]);
//     }
//
//     if (fp & 2) {
//       std::swap(p[k + 1], p[k + 2]);
//     }
//
//     if (fp & 4) {
//       std::swap(p[k], p[k + 1]);
//     }
//   }
//
//   double cf[][3] = {
//     {-1.25, 0.25, 0} /* 0 */,  {-1.25, 0, 0.25} /* 1 */,   {-1.5, -0.25,
//     -0.25} /* 2 */, {0.25, -1.25, 0} /* 3 */,  {0, -1.25, 0.25} /* 4 */,
//     {-0.25, -1.5, -0.25} /* 5 */, {0.25, 0, -1.25} /* 6 */,  {0, 0.25, -1.25}
//     /* 7 */,   {-0.25, -0.25, -1.5} /* 8 */, {1.5, 1.25, 1.25} /* 9 */,
//     {1.25, 1.5, 1.25} /* 10 */, {1.25, 1.25, 1.5} /* 11 */,
//     {-15, 15, 0} /* 12 */,     {-15, 0, 15} /* 13 */,      {-30, -15, -15} /*
//     14 */};
//   int Bii[][2] = {{0, 0} /* 0 */,  {0, 1} /* 1 */,  {0, 2} /* 2 */,  {0, 3}
//   /* 3 */,
//                   {1, 0} /* 4 */,  {1, 1} /* 5 */,  {1, 2} /* 6 */,  {1, 3}
//                   /* 7 */, {2, 0} /* 8 */,  {2, 1} /* 9 */,  {2, 2} /* 10 */,
//                   {2, 3} /* 11 */, {3, 0} /* 12 */, {3, 1} /* 13 */, {3, 2}
//                   /* 14 */, {3, 3} /* 15 */};
//   int fe[] = {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14};
//   int k6[] = {0, 5, 10};
//
//   // here we build all the monomials phi_i lambda_j
//   for (int l = 0; l < 16; ++l) {
//     int i = Bii[l][0];
//     int j = Bii[l][1];
//     Pm[l] = phi[i] * lambda[j] / cc;
//   }
//
//   // caution to the numbering of the nodes of the face to the numbering of
//   the basis function
//   // cci
//
//   double sg[15] = {fo[0], fo[0], fo[0], fo[1], fo[1], fo[1], fo[2], fo[2],
//                    fo[2], fo[3], fo[3], fo[3], 1.,    1.,    1.};
//
//   if (whatd & Fop_D0) {
//     for (int pdf = 0; pdf < 15; ++pdf) {
//       int df = p[pdf];
//       R3 fd(0., 0., 0.);
//       if (df < 12) {
//         fd = Pm[fe[df]];    // edge function ..
//       }
//
//       for (int k = 0; k < 3; ++k) {
//         fd += cf[df][k] * Pm[k6[k]];
//       }
//
//       fd *= sg[pdf];
//       val(pdf, 0, op_id) = fd.x;
//       val(pdf, 1, op_id) = fd.y;
//       val(pdf, 2, op_id) = fd.z;
//     }
//   }
//
//   if (whatd & Fop_D1) {
//     R3 DL[4];
//     K.Gradlambda(DL);
//     R3 Dphix(1, 0, 0);
//     R3 Dphiy(0, 1, 0);
//     R3 Dphiz(0, 0, 1);
//     R3 DxPm[16];
//     R3 DyPm[16];
//     R3 DzPm[16];
//
//     for (int l = 0; l < 16; ++l) {
//       // diff phi[i]*lambda[j]/cc;
//       int i = Bii[l][0];
//       int j = Bii[l][1];
//       R Lj = lambda[j];
//       R3 DLj = DL[j];
//       R3 DF1 = (Dphix * Lj + phi[i].x * DLj) / cc;
//       R3 DF2 = (Dphiy * Lj + phi[i].y * DLj) / cc;
//       R3 DF3 = (Dphiz * Lj + phi[i].z * DLj) / cc;
//
//       DxPm[l] = R3(DF1.x, DF2.x, DF3.x);
//       DyPm[l] = R3(DF1.y, DF2.y, DF3.y);
//       DzPm[l] = R3(DF1.z, DF2.z, DF3.z);
//     }
//
//     if (whatd & Fop_dx) {
//       for (int pdf = 0; pdf < 15; ++pdf) {
//         int df = p[pdf];
//         R3 fd(0., 0., 0.);
//         if (df < 12) {
//           fd = DxPm[fe[df]];    // edge function ..
//         }
//
//         for (int k = 0; k < 3; ++k) {
//           fd += cf[df][k] * DxPm[k6[k]];
//         }
//
//         fd *= sg[df];
//         val(pdf, 0, op_dx) = fd.x;
//         val(pdf, 1, op_dx) = fd.y;
//         val(pdf, 2, op_dx) = fd.z;
//       }
//     }
//
//     if (whatd & Fop_dy) {
//       for (int pdf = 0; pdf < 15; ++pdf) {
//         int df = p[pdf];
//         R3 fd(0., 0., 0.);
//         if (df < 12) {
//           fd = DyPm[fe[df]];    // edge function ..
//         }
//
//         for (int k = 0; k < 3; ++k) {
//           fd += cf[df][k] * DyPm[k6[k]];
//         }
//
//         fd *= sg[df];
//         val(pdf, 0, op_dy) = fd.x;
//         val(pdf, 1, op_dy) = fd.y;
//         val(pdf, 2, op_dy) = fd.z;
//       }
//     }
//
//     if (whatd & Fop_dz) {
//       for (int pdf = 0; pdf < 15; ++pdf) {
//         int df = p[pdf];
//         R3 fd(0., 0., 0.);
//         if (df < 12) {
//           fd = DzPm[fe[df]];    // edge function ..
//         }
//
//         for (int k = 0; k < 3; ++k) {
//           fd += cf[df][k] * DzPm[k6[k]];
//         }
//
//         fd *= sg[df];
//         val(pdf, 0, op_dz) = fd.x;
//         val(pdf, 1, op_dz) = fd.y;
//         val(pdf, 2, op_dz) = fd.z;
//       }
//     }
//
//     if (whatd & Fop_D2) {
//       cout << " to do FH RT2 dxx, dyy, dzz, dxy, dxz, dyz " << endl;
//     }
//   }
// }    // end basis functions declaration
