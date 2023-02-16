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

#include "FESpace.hpp"
#include "../num/print_container.hpp"

struct InitTypeOfBDM2_2d {
    int k;    // order poly on edge
    int ndfi; // nb of internal dof
    int npe;  // nb point on edge
    int ndf;  // nb dof

    const QuadratureFormular1d QFE;
    const GQuadratureFormular<R2> &QFK;

    mutable std::vector<std::vector<double>> bf_ref;

    InitTypeOfBDM2_2d()
        : k(2), ndfi(3), npe(3), ndf(3 * npe + ndfi),
          QFE(-1 + 2 * npe, npe, GaussLegendre(npe), true),
          QFK(QuadratureFormular_T_5) {
        bf_ref.resize(ndf, std::vector<double>(2));
    }
};

class TypeOfFE_BDM2_2d : public InitTypeOfBDM2_2d, public GTypeOfFE<Mesh2> {

    typedef Mesh2 Mesh;

  public:
    static int Data[];
    TypeOfFE_BDM2_2d()
        : InitTypeOfBDM2_2d(), GTypeOfFE<Mesh2>(12, 2, Data,
                                                2 * 3 * 3 * QFE.n +
                                                    QFK.n * 3 * 4,
                                                3 * QFE.n + QFK.n) {

        GTypeOfFE<Mesh>::basisFctType    = BasisFctType::BDM2;
        GTypeOfFE<Mesh>::polynomialOrder = 2;

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

                ipj_Pi_h[kkk++] = IPJ(3 * e, i, 0);
                ipj_Pi_h[kkk++] = IPJ(3 * e, i, 1);
                ipj_Pi_h[kkk++] = IPJ(3 * e + 1, i, 0);
                ipj_Pi_h[kkk++] = IPJ(3 * e + 1, i, 1);
                ipj_Pi_h[kkk++] = IPJ(3 * e + 2, i, 0);
                ipj_Pi_h[kkk++] = IPJ(3 * e + 2, i, 1);

                Pt_Pi_h[i++] = B * (QFE[p].X()) + A * (1. - QFE[p].X());
            }
        }
        for (int p = 0; p < QFK.n; ++p) {
            ipj_Pi_h[kkk++] = IPJ(9, i, 0);
            ipj_Pi_h[kkk++] = IPJ(9, i, 1);
            ipj_Pi_h[kkk++] = IPJ(9, i, 0);
            ipj_Pi_h[kkk++] = IPJ(9, i, 1);
            ipj_Pi_h[kkk++] = IPJ(10, i, 0);
            ipj_Pi_h[kkk++] = IPJ(10, i, 1);
            ipj_Pi_h[kkk++] = IPJ(10, i, 0);
            ipj_Pi_h[kkk++] = IPJ(10, i, 1);
            ipj_Pi_h[kkk++] = IPJ(11, i, 0);
            ipj_Pi_h[kkk++] = IPJ(11, i, 1);
            ipj_Pi_h[kkk++] = IPJ(11, i, 0);
            ipj_Pi_h[kkk++] = IPJ(11, i, 1);
            Pt_Pi_h[i++]    = QFK[p];
        }
        assert(kkk == this->ipj_Pi_h.N());
        assert(i == this->Pt_Pi_h.N());
    }

    void get_Coef_Pi_h(const GbaseFElement<Mesh2> &K, KN_<double> &v) const {
        const Element &T = K.T;
        int k            = 0;

        for (int i = 0; i < 3; i++) {
            R2 E(-T.Edge(i).perp());
            R eOrientation = T.EdgeOrientation(i);

            for (int p = 0; p < QFE.n; ++p) {
                R l0 = 1 - QFE[p].X();
                R l1 = QFE[p].X();

                R lambda0 = -eOrientation * l0 * (2 * l0 - 1) * QFE[p].a;
                R lambda1 = -eOrientation * l1 * (2 * l1 - 1) * QFE[p].a;
                R lambda2 = -eOrientation * 4 * l0 * l1 * QFE[p].a;

                if (eOrientation < 0) {
                    std::swap(lambda1, lambda0);
                }
                v[k++] = lambda0 * E.x;
                v[k++] = lambda0 * E.y;
                v[k++] = lambda1 * E.x;
                v[k++] = lambda1 * E.y;
                v[k++] = lambda2 * E.x;
                v[k++] = lambda2 * E.y;
            }
        }

        R2 B[2] = {T.Edge(1), T.Edge(2)};
        B[0]    = B[0].perp();
        B[1]    = B[1].perp();

        double CK = 0.5;
        for (int p = 0; p < QFK.n; ++p) {
            double w = QFK[p].a * CK;

            std::vector<std::vector<double>> lambda{{-QFK[p].y, QFK[p].x},
                                                    {QFK[p].y, 1 - QFK[p].x},
                                                    {1 - QFK[p].y, QFK[p].x}};

            for (int e = 0; e < 3; ++e) {
                v[k++] = lambda[e][0] * B[0].x * w;
                v[k++] = lambda[e][0] * B[0].y * w;
                v[k++] = lambda[e][1] * B[1].x * w;
                v[k++] = lambda[e][1] * B[1].y * w;
            }
        }
        assert(k == this->ipj_Pi_h.N());
    }

    void FB(const What_d whatd, const Element &K, const Rd &PHat,
            RNMK_ &val) const;
    void FB(const What_d whatd, const int c0, const Element &K, const R2 &PHat,
            RN_ &bfMat);

    bf_type referenceBasisFunction(int i) const {
        switch (i) {
        case 0:
            return [](double *x, int c0) {
                return (c0 == 0) ? 9 * x[0] * (1 - 2 * x[0])
                                 : 3 * x[1] * (4 * x[0] - 1);
            };
            break;
        case 1:
            return [](double *x, int c0) {
                return (c0 == 0) ? 3 * x[0] * (4 * x[1] - 1)
                                 : 9 * x[1] * (1 - 2 * x[1]);
            };
            break;
        case 2:
            return [](double *x, int c0) {
                return (c0 == 0) ? 0.5 * 3 * x[0] * (-x[0] - 6 * x[1] + 2)
                                 : 0.5 * 3 * x[1] * (-6 * x[0] - x[1] + 2);
            };
            break;
        case 4:
            return [](double *x, int c0) {
                return (c0 == 0)
                           ? -18 * x[0] * x[0] - 48 * x[0] * x[1] + 27 * x[0] -
                                 30 * x[1] * x[1] + 36 * x[1] - 9
                           : 3 * x[1] * (4 * x[0] + 4 * x[1] - 3);
            };
            break;
        case 3:
            return [](double *x, int c0) {
                return (c0 == 0) ? -12 * x[0] * x[1] + 3 * x[0] -
                                       30 * x[1] * x[1] + 24 * x[1] - 3
                                 : 9 * x[1] * (2 * x[1] - 1);
            };
            break;
        case 5:
            return [](double *x, int c0) {
                return (c0 == 0) ? -1.5 * x[0] * x[0] + 15 * x[0] * x[1] +
                                       15 * x[1] * x[1] - 15 * x[1] + 1.5
                                 : 1.5 * x[1] * (-6 * x[0] - 5 * x[1] + 4);
            };
            break;
        case 6:
            return [](double *x, int c0) {
                return (c0 == 0)
                           ? 3 * x[0] * (-4 * x[0] - 4 * x[1] + 3)
                           : 30 * x[0] * x[0] + 48 * x[0] * x[1] - 36 * x[0] +
                                 18 * x[1] * x[1] - 27 * x[1] + 9;
            };
            break;
        case 7:
            return [](double *x, int c0) {
                return (c0 == 0) ? 9 * x[0] * (1 - 2 * x[0])
                                 : 30 * x[0] * x[0] + 12 * x[0] * x[1] -
                                       24 * x[0] - 3 * x[1] + 3;
            };
            break;
        case 8:
            return [](double *x, int c0) {
                return (c0 == 0) ? 1.5 * x[0] * (5 * x[0] + 6 * x[1] - 4)
                                 : -15 * x[0] * x[0] - 15 * x[0] * x[1] +
                                       15 * x[0] + 1.5 * x[1] * x[1] - 1.5;
            };
            break;
        case 9:
            return [](double *x, int c0) {
                return (c0 == 0) ? 12 * x[0] * (-x[0] - 4 * x[1] + 1)
                                 : 12 * x[1] * (4 * x[0] + x[1] - 1);
            };
            break;
        case 10:
            return [](double *x, int c0) {
                return (c0 == 0) ? 12 * x[0] * (x[0] + 2 * x[1] - 1)
                                 : 12 * x[1] * (-4 * x[0] - 3 * x[1] + 3);
            };
            break;
        case 11:
            return [](double *x, int c0) {
                return (c0 == 0) ? 12 * x[0] * (-3 * x[0] - 4 * x[1] + 3)
                                 : 12 * x[1] * (2 * x[0] + x[1] - 1);
            };
            break;
        default:
            return [](double *x, int c0) { return 0.; };
            break;
        }
    }
};

int TypeOfFE_BDM2_2d::Data[] = {3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6,  6,  //
                                0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,  2,  //
                                0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3,  3,  //
                                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, //
                                0, 1, 1, 0,                           //
                                0, 0, 12};

void TypeOfFE_BDM2_2d::FB(const What_d whatd, const Element &K, const R2 &PHat,
                          RNMK_ &bfMat) const {

    double x = PHat.x, y = PHat.y;
    assert(bfMat.N() >= 12);
    // assert(0 <= c0 && c0 < 2);
    bfMat = 0.;
    KNM<R> J(2, 2), invJ(2, 2);
    KNM<R> Dphi_ref(2, 2);
    jacobianLinearTransformation<2>(J, K);
    double inv_det_J = inverseDeterminant<2>(J);
    inverseJacobian<2>(J, inv_det_J, invJ);

    double cc = 1.;

    if ((whatd & Fop_D0)) {
        int df = 0;
        int e  = 0;
        { // edge 0
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            int df2             = df + 2;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }
            bf_ref[df0][0] = cc * 9 * x * (1 - 2 * x);
            bf_ref[df0][1] = cc * 3 * y * (4 * x - 1);
            df++;
            bf_ref[df1][0] = cc * 3 * x * (4 * y - 1);
            bf_ref[df1][1] = cc * 9 * y * (1 - 2 * y);
            df++;
            bf_ref[df][0] = cc * 0.5 * 3 * x * (-x - 6 * y + 2);
            bf_ref[df][1] = cc * 0.5 * 3 * y * (-6 * x - y + 2);
            df++;
            e++;
        }

        { // edge 1
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            int df2             = df + 2;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }

            bf_ref[df0][0] =
                cc * (-12 * x * y + 3 * x - 30 * y * y + 24 * y - 3);
            bf_ref[df0][1] = cc * 9 * y * (2 * y - 1);
            df++;
            bf_ref[df1][0] = cc * (-18 * x * x - 48 * x * y + 27 * x -
                                   30 * y * y + 36 * y - 9);
            bf_ref[df1][1] = cc * 3 * y * (4 * x + 4 * y - 3);
            df++;
            bf_ref[df][0] =
                cc * (-1.5 * x * x + 15 * x * y + 15 * y * y - 15 * y + 1.5);
            bf_ref[df][1] = cc * 1.5 * y * (-6 * x - 5 * y + 4);
            df++;
            e++;
        }
        { // edge 2
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }
            bf_ref[df0][0] = cc * 3 * x * (-4 * x - 4 * y + 3);
            bf_ref[df0][1] = cc * (30 * x * x + 48 * x * y - 36 * x +
                                   18 * y * y - 27 * y + 9);
            df++;
            bf_ref[df1][0] = cc * 9 * x * (1 - 2 * x);
            bf_ref[df1][1] =
                cc * (30 * x * x + 12 * x * y - 24 * x - 3 * y + 3);
            df++;
            bf_ref[df][0] = cc * (1.5 * x * (5 * x + 6 * y - 4));
            bf_ref[df][1] =
                cc * (-15 * x * x - 15 * x * y + 15 * x + 1.5 * y * y - 1.5);
            df++;
            e++;
        }
        bf_ref[df][0] = 12 * x * (-x - 4 * y + 1);
        bf_ref[df][1] = 12 * y * (4 * x + y - 1);
        df++;
        bf_ref[df][0] = 12 * x * (x + 2 * y - 1);
        bf_ref[df][1] = 12 * y * (-4 * x - 3 * y + 3);
        df++;
        bf_ref[df][0] = 12 * x * (-3 * x - 4 * y + 3);
        bf_ref[df][1] = 12 * y * (2 * x + y - 1);

        KN_<double> bfMat0 = bfMat('.', 0, op_id);
        piolatTransformation<2>(J(0, '.'), inv_det_J, bf_ref, bfMat0);

        KN_<double> bfMat1 = bfMat('.', 1, op_id);
        piolatTransformation<2>(J(1, '.'), inv_det_J, bf_ref, bfMat1);
    }

    if ((whatd & Fop_D1) || (whatd & Fop_D2)) {
        int df = 0;
        int e  = 0;

        { // edge 0
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }

            Dphi_ref(0, 0) = 9 - 36 * x;
            Dphi_ref(0, 1) = 0;
            Dphi_ref(1, 0) = 12 * y;
            Dphi_ref(1, 1) = 12 * x - 3;
            piolatTransformationGradient<2>(df0, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = 12 * y - 3;
            Dphi_ref(0, 1) = 12 * x;
            Dphi_ref(1, 0) = 0;
            Dphi_ref(1, 1) = 9 - 36 * y;
            piolatTransformationGradient<2>(df1, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = 1.5 * (-2 * x - 6 * y + 2);
            Dphi_ref(0, 1) = -9 * x;
            Dphi_ref(1, 0) = -6 * y;
            Dphi_ref(1, 1) = 1.5 * (-2 * y - 6 * x + 2);
            piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;
            e++;
        }

        { // edge 1
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }
            Dphi_ref(0, 0) = -36 * x - 48 * y + 27;
            Dphi_ref(0, 1) = -48 * x - 60 * y + 36;
            Dphi_ref(1, 0) = 12 * y;
            Dphi_ref(1, 1) = 12 * x + 24 * y - 9;
            piolatTransformationGradient<2>(df0, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = -12 * y + 3;
            Dphi_ref(0, 1) = -12 * x - 60 * y + 24;
            Dphi_ref(1, 0) = 0;
            Dphi_ref(1, 1) = 36 * y - 9;
            piolatTransformationGradient<2>(df1, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = -3 * x + 15 * y;
            Dphi_ref(0, 1) = 15 * x + 30 * y - 15;
            Dphi_ref(1, 0) = -9 * y;
            Dphi_ref(1, 1) = 1.5 * (-6 * x - 10 * y + 4);
            piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;
            e++;
        }
        { // edge 2
            double eOrientation = K.EdgeOrientation(e);
            int df0             = df;
            int df1             = df + 1;
            if (eOrientation < 0) {
                std::swap(df0, df1);
            }
            Dphi_ref(0, 0) = -12 * x;
            Dphi_ref(0, 1) = -12 * x;
            Dphi_ref(1, 0) = 48 * x - 36 * y - 27;
            Dphi_ref(1, 1) = 48 * x - 36 * y - 27;
            piolatTransformationGradient<2>(df0, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = 9 - 36 * x;
            Dphi_ref(0, 1) = 0.;
            Dphi_ref(1, 0) = 60 * x + 12 * y - 24;
            Dphi_ref(1, 1) = 12 * x - 3;
            piolatTransformationGradient<2>(df1, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;

            Dphi_ref(0, 0) = 1.5 * (10 * x + 6 * y - 4);
            Dphi_ref(0, 1) = 9 * x;
            Dphi_ref(1, 0) = -30 * x - 15 * y + 15;
            Dphi_ref(1, 1) = -15 * x + 3 * y;
            piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                            bfMat);
            df++;
            e++;
        }
        Dphi_ref(0, 0) = 12 * (-2 * x - 4 * y + 1);
        Dphi_ref(0, 1) = 12 * (-4 * x);
        Dphi_ref(1, 0) = 12 * (4 * y);
        Dphi_ref(1, 1) = 12 * (2 * y + 4 * x - 1);
        piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                        bfMat);
        df++;

        Dphi_ref(0, 0) = 12 * (2 * x + 2 * x - 1);
        Dphi_ref(0, 1) = 12 * (2 * x);
        Dphi_ref(1, 0) = 12 * (-4 * y);
        Dphi_ref(1, 1) = 12 * (-6 * y - 4 * x + 3);
        piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                        bfMat);
        df++;

        Dphi_ref(0, 0) = 12 * (-6 * x - 4 * y + 3);
        Dphi_ref(0, 1) = 12 * (-4 * x);
        Dphi_ref(1, 0) = 12 * (2 * y);
        Dphi_ref(1, 1) = 12 * (2 * y + 2 * x - 1);
        piolatTransformationGradient<2>(df, J, invJ, inv_det_J, Dphi_ref,
                                        bfMat);
        df++;
        assert(df == 12);
    }
}

// void TypeOfFE_BDM2_2d::FB(const What_d whatd, const int c0, const Element &K,
//                           const R2 &PHat, RN_ &bfMat) {

//    double x = PHat.x, y = PHat.y;
//    assert(bfMat.N() >= 12);
//    assert(0 <= c0 && c0 < 2);
//    bfMat = 0.;
//    KNM<R> J(2, 2), invJ(2, 2);
//    jacobianLinearTransformation<2>(J, K);
//    double inv_det_J = inverseDeterminant<2>(J);
//    inverseJacobian<2>(J, inv_det_J, invJ);

//    std::vector<std::vector<double>> bf_ref;
//    bf_ref.resize(12, std::vector<double>(2));
//    bfMat = 0.;
//    if (whatd == op_id) {
//       int df        = 0;
//       bf_ref[df][0] = 9 * x * (1 - 2 * x);
//       bf_ref[df][1] = 3 * y * (4 * x - 1);
//       df++;
//       bf_ref[df][0] = 3 * x * (4 * y - 1);
//       bf_ref[df][1] = 9 * y * (1 - 2 * y);
//       df++;
//       bf_ref[df][0] = 0.5 * 3 * x * (-x - 6 * y + 2);
//       bf_ref[df][1] = 0.5 * 3 * y * (-6 * x - y + 2);
//       df++;
//       bf_ref[df][0] = -12 * x * y + 3 * x - 30 * y * y + 24 * y - 3;
//       bf_ref[df][1] = 9 * y * (2 * y - 1);
//       df++;
//       bf_ref[df][0] =
//           -18 * x * x - 48 * x * y + 27 * x - 30 * y * y + 36 * y - 9;
//       bf_ref[df][1] = 3 * y * (4 * x + 4 * y - 3);
//       df++;
//       bf_ref[df][0] = -1.5 * x * x + 15 * x * y + 15 * y * y - 15 * y + 1.5;
//       bf_ref[df][1] = 1.5 * y * (-6 * x - 5 * y + 4);
//       df++;
//       bf_ref[df][0] = 3 * x * (-4 * x - 4 * y + 3);
//       bf_ref[df][1] =
//           30 * x * x + 48 * x * y - 36 * x + 18 * y * y - 27 * y + 9;
//       df++;
//       bf_ref[df][0] = 9 * x * (1 - 2 * x);
//       bf_ref[df][1] = 30 * x * x + 12 * x * y - 24 * x - 3 * y + 3;
//       df++;
//       bf_ref[df][0] = 1.5 * x * (5 * x + 6 * y - 4);
//       bf_ref[df][1] = -15 * x * x - 15 * x * y + 15 * x + 1.5 * y * y - 1.5;
//       df++;
//       bf_ref[df][0] = 12 * x * (-x - 4 * y + 1);
//       bf_ref[df][1] = 12 * y * (4 * x + y - 1);
//       df++;
//       bf_ref[df][0] = 12 * x * (x + 2 * y - 1);
//       bf_ref[df][1] = 12 * y * (-4 * x - 3 * y + 3);
//       df++;
//       bf_ref[df][0] = 12 * x * (-3 * x - 4 * y + 3);
//       bf_ref[df][1] = 12 * y * (2 * x + y - 1);

//       piolatTransformation<2>(J(c0, '.'), inv_det_J, bf_ref, bfMat);
//    }
// }

static TypeOfFE_BDM2_2d myBDM2_2d;
GTypeOfFE<Mesh2> &BDM2_2d(myBDM2_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::BDM2 = myBDM2_2d;