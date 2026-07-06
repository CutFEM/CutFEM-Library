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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "QuadratureFormular.hpp"

static const double gauss_n2_1 = (1 - sqrt(1. / 3.)) / 2;
static const double gauss_n2_2 = 1 - gauss_n2_1;

const double gauss_n3_0 = 0.5;
const double gauss_n3_1 = (1 - sqrt(3. / 5.)) / 2;
const double gauss_n3_2 = 1 - gauss_n3_1;

const double pgauss_n3_0 = 8. / 18.;
const double pgauss_n3_1 = 5. / 18.;
const double pgauss_n3_2 = 5. / 18.;

const double pgauss1_n4_a = (18. + sqrt(30.)) / 36.;
const double pgauss1_n4_b = (18. - sqrt(30.)) / 36.;
const double gauss1_n4_a  = sqrt(525. - 70. * sqrt(30.)) / 35.;
const double gauss1_n4_b  = sqrt(525. + 70. * sqrt(30.)) / 35.;

const double pgauss1_n5_0 = 128. / 225.;
const double pgauss1_n5_a = (322 + 13. * sqrt(70.)) / 900.;
const double pgauss1_n5_b = (322 - 13. * sqrt(70.)) / 900.;
const double gauss1_n5_a  = sqrt(245. - 14. * sqrt(70.)) / 21.;
const double gauss1_n5_b  = sqrt(245. + 14. * sqrt(70.)) / 21.;

// Computation of nodes and weights for a Gauss-Legendre quadrature formula on
// [0 1]
typedef GQuadraturePoint<R1> QP1;
QP1 *GaussLegendre(int nn) {
    QP1 *p = new QP1[nn];
    int n  = nn;
    double r, r1, p1, p2, p3, dp3;
    double eps      = 1e-16;
    const double pi = M_PI;

    for (int i = 0, ii = n; i <= (n + 1) / 2 - 1; i++) {
        ii -= 1;
        r = cos(pi * (4 * i + 3) / (4 * n + 2));
        do {
            p2 = 0;
            p3 = 1;
            for (int j = 0; j <= n - 1; j++) {
                p1 = p2;
                p2 = p3;
                p3 = ((2 * j + 1) * r * p2 - j * p1) / (j + 1);
            }

            dp3 = n * (r * p3 - p2) / (r * r - 1);
            r1  = r;
            r   = r - p3 / dp3;
        } while (fabs(r - r1) >= eps * (1 + fabs(r)) * 100);
        p[i].X()  = 0.5 + r / 2;
        p[ii].X() = 0.5 - r / 2;
        p[i].a = p[ii].a = 1. / ((1. - r * r) * dp3 * dp3);
    }
    return p;
}

const QuadratureFormular1d QF_GaussLegendre1(1, QuadratureFormular1d::QP(1, 0.5));
const QuadratureFormular1d QF_GaussLegendre2(3, QuadratureFormular1d::QP(0.5, gauss_n2_1),
                                             QuadratureFormular1d::QP(0.5, gauss_n2_2));
const QuadratureFormular1d QF_GaussLegendre3(5, QuadratureFormular1d::QP(pgauss_n3_0, gauss_n3_0),
                                             QuadratureFormular1d::QP(pgauss_n3_1, gauss_n3_1),
                                             QuadratureFormular1d::QP(pgauss_n3_2, gauss_n3_2));
const QuadratureFormular1d QF_GaussLegendre4(7, QuadratureFormular1d::QP(pgauss1_n4_a / 2., (1. + gauss1_n4_a) / 2.),
                                             QuadratureFormular1d::QP(pgauss1_n4_a / 2., (1. - gauss1_n4_a) / 2.),
                                             QuadratureFormular1d::QP(pgauss1_n4_b / 2., (1. + gauss1_n4_b) / 2.),
                                             QuadratureFormular1d::QP(pgauss1_n4_b / 2., (1. - gauss1_n4_b) / 2.));
const QuadratureFormular1d QF_GaussLegendre5(-1 + 2 * 5, 5, GaussLegendre(5), true);
const QuadratureFormular1d QF_GaussLegendre6(-1 + 2 * 6, 6, GaussLegendre(6), true);
const QuadratureFormular1d QF_GaussLegendre7(-1 + 2 * 7, 7, GaussLegendre(7), true);
const QuadratureFormular1d QF_GaussLegendre8(-1 + 2 * 8, 8, GaussLegendre(8), true);
const QuadratureFormular1d QF_GaussLegendre9(-1 + 2 * 9, 9, GaussLegendre(9), true);
const QuadratureFormular1d QF_GaussLegendre10(-1 + 2 * 10, 10, GaussLegendre(10), true);

// explict instantiation
template <> const GQuadratureFormular<R1> *QF_Simplex<R1>(int exact) {
    switch (exact) {
    case 1:
        return &QF_GaussLegendre1;
    case 2:
        return &QF_GaussLegendre2;
    case 3:
        return &QF_GaussLegendre2;
    case 4:
        return &QF_GaussLegendre3;
    case 5:
        return &QF_GaussLegendre3;
    case 6:
        return &QF_GaussLegendre4;
    case 7:
        return &QF_GaussLegendre4;
    case 8:
        return &QF_GaussLegendre5;
    case 9:
        return &QF_GaussLegendre5;
    case 10:
        return &QF_GaussLegendre6;
    case 11:
        return &QF_GaussLegendre6;
    case 12:
        return &QF_GaussLegendre7;
    case 13:
        return &QF_GaussLegendre7;
    case 14:
        return &QF_GaussLegendre8;
    case 15:
        return &QF_GaussLegendre8;
    case 16:
        return &QF_GaussLegendre9;
    case 17:
        return &QF_GaussLegendre9;
    default:
        return &QF_GaussLegendre9;
    }
};

typedef GQuadraturePoint<R2> QP2;
QP2 *GaussLegendre2D(int exact) {

    const GQuadratureFormular<R1> &quad1D = *QF_Simplex<R1>(exact);
    int n1                                = quad1D.n;
    int n2                                = n1 * n1;
    QP2 *p                                = new QP2[n2];
    int ii                                = 0;
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            p[ii].x = quad1D[i].x;
            p[ii].y = quad1D[j].x;
            p[ii].a = quad1D[i].a * quad1D[j].a;

            ++ii;
        }
    }
    return p;
}

const QuadratureFormular2d QF_GaussLegendreQuad1(1, 1, GaussLegendre2D(1), true);
const QuadratureFormular2d QF_GaussLegendreQuad2(3, 4, GaussLegendre2D(3), true);
const QuadratureFormular2d QF_GaussLegendreQuad3(5, 9, GaussLegendre2D(5), true);
const QuadratureFormular2d QF_GaussLegendreQuad4(7, 16, GaussLegendre2D(7), true);
const QuadratureFormular2d QF_GaussLegendreQuad5(9, 25, GaussLegendre2D(9), true);
const QuadratureFormular2d QF_GaussLegendreQuad9(16, 81, GaussLegendre2D(16), true);

// explict instantiation
const GQuadratureFormular<R2> *QF_Quad(int exact) {
    switch (exact) {
    case 1:
        return &QF_GaussLegendreQuad1;
    case 2:
        return &QF_GaussLegendreQuad2;
    case 3:
        return &QF_GaussLegendreQuad2;
    case 4:
        return &QF_GaussLegendreQuad3;
    case 5:
        return &QF_GaussLegendreQuad3;
    case 6:
        return &QF_GaussLegendreQuad4;
    case 7:
        return &QF_GaussLegendreQuad4;
    case 8:
        return &QF_GaussLegendreQuad5;
    case 9:
        return &QF_GaussLegendreQuad5;
    case 16:
        return &QF_GaussLegendreQuad9;
    // case 10 : return &QF_GaussLegendre6;
    // case 11 : return &QF_GaussLegendre6;
    // case 12 : return &QF_GaussLegendre7;
    // case 13 : return &QF_GaussLegendre7;
    // case 14 : return &QF_GaussLegendre8;
    // case 15 : return &QF_GaussLegendre8;
    // case 17 : return &QF_GaussLegendre9;
    default:
        return &QF_GaussLegendreQuad3;
    }
};

typedef GQuadraturePoint<R3> QP3;
QP3 *GaussLegendre3D(int exact) {

    const GQuadratureFormular<R1> &quad1D = *QF_Simplex<R1>(exact);
    int n1                                = quad1D.n;
    int n3                                = n1 * n1 * n1;
    QP3 *p                                = new QP3[n3];

    int ii = 0;
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            for (int k = 0; k < n1; ++k) {
                p[ii].x = quad1D[i].x;
                p[ii].y = quad1D[j].x;
                p[ii].z = quad1D[k].x;
                p[ii].a = quad1D[i].a * quad1D[j].a * quad1D[k].a;
                ++ii;
            }
        }
    }
    return p;
}

const QuadratureFormular3d QF_GaussLegendreHexa1(1, 1, GaussLegendre3D(1), true);
const QuadratureFormular3d QF_GaussLegendreHexa2(3, 8, GaussLegendre3D(3), true);
const QuadratureFormular3d QF_GaussLegendreHexa3(5, 27, GaussLegendre3D(5), true);
const QuadratureFormular3d QF_GaussLegendreHexa4(7, 64, GaussLegendre3D(7), true);
const QuadratureFormular3d QF_GaussLegendreHexa5(9, 125, GaussLegendre3D(9), true);

// explict instantiation
const GQuadratureFormular<R3> *QF_Hexa(int exact) {
    switch (exact) {
    case 1:
        return &QF_GaussLegendreHexa1;
    case 2:
        return &QF_GaussLegendreHexa2;
    case 3:
        return &QF_GaussLegendreHexa2;
    case 4:
        return &QF_GaussLegendreHexa3;
    case 5:
        return &QF_GaussLegendreHexa3;
    case 6:
        return &QF_GaussLegendreHexa4;
    case 7:
        return &QF_GaussLegendreHexa4;
    case 8:
        return &QF_GaussLegendreHexa5;
    case 9:
        return &QF_GaussLegendreHexa5;
    // case 10 : return &QF_GaussLegendre6;
    // case 11 : return &QF_GaussLegendre6;
    // case 12 : return &QF_GaussLegendre7;
    // case 13 : return &QF_GaussLegendre7;
    // case 14 : return &QF_GaussLegendre8;
    // case 15 : return &QF_GaussLegendre8;
    // case 16 : return &QF_GaussLegendre9;
    // case 17 : return &QF_GaussLegendre9;
    default:
        return &QF_GaussLegendreHexa3;
    }
};

//---------------------------------------------------------------------------------
// Quadrature formula for the time integrals. Those are convenient because
// of the quadrature points at the extremities of the interval
//---------------------------------------------------------------------------------
const double pLob_n3_0 = 1. / 6;
const double pLob_n3_1 = 2. / 3;

const double Lob_n5_2  = 1.4472135954999579 * 0.5;
const double Lob_n5_1  = (-0.4472135954999579 + 1) * 0.5;
const double pLob_n5_1 = 5. / 12;
const double pLob_n5_0 = 1. / 12;

const double Lob_n7_1  = 0.172673164646011;
const double Lob_n7_2  = 0.5;
const double Lob_n7_3  = 0.827326835353989;
const double pLob_n7_0 = 0.05;
const double pLob_n7_1 = 0.272222222222222;
const double pLob_n7_2 = 0.355555555555556;

//---------------------------------------------------------------------------------
const QuadratureFormular1d QF_LumpP1_1D(1, QuadratureFormular1d::QP(0.5, 0.), QuadratureFormular1d::QP(0.5, 1.));
//---------------------------------------------------------------------------------
static GQuadraturePoint<R1> P_QF_Euler[2] = {QuadratureFormular1d::QP(0., R1(0.)),
                                             QuadratureFormular1d::QP(1, R1(1.))};
const QuadratureFormular1d QF_Euler(1, 2, P_QF_Euler);
//---------------------------------------------------------------------------------
// static GQuadraturePoint<R1> P_QF_Lobatto1[2] = {QuadratureFormular1d::QP(0.5, R1(0.)),
//                                                 QuadratureFormular1d::QP(0.5, R1(1.))};
// const QuadratureFormular1d QF_Lobatto1(1, 2, P_QF_Lobatto1);

//---------------------------------------------------------------------------------
static GQuadraturePoint<R1> P_QF_Lobatto2[2] = {
QuadratureFormular1d::QP(0.5000000000000000, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.5000000000000000, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto2(2, 2, P_QF_Lobatto2);

static GQuadraturePoint<R1> P_QF_Lobatto3[3] = {QuadratureFormular1d::QP(pLob_n3_0, R1(0.)),
                                                QuadratureFormular1d::QP(pLob_n3_1, R1(0.5)),
                                                QuadratureFormular1d::QP(pLob_n3_0, R1(1.))};
const QuadratureFormular1d QF_Lobatto3(3, 3, P_QF_Lobatto3);

static GQuadraturePoint<R1> P_QF_Lobatto4[4] = {QuadratureFormular1d::QP(0.08333333333333333, R1(0.)),
                                                QuadratureFormular1d::QP(0.4166666666666667, R1(0.27639320225002106)),
                                                QuadratureFormular1d::QP(0.4166666666666667, R1(0.7236067977499789)),
                                                QuadratureFormular1d::QP(0.08333333333333333, R1(1.))};
const QuadratureFormular1d QF_Lobatto4(4, 4, P_QF_Lobatto4);


//---------------------------------------------------------------------------------
// static GQuadraturePoint<R1> P_QF_Lobatto5[4] = {
//   QuadratureFormular1d::QP(pLob_n5_0,R1(0.)),
//   QuadratureFormular1d::QP(pLob_n5_1,R1(Lob_n5_1)),
//   QuadratureFormular1d::QP(pLob_n5_1,R1(Lob_n5_2)),
//   QuadratureFormular1d::QP(pLob_n5_0,R1(1.))};
// const QuadratureFormular1d QF_Lobatto5(5,4,P_QF_Lobatto5);
static GQuadraturePoint<R1> P_QF_Lobatto5[5] = {
    QuadratureFormular1d::QP(0.05, R1(0.)), QuadratureFormular1d::QP(0.2722222222222222, R1(0.1726731646460114)),
    QuadratureFormular1d::QP(0.35555555555555557, R1(0.5)),
    QuadratureFormular1d::QP(0.2722222222222222, R1(0.8273268353539887)), QuadratureFormular1d::QP(0.05, R1(1.))};
const QuadratureFormular1d QF_Lobatto5(5, 5, P_QF_Lobatto5);

static GQuadraturePoint<R1> P_QF_Lobatto6[6] = {QuadratureFormular1d::QP(0.03333333333333333, R1(0.)),
                                                QuadratureFormular1d::QP(0.1892374781489235, R1(0.11747233803526763)),
                                                QuadratureFormular1d::QP(0.2774291885177432, R1(0.3573842417596774)),
                                                QuadratureFormular1d::QP(0.2774291885177432, R1(0.6426157582403226)),
                                                QuadratureFormular1d::QP(0.1892374781489235, R1(0.8825276619647324)),
                                                QuadratureFormular1d::QP(0.03333333333333333, R1(1.))};
const QuadratureFormular1d QF_Lobatto6(6, 6, P_QF_Lobatto6);

static GQuadraturePoint<R1> P_QF_Lobatto7[7] = {QuadratureFormular1d::QP(0.023809523809523808, R1(0.)),
                                                QuadratureFormular1d::QP(0.13841302368078298, R1(0.08488805186071652)),
                                                QuadratureFormular1d::QP(0.2158726906049313, R1(0.2655756032646429)),
                                                QuadratureFormular1d::QP(0.2438095238095238, R1(0.5)),
                                                QuadratureFormular1d::QP(0.2158726906049313, R1(0.7344243967353571)),
                                                QuadratureFormular1d::QP(0.13841302368078298, R1(0.9151119481392835)),
                                                QuadratureFormular1d::QP(0.023809523809523808, R1(1.))};
const QuadratureFormular1d QF_Lobatto7(7, 7, P_QF_Lobatto7);

static GQuadraturePoint<R1> P_QF_Lobatto8[8] = {QuadratureFormular1d::QP(0.0178571428571429, R1(0.0000000000000000)),
                                                QuadratureFormular1d::QP(0.1053521135717530, R1(0.0641299257451967)),
                                                QuadratureFormular1d::QP(0.1705613462417522, R1(0.2041499092834289)),
                                                QuadratureFormular1d::QP(0.2062293973293519, R1(0.3953503910487606)),
                                                QuadratureFormular1d::QP(0.2062293973293519, R1(0.6046496089512394)),
                                                QuadratureFormular1d::QP(0.1705613462417522, R1(0.7958500907165711)),
                                                QuadratureFormular1d::QP(0.1053521135717530, R1(0.9358700742548033)),
                                                QuadratureFormular1d::QP(0.0178571428571429, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto8(8, 8, P_QF_Lobatto8);

        


static GQuadraturePoint<R1> P_QF_Lobatto9[9] = {QuadratureFormular1d::QP(0.013888888888888888, R1(0.)),
                                                QuadratureFormular1d::QP(0.08274768078040276, R1(0.05012100229426991)),
                                                QuadratureFormular1d::QP(0.13726935625008085, R1(0.16140686024463113)),
                                                QuadratureFormular1d::QP(0.17321425548652317, R1(0.3184412680869109)),
                                                QuadratureFormular1d::QP(0.18575963718820862, R1(0.5)),
                                                QuadratureFormular1d::QP(0.17321425548652317, R1(0.6815587319130891)),
                                                QuadratureFormular1d::QP(0.13726935625008085, R1(0.8385931397553689)),
                                                QuadratureFormular1d::QP(0.08274768078040276, R1(0.94987899770573)),
                                                QuadratureFormular1d::QP(0.013888888888888888, R1(1.))};
const QuadratureFormular1d QF_Lobatto9(9, 9, P_QF_Lobatto9);

static GQuadraturePoint<R1> P_QF_Lobatto10[10] = {
    QuadratureFormular1d::QP(0.011111111111111112, R1(0.)),
    QuadratureFormular1d::QP(0.06665299542553506, R1(0.04023304591677057)),
    QuadratureFormular1d::QP(0.11244467103156322, R1(0.1306130674472475)),
    QuadratureFormular1d::QP(0.1460213418398419, R1(0.26103752509477773)),
    QuadratureFormular1d::QP(0.16376988059194872, R1(0.4173605211668065)),
    QuadratureFormular1d::QP(0.16376988059194872, R1(0.5826394788331936)),
    QuadratureFormular1d::QP(0.1460213418398419, R1(0.7389624749052223)),
    QuadratureFormular1d::QP(0.11244467103156322, R1(0.8693869325527526)),
    QuadratureFormular1d::QP(0.06665299542553506, R1(0.9597669540832294)),
    QuadratureFormular1d::QP(0.011111111111111112, R1(1.))};
const QuadratureFormular1d QF_Lobatto10(10, 10, P_QF_Lobatto10);

static GQuadraturePoint<R1> P_QF_Lobatto11[11] = {
QuadratureFormular1d::QP(0.0090909090909091, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0548061366334974, R1(0.0329992847959704)),
QuadratureFormular1d::QP(0.0935849408901526, R1(0.1077582631684278)),
QuadratureFormular1d::QP(0.1240240521320141, R1(0.2173823365018975)),
QuadratureFormular1d::QP(0.1434395623895040, R1(0.3521209322065303)),
QuadratureFormular1d::QP(0.1501087977278454, R1(0.5000000000000000)),
QuadratureFormular1d::QP(0.1434395623895040, R1(0.6478790677934697)),
QuadratureFormular1d::QP(0.1240240521320141, R1(0.7826176634981026)),
QuadratureFormular1d::QP(0.0935849408901526, R1(0.8922417368315723)),
QuadratureFormular1d::QP(0.0548061366334974, R1(0.9670007152040296)),
QuadratureFormular1d::QP(0.0090909090909091, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto11(11, 11, P_QF_Lobatto11);


static GQuadraturePoint<R1> P_QF_Lobatto12[12] = {
QuadratureFormular1d::QP(0.0075757575757576, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0458422587065980, R1(0.0275503638885589)),
QuadratureFormular1d::QP(0.0789873527821850, R1(0.0903603391779967)),
QuadratureFormular1d::QP(0.1062542088805106, R1(0.1835619234840697)),
QuadratureFormular1d::QP(0.1256378015996006, R1(0.3002345295173255)),
QuadratureFormular1d::QP(0.1357026204553481, R1(0.4317235335725362)),
QuadratureFormular1d::QP(0.1357026204553481, R1(0.5682764664274638)),
QuadratureFormular1d::QP(0.1256378015996006, R1(0.6997654704826745)),
QuadratureFormular1d::QP(0.1062542088805106, R1(0.8164380765159303)),
QuadratureFormular1d::QP(0.0789873527821850, R1(0.9096396608220033)),
QuadratureFormular1d::QP(0.0458422587065980, R1(0.9724496361114411)),
QuadratureFormular1d::QP(0.0075757575757576, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto12(12, 12, P_QF_Lobatto12);

static GQuadraturePoint<R1> P_QF_Lobatto13[13] = {
QuadratureFormular1d::QP(0.0064102564102564, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0389008433734095, R1(0.0233450766789181)),
QuadratureFormular1d::QP(0.0674909633448041, R1(0.0768262176740638)),
QuadratureFormular1d::QP(0.0918234326017750, R1(0.1569057654591213)),
QuadratureFormular1d::QP(0.1103838967830551, R1(0.2585450894543319)),
QuadratureFormular1d::QP(0.1220078951533382, R1(0.3753565349468800)),
QuadratureFormular1d::QP(0.1259654246667234, R1(0.5000000000000000)),
QuadratureFormular1d::QP(0.1220078951533382, R1(0.6246434650531200)),
QuadratureFormular1d::QP(0.1103838967830551, R1(0.7414549105456680)),
QuadratureFormular1d::QP(0.0918234326017750, R1(0.8430942345408787)),
QuadratureFormular1d::QP(0.0674909633448041, R1(0.9231737823259362)),
QuadratureFormular1d::QP(0.0389008433734095, R1(0.9766549233210819)),
QuadratureFormular1d::QP(0.0064102564102564, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto13(13, 13, P_QF_Lobatto13);


static GQuadraturePoint<R1> P_QF_Lobatto14[14] = {
    QuadratureFormular1d::QP(0.00549450549450549451, R1(0.)),
    QuadratureFormular1d::QP(0.03341864224884064232, R1(0.02003247736636954932)), 
    QuadratureFormular1d::QP(0.05829332794935582577, R1(0.06609947308482637450)),
    QuadratureFormular1d::QP(0.08001092588147607121, R1(0.13556570045433692971)), 
    QuadratureFormular1d::QP(0.09741307468670805932, R1(0.22468029853567647234)), 
    QuadratureFormular1d::QP(0.10956312650488537744, R1(0.32863799332864357748)),
    QuadratureFormular1d::QP(0.11580639723422852944, R1(0.44183406555814806617)), 
    QuadratureFormular1d::QP(0.11580639723422852944, R1(0.55816593444185193383)), 
    QuadratureFormular1d::QP(0.10956312650488537744, R1(0.67136200667135642252)),
    QuadratureFormular1d::QP(0.09741307468670805932, R1(0.77531970146432352766)), 
    QuadratureFormular1d::QP(0.08001092588147607121, R1(0.86443429954566307029)), 
    QuadratureFormular1d::QP(0.05829332794935582577, R1(0.93390052691517362550)),
    QuadratureFormular1d::QP(0.03341864224884064232, R1(0.97996752263363045068)), 
    QuadratureFormular1d::QP(0.00549450549450549451, R1(1.))};
const QuadratureFormular1d QF_Lobatto14(14, 14, P_QF_Lobatto14);

// static GQuadraturePoint<R1> P_QF_Lobatto15[9] = {QuadratureFormular1d::QP(0.013888888888889, R1(0.)),
//                                                  QuadratureFormular1d::QP(0.082747680780403, R1(0.050121002294270)),
//                                                  QuadratureFormular1d::QP(0.137269356250081, R1(0.161406860244631)),
//                                                  QuadratureFormular1d::QP(0.173214255486523, R1(0.318441268086911)),
//                                                  QuadratureFormular1d::QP(0.185759637188209, R1(0.500000000000000)),
//                                                  QuadratureFormular1d::QP(0.173214255486523, R1(0.681558731913089)),
//                                                  QuadratureFormular1d::QP(0.137269356250081, R1(0.838593139755369)),
//                                                  QuadratureFormular1d::QP(0.082747680780403, R1(0.949878997705730)),
//                                                  QuadratureFormular1d::QP(0.013888888888889, R1(1.))};
// const QuadratureFormular1d QF_Lobatto15(15, 9, P_QF_Lobatto15);

static GQuadraturePoint<R1> P_QF_Lobatto15[15] = {
QuadratureFormular1d::QP(0.0047619047619048, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0290149465143007, R1(0.0173770367480807)),
QuadratureFormular1d::QP(0.0508300351628590, R1(0.0574589778885118)),
QuadratureFormular1d::QP(0.0702558499012140, R1(0.1182401550240924)),
QuadratureFormular1d::QP(0.0863948236268004, R1(0.1968733972650771)),
QuadratureFormular1d::QP(0.0984936179823067, R1(0.2896809726431637)),
QuadratureFormular1d::QP(0.1059867929634105, R1(0.3923230223181029)),
QuadratureFormular1d::QP(0.1085240581744078, R1(0.5000000000000000)),
QuadratureFormular1d::QP(0.1059867929634105, R1(0.6076769776818971)),
QuadratureFormular1d::QP(0.0984936179823067, R1(0.7103190273568363)),
QuadratureFormular1d::QP(0.0863948236268004, R1(0.8031266027349229)),
QuadratureFormular1d::QP(0.0702558499012140, R1(0.8817598449759076)),
QuadratureFormular1d::QP(0.0508300351628590, R1(0.9425410221114882)),
QuadratureFormular1d::QP(0.0290149465143007, R1(0.9826229632519192)),
QuadratureFormular1d::QP(0.0047619047619048, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto15(15, 15, P_QF_Lobatto15);


static GQuadraturePoint<R1> P_QF_Lobatto16[16] = {
QuadratureFormular1d::QP(0.0041666666666667, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0254251805029600, R1(0.0152159768648911)),
QuadratureFormular1d::QP(0.0446968486629655, R1(0.0503997334532639)),
QuadratureFormular1d::QP(0.0621276910662570, R1(0.1039958540690925)),
QuadratureFormular1d::QP(0.0770134904035821, R1(0.1738056485587535)),
QuadratureFormular1d::QP(0.0887459566958521, R1(0.2569702890564312)),
QuadratureFormular1d::QP(0.0968450119126018, R1(0.3500847655496184)),
QuadratureFormular1d::QP(0.1009791540891149, R1(0.4493368632390253)),
QuadratureFormular1d::QP(0.1009791540891149, R1(0.5506631367609747)),
QuadratureFormular1d::QP(0.0968450119126018, R1(0.6499152344503816)),
QuadratureFormular1d::QP(0.0887459566958521, R1(0.7430297109435688)),
QuadratureFormular1d::QP(0.0770134904035821, R1(0.8261943514412465)),
QuadratureFormular1d::QP(0.0621276910662570, R1(0.8960041459309076)),
QuadratureFormular1d::QP(0.0446968486629655, R1(0.9496002665467360)),
QuadratureFormular1d::QP(0.0254251805029600, R1(0.9847840231351089)),
QuadratureFormular1d::QP(0.0041666666666667, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto16(16, 16, P_QF_Lobatto16);

static GQuadraturePoint<R1> P_QF_Lobatto17[17] = {
QuadratureFormular1d::QP(0.0036764705882353, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0224609702716271, R1(0.0134339116842909)),
QuadratureFormular1d::QP(0.0395991352518435, R1(0.0445600020422132)),
QuadratureFormular1d::QP(0.0552964545035141, R1(0.0921518743891149)),
QuadratureFormular1d::QP(0.0689938731009633, R1(0.1544855096861577)),
QuadratureFormular1d::QP(0.0801973309988108, R1(0.2293073003349492)),
QuadratureFormular1d::QP(0.0885021267578289, R1(0.3139127832172615)),
QuadratureFormular1d::QP(0.0936081698388096, R1(0.4052440132408413)),
QuadratureFormular1d::QP(0.0953309373767347, R1(0.5000000000000000)),
QuadratureFormular1d::QP(0.0936081698388096, R1(0.5947559867591588)),
QuadratureFormular1d::QP(0.0885021267578289, R1(0.6860872167827385)),
QuadratureFormular1d::QP(0.0801973309988108, R1(0.7706926996650507)),
QuadratureFormular1d::QP(0.0689938731009633, R1(0.8455144903138423)),
QuadratureFormular1d::QP(0.0552964545035141, R1(0.9078481256108851)),
QuadratureFormular1d::QP(0.0395991352518435, R1(0.9554399979577868)),
QuadratureFormular1d::QP(0.0224609702716271, R1(0.9865660883157091)),
QuadratureFormular1d::QP(0.0036764705882353, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto17(17, 17, P_QF_Lobatto17);


static GQuadraturePoint<R1> P_QF_Lobatto18[18] = {
QuadratureFormular1d::QP(0.0032679738562092, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0199853144054570, R1(0.0119472212939007)),
QuadratureFormular1d::QP(0.0353185834428169, R1(0.0396754073262330)),
QuadratureFormular1d::QP(0.0495081358587514, R1(0.0822032323909549)),
QuadratureFormular1d::QP(0.0621052665664836, R1(0.1381603353583787)),
QuadratureFormular1d::QP(0.0727059807869011, R1(0.2057475828406691)),
QuadratureFormular1d::QP(0.0809697586188012, R1(0.2827924815439380)),
QuadratureFormular1d::QP(0.0866310547447282, R1(0.3668186735608595)),
QuadratureFormular1d::QP(0.0895079317198515, R1(0.4551254532576740)),
QuadratureFormular1d::QP(0.0895079317198515, R1(0.5448745467423260)),
QuadratureFormular1d::QP(0.0866310547447282, R1(0.6331813264391405)),
QuadratureFormular1d::QP(0.0809697586188012, R1(0.7172075184560620)),
QuadratureFormular1d::QP(0.0727059807869011, R1(0.7942524171593308)),
QuadratureFormular1d::QP(0.0621052665664836, R1(0.8618396646416213)),
QuadratureFormular1d::QP(0.0495081358587514, R1(0.9177967676090450)),
QuadratureFormular1d::QP(0.0353185834428169, R1(0.9603245926737669)),
QuadratureFormular1d::QP(0.0199853144054570, R1(0.9880527787060993)),
QuadratureFormular1d::QP(0.0032679738562092, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto18(18, 18, P_QF_Lobatto18);

static GQuadraturePoint<R1> P_QF_Lobatto19[19] = {
QuadratureFormular1d::QP(0.0029239766081871, R1(0.0000000000000000)),
QuadratureFormular1d::QP(0.0178966825930883, R1(0.0106941168889599)),
QuadratureFormular1d::QP(0.0316909458813149, R1(0.0355492359237068)),
QuadratureFormular1d::QP(0.0445658785496035, R1(0.0737697111016770)),
QuadratureFormular1d::QP(0.0561576707386525, R1(0.1242528987236935)),
QuadratureFormular1d::QP(0.0661336402243754, R1(0.1855459313673897)),
QuadratureFormular1d::QP(0.0742069712979695, R1(0.2558853571596432)),
QuadratureFormular1d::QP(0.0801454620220306, R1(0.3332475760877507)),
QuadratureFormular1d::QP(0.0837782922635714, R1(0.4154069882953592)),
QuadratureFormular1d::QP(0.0850009596424136, R1(0.5000000000000000)),
QuadratureFormular1d::QP(0.0837782922635714, R1(0.5845930117046407)),
QuadratureFormular1d::QP(0.0801454620220306, R1(0.6667524239122493)),
QuadratureFormular1d::QP(0.0742069712979695, R1(0.7441146428403568)),
QuadratureFormular1d::QP(0.0661336402243754, R1(0.8144540686326103)),
QuadratureFormular1d::QP(0.0561576707386525, R1(0.8757471012763065)),
QuadratureFormular1d::QP(0.0445658785496035, R1(0.9262302888983230)),
QuadratureFormular1d::QP(0.0316909458813149, R1(0.9644507640762932)),
QuadratureFormular1d::QP(0.0178966825930883, R1(0.9893058831110401)),
QuadratureFormular1d::QP(0.0029239766081871, R1(1.0000000000000000))};
const QuadratureFormular1d QF_Lobatto19(19, 19, P_QF_Lobatto19);

static GQuadraturePoint<R1> P_QF_Lobatto20[20] = {
    QuadratureFormular1d::QP(0.00263157894736842, R1(0.)),
    QuadratureFormular1d::QP(0.0161185615942445, R1(0.00962814755304292)), 
    QuadratureFormular1d::QP(0.0285909010637834, R1(0.0320327505936673)),
    QuadratureFormular1d::QP(0.0403158819980598, R1(0.0665610109550249)), 
    QuadratureFormular1d::QP(0.0509957498497254, R1(0.112315869523972)), 
    QuadratureFormular1d::QP(0.0603546138143374, R1(0.168111798854844)),
    QuadratureFormular1d::QP(0.0681502411793621, R1(0.232503567984057)), 
    QuadratureFormular1d::QP(0.0741807770354584, R1(0.303823408143045)), 
    QuadratureFormular1d::QP(0.0782900513237377, R1(0.380224147038507)),
    QuadratureFormular1d::QP(0.0803716431939229, R1(0.459727031380589)), 
    QuadratureFormular1d::QP(0.0803716431939229, R1(0.540272968619411)), 
    QuadratureFormular1d::QP(0.0782900513237377, R1(0.619775852961493)),
    QuadratureFormular1d::QP(0.0741807770354584, R1(0.696176591856955)), 
    QuadratureFormular1d::QP(0.0681502411793621, R1(0.767496432015943)),
    QuadratureFormular1d::QP(0.0603546138143374, R1(0.831888201145156)),
    QuadratureFormular1d::QP(0.0509957498497254, R1(0.887684130476028)),
    QuadratureFormular1d::QP(0.0403158819980598, R1(0.933438989044975)),
    QuadratureFormular1d::QP(0.0285909010637834, R1(0.967967249406333)),
    QuadratureFormular1d::QP(0.0161185615942445, R1(0.990371852446957)),
    QuadratureFormular1d::QP(0.00263157894736842, R1(1.))};
const QuadratureFormular1d QF_Lobatto20(20, 20, P_QF_Lobatto20);


// explict instantiation
int exactLobatto_nPt(int n) {
    switch (n) {
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 3;
    case 4:
        return 4;
    case 5:
        return 5;
    case 6:
        return 6;
    case 7:
        return 7;
    case 8:
        return 8;
    case 9:
        return 9;
    case 10:
        return 10;
    default: {
        std::cout << " Lobatto with " << n << " not implemented " << std::endl;
        exit(EXIT_FAILURE);
        return 0;
    }
    }
};

// explict instantiation
const QuadratureFormular1d *Lobatto(int exact) {

    if ((exact < 2) || (exact > 20))   
        std::cout << "PROBLEM: the number of quadrature points must be between 2 and 20, choosing N=3\n";

    switch (exact) {
    // case 1:
    //     return &QF_Lobatto1;
    case 2:
        return &QF_Lobatto2;
    case 3:
        return &QF_Lobatto3;
    case 4:
        return &QF_Lobatto4;
    case 5:
        return &QF_Lobatto5;
    case 6:
        return &QF_Lobatto6;
    case 7:
        return &QF_Lobatto7;
    case 8:
        return &QF_Lobatto8;
    case 9:
        return &QF_Lobatto9;
    case 10:
        return &QF_Lobatto10;
    case 11:
        return &QF_Lobatto11;
    case 12:
        return &QF_Lobatto12;
    case 13:
        return &QF_Lobatto13;
    case 14:
        return &QF_Lobatto14;
    case 15:
        return &QF_Lobatto15;
    case 16:
        return &QF_Lobatto16;
    case 17:
        return &QF_Lobatto17;
    case 18:
        return &QF_Lobatto18;
    case 19:
        return &QF_Lobatto19;
    case 20:
        return &QF_Lobatto20;
    // case 15:
    //    return &QF_Lobatto15;
    default:
        return &QF_Lobatto3;
    }
};

// template<class QF,int ON>
//   QF  * QF_exact(int exact,QF * p=0)
//   {
//     exact=max(0,exact);
//     const int N=100;
//     assert(exact<N&& exact>=0);
//     static QF ** a=0;
//     if(a==0)
//       { //
// 	a = new  QF*[N];
// 	assert(a);
// 	for(int i=0;i<N;++i)
// 	  a[i]=0;
//       }
//     assert(a && exact >=0 && exact < N);
//     if ( p  )
//       {
// 	//cout << endl << " QF " << exact << " " << p->exact << " " << p->n <<
// endl;; 	for( int i=0;i<=exact;++i)
// 	  {
// 	    if( a[i]== 0 || a[i]->n > p->n)
// 	      a[i]= p;
// 	    //  cout << " QF: on " << ON << " exact P_" << i << " : "<< a[i]->n
// << endl;
// 	  }
//       }
//     else
//       p=a[exact];
//     return p;
//   }

// template<class Rd>
// GQuadratureFormular<Rd> * QF_Simplex(int exact)
//  {
//    return  QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact);
//  }

//---------------------------------------------------------------------------------
static GQuadraturePoint<R2> P_QF_NODE_TRIANGLE[3] = {QuadratureFormular2d::QP(1, R2(0., 0.)),
                                                     QuadratureFormular2d::QP(1, R2(1., 0.)),
                                                     QuadratureFormular2d::QP(1, R2(0., 1.))};
const QuadratureFormular2d QF_NODE_TRIANGLE(0, 3, P_QF_NODE_TRIANGLE);

// ----------------------------------------------------------------------
static GQuadraturePoint<R2> P_QuadratureFormular_T_1[1] = {GQuadraturePoint<R2>(1., R2(1. / 3., 1. / 3.))};

const GQuadratureFormular<R2> QuadratureFormular_T_1(1, 1, P_QuadratureFormular_T_1);
// ----------------------------------------------------------------------

static GQuadraturePoint<R2> P_QuadratureFormular_T_2[3] = {GQuadraturePoint<R2>(1. / 3., R2(0.5, 0.5)),
                                                           GQuadraturePoint<R2>(1. / 3., R2(0.0, 0.5)),
                                                           GQuadraturePoint<R2>(1. / 3., R2(0.5, 0.0))};

GQuadratureFormular<R2> const QuadratureFormular_T_2(2, 3, P_QuadratureFormular_T_2);
// ----------------------------------------------------------------------

const double sqrt15 = 3.87298334620741688517926539978;
const double t_T5 = 1.E0 / 3.E0, A_T5 = 0.225E0;
const double r_T5 = (6 - sqrt15) / 21, s_T5 = (9 + 2 * sqrt15) / 21, B_T5 = (155 - sqrt15) / 1200;
const double u_T5 = (6 + sqrt15) / 21, v_T5 = (9 - 2 * sqrt15) / 21, C_T5 = (155 + sqrt15) / 1200;

static GQuadraturePoint<R2> P_QuadratureFormular_T_5[] = {
    GQuadraturePoint<R2>(A_T5, R2(t_T5, t_T5)), GQuadraturePoint<R2>(B_T5, R2(r_T5, r_T5)),
    GQuadraturePoint<R2>(B_T5, R2(r_T5, s_T5)), GQuadraturePoint<R2>(B_T5, R2(s_T5, r_T5)),
    GQuadraturePoint<R2>(C_T5, R2(u_T5, u_T5)), GQuadraturePoint<R2>(C_T5, R2(u_T5, v_T5)),
    GQuadraturePoint<R2>(C_T5, R2(v_T5, u_T5))};
const GQuadratureFormular<R2> QuadratureFormular_T_5(5, 7, P_QuadratureFormular_T_5);

// ----------------------------------------------------------------------
static GQuadraturePoint<R2> P_QuadratureFormular_T_7[] = {
    GQuadraturePoint<R2>(0.0102558174092 / 2, R2(1.0000000000000, 0.0000000000000)),
    GQuadraturePoint<R2>(0.0102558174092 / 2, R2(0.0000000000000, 0.0000000000000)),
    GQuadraturePoint<R2>(0.0102558174092 / 2, R2(0.0000000000000, 1.0000000000000)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.7839656651012, 0.0421382841642)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.1738960507345, 0.7839656651012)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.1738960507345, 0.0421382841642)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.0421382841642, 0.1738960507345)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.7839656651012, 0.1738960507345)),
    GQuadraturePoint<R2>(0.1116047046647 / 2, R2(0.0421382841642, 0.7839656651012)),
    GQuadraturePoint<R2>(0.1679775595335 / 2, R2(0.4743880861752, 0.4743880861752)),
    GQuadraturePoint<R2>(0.1679775595335 / 2, R2(0.4743880861752, 0.0512238276497)),
    GQuadraturePoint<R2>(0.1679775595335 / 2, R2(0.0512238276497, 0.4743880861752)),
    GQuadraturePoint<R2>(0.2652238803946 / 2, R2(0.2385615300181, 0.5228769399639)),
    GQuadraturePoint<R2>(0.2652238803946 / 2, R2(0.5228769399639, 0.2385615300181)),
    GQuadraturePoint<R2>(0.2652238803946 / 2, R2(0.2385615300181, 0.2385615300181))};
const GQuadratureFormular<R2> QuadratureFormular_T_7(7, 15, P_QuadratureFormular_T_7);

// awk '/21:/,/28:/ {print "GQuadraturePoint<R2>(" $3 "/2," $1"," $2"),"}'
// coords.txt
static GQuadraturePoint<R2> P_QuadratureFormular_T_9[] = {
    GQuadraturePoint<R2>(0.0519871420646 / 2, R2(0.0451890097844, 0.0451890097844)),
    GQuadraturePoint<R2>(0.0519871420646 / 2, R2(0.0451890097844, 0.9096219804312)),
    GQuadraturePoint<R2>(0.0519871420646 / 2, R2(0.9096219804312, 0.0451890097844)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.7475124727339, 0.0304243617288)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.2220631655373, 0.0304243617288)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.7475124727339, 0.2220631655373)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.2220631655373, 0.7475124727339)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.0304243617288, 0.7475124727339)),
    GQuadraturePoint<R2>(0.0707034101784 / 2, R2(0.0304243617288, 0.2220631655373)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.1369912012649, 0.2182900709714)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.6447187277637, 0.2182900709714)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.1369912012649, 0.6447187277637)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.2182900709714, 0.6447187277637)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.2182900709714, 0.1369912012649)),
    GQuadraturePoint<R2>(0.0909390760952 / 2, R2(0.6447187277637, 0.1369912012649)),
    GQuadraturePoint<R2>(0.1032344051380 / 2, R2(0.0369603304334, 0.4815198347833)),
    GQuadraturePoint<R2>(0.1032344051380 / 2, R2(0.4815198347833, 0.0369603304334)),
    GQuadraturePoint<R2>(0.1032344051380 / 2, R2(0.4815198347833, 0.4815198347833)),
    GQuadraturePoint<R2>(0.1881601469167 / 2, R2(0.4036039798179, 0.1927920403641)),
    GQuadraturePoint<R2>(0.1881601469167 / 2, R2(0.4036039798179, 0.4036039798179)),
    GQuadraturePoint<R2>(0.1881601469167 / 2, R2(0.1927920403641, 0.4036039798179))};
const GQuadratureFormular<R2> QuadratureFormular_T_9(9, 21, P_QuadratureFormular_T_9);

template <> const GQuadratureFormular<R2> *QF_Simplex<R2>(int exact) {
    switch (exact) {
    case 0:
        return &QF_NODE_TRIANGLE;
    case 1:
        return &QuadratureFormular_T_1;
    case 2:
        return &QuadratureFormular_T_2;
    case 3:
        return &QuadratureFormular_T_5;
    case 4:
        return &QuadratureFormular_T_5;
    case 5:
        return &QuadratureFormular_T_5;
    case 6:
        return &QuadratureFormular_T_7;
    case 7:
        return &QuadratureFormular_T_7;
    case 8:
        return &QuadratureFormular_T_9;
    case 9:
        return &QuadratureFormular_T_9;
    default:
        return &QuadratureFormular_T_9;
    }
}
template <> const GQuadratureFormular<R3> *QF_Simplex<R3>(int exact);

/*
void QuadratureFormular1d::Check()
{
    QF_exact<QuadratureFormular1d,2>(exact,this);
    int err=0;
    for(int m=0;m<=exact;++m)
    {
        double ve = 1./ ( m+1), v=0.;
        for (int i=0;i<n;++i)
            v += p[i].a*pow(p[i].x,m);
        if (abs( ve-v)/ve > 1.e-8)
        {
            cout << " erreur QuadratureFormular1d  n= " << n  << " exact = " <<
exact << endl; cout << " int x^" <<m << " == " << ve << " ! = " << v  << endl;
            err++;
        }
    }
    assert(err==0);
}
*/

/*
Region: Simplex
Dimension: 3
Degree: 3
Points: 5
Structure: Fully symmetric
Rule struct: 0 0 0 0 0 0 1 1 0 0 0
Generator: [ Fully symmetric ]
( 0.25, 0.25, 0.25, )
Corresponding weight:
-0.133333333333333333333333333333333,
Generator: [ Fully symmetric ]
( 0.166666666666666666666666666666666,
0.166666666666666666666666666666666,
0.166666666666666666666666666666666,
)
Corresponding weight:
0.075,

Region: Simplex
Dimension: 3
Degree: 4
Points: 11
Structure: Fully symmetric
Rule struct: 0 0 0 0 0 0 1 1 1 0 0
Generator: [ Fully symmetric ]
( 0.25, 0.25, 0.25, )
Corresponding weight:
-0.0131555555555555555555555555555555,
Generator: [ Fully symmetric ]
( 0.0714285714285714285714285714285714,
0.0714285714285714285714285714285714,
0.0714285714285714285714285714285714,
)
Corresponding weight:
10 ^ -3 x 7.62222222222222222222222222222222,

Generator: [ Fully symmetric ]
( 0.399403576166799204996102147461640,
0.399403576166799204996102147461640,
0.100596423833200795003897852538359,
)
Corresponding weight:
0.0248888888888888888888888888888888,


Region: Simplex
Dimension: 3
Degree: 6
Points: 24
Structure: Fully symmetric
Rule struct: 0 0 0 0 0 0 0 3 0 1 0
Generator: [ Fully symmetric ]
( 0.214602871259152029288839219386284,
0.214602871259152029288839219386284,
0.214602871259152029288839219386284,
)
Corresponding weight:
10 ^ -3 x 6.65379170969458201661510459291332,
Generator: [ Fully symmetric ]
( 0.0406739585346113531155794489564100,
0.0406739585346113531155794489564100,
0.0406739585346113531155794489564100,
)
Corresponding weight:
10 ^ -3 x 1.67953517588677382466887290765614,

Generator: [ Fully symmetric ]
( 0.322337890142275510343994470762492,
0.322337890142275510343994470762492,
0.322337890142275510343994470762492,
)
Corresponding weight:
10 ^ -3 x 9.22619692394245368252554630895433,

Generator: [ Fully symmetric ]
( 0.0636610018750175252992355276057269,
0.0636610018750175252992355276057269,
0.269672331458315808034097805727606,
)
Corresponding weight:
10 ^ -3 x 8.03571428571428571428571428571428,



Region: Simplex
Dimension: 3
Degree: 7
Points: 31
Structure: Fully symmetric
Rule struct: 0 1 0 0 0 0 1 3 0 1 0
Generator: [ Fully symmetric ]
( 0.5, 0.5, 0., )
Corresponding weight:
10 ^ -4 x 9.70017636684303350970017636684303,
Generator: [ Fully symmetric ]
( 0.25, 0.25, 0.25, )
Corresponding weight:
0.0182642234661088202912015685649462,

Generator: [ Fully symmetric ]
( 0.0782131923303180643739942508375545,
0.0782131923303180643739942508375545,
0.0782131923303180643739942508375545,
)
Corresponding weight:
0.0105999415244136869164138748545257,

Generator: [ Fully symmetric ]
( 0.121843216663905174652156372684818,
0.121843216663905174652156372684818,
0.121843216663905174652156372684818,
)
Corresponding weight:
-0.0625177401143318516914703474927900,

Generator: [ Fully symmetric ]
( 0.332539164446420624152923823157707,
0.332539164446420624152923823157707,
0.332539164446420624152923823157707,
)
Corresponding weight:
10 ^ -3 x 4.89142526307349938479576303671027,

Generator: [ Fully symmetric ]
( 0.1, 0.1, 0.2, )
Corresponding weight:
0.0275573192239858906525573192239858,



Region: Simplex
Dimension: 3
Degree: 8
Points: 43
Structure: Fully symmetric
Rule struct: 0 0 0 0 0 0 1 3 1 2 0
Generator: [ Fully symmetric ]
( 0.25, 0.25, 0.25, )
Corresponding weight:
-0.0205001886586399158405865177642941,
Generator: [ Fully symmetric ]
( 0.206829931610673204083980900024961,
0.206829931610673204083980900024961,
0.206829931610673204083980900024961,
)
Corresponding weight:
0.0142503058228669012484397415358704,

Generator: [ Fully symmetric ]
( 0.0821035883105467230906058078714215,
0.0821035883105467230906058078714215,
0.0821035883105467230906058078714215,
)
Corresponding weight:
10 ^ -3 x 1.96703331313390098756280342445466,

Generator: [ Fully symmetric ]
( 10 ^ -3 x 5.78195050519799725317663886414270,
10 ^ -3 x 5.78195050519799725317663886414270,
10 ^ -3 x 5.78195050519799725317663886414270,
)
Corresponding weight:
10 ^ -4 x 1.69834109092887379837744566704016,

Generator: [ Fully symmetric ]
( 0.0505327400188942244256245285579071,
0.0505327400188942244256245285579071,
0.449467259981105775574375471442092,
)
Corresponding weight:
10 ^ -3 x 4.57968382446728180074351446297276,

Generator: [ Fully symmetric ]
( 0.229066536116811139600408854554753,
0.229066536116811139600408854554753,
0.0356395827885340437169173969506114,
)
Corresponding weight:
10 ^ -3 x 5.70448580868191850680255862783040,

Generator: [ Fully symmetric ]
( 0.0366077495531974236787738546327104,
0.0366077495531974236787738546327104,
0.190486041934633455699433285315099,
)
Corresponding weight:
10 ^ -3 x 2.14051914116209259648335300092023,

Region: Simplex
Dimension: 3
Degree: 9
Points: 53
Structure: Fully symmetric
Rule struct: 0 0 0 0 0 0 1 4 0 3 0
Generator: [ Fully symmetric ]
( 0.25, 0.25, 0.25, )
Corresponding weight:
-0.137799038326108641245781827743502,
Generator: [ Fully symmetric ]
( 0.0483510385497367408791638973335652,
0.0483510385497367408791638973335652,
0.0483510385497367408791638973335652,
)
Corresponding weight:
10 ^ -3 x 1.86533656908528954751732956352791,

Generator: [ Fully symmetric ]
( 0.324579280117882365858796772348839,
0.324579280117882365858796772348839,
0.324579280117882365858796772348839,
)
Corresponding weight:
10 ^ -3 x 4.30942396949340069685481460143672,

Generator: [ Fully symmetric ]
( 0.114616540223995219683790300317451,
0.114616540223995219683790300317451,
0.114616540223995219683790300317451,
)
Corresponding weight:
-0.0901847664812015251273931145424214,

Generator: [ Fully symmetric ]
( 0.225489951911513918847703826550363,
0.225489951911513918847703826550363,
0.225489951911513918847703826550363,
)
Corresponding weight:
0.0446725762025114446937850368088403,

Generator: [ Fully symmetric ]
( 0.131627809246869809838976179197822,
0.131627809246869809838976179197822,
0.0836647016171849679665385541714202,
)
Corresponding weight:
0.0347004058845507618227002576570271,

Generator: [ Fully symmetric ]
( 0.433951461411406772411426685800524,
0.433951461411406772411426685800524,
0.107769859549428611825218473500043,
)
Corresponding weight:
10 ^ -3 x 3.35258390266064697010720483063291,

Generator: [ Fully symmetric ]
( 10 ^ -3 x -1.37627731813820071002030321419028,
10 ^ -3 x -1.37627731813820071002030321419028,
0.276553472636807342120192082186261,
)
Corresponding weight:
10 ^ -4 x 4.31628875556996929641889902726182,


 */
// template<class Rd>
// void GQuadratureFormular<Rd>::Verification()
// {
//   // QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact,this);
//   // const int d=Rd::d;
//   // const double tol=1.e-12;
//   // double err=0;
//   // double a[d+1],h;
//
//   // for (int k=0;k<=exact;  k++)
//   //   {
//
//   //     //  formule magic: int_K  \lamda^k =   d! k! / ( d+k)!
//
//   //     double sa[d+1];
//   //     for(int l=0;l<=d;++l)
//   // 	sa[l]=0.;
//
//   //     for (int j=0;j<n;j++)
//   // 	{
//   // 	  h = p[j];
//   // 	  Rd P = p[j];
//
//   // 	  for(int l=0;l<d;++l)
//   // 	    a[l]=p[j][l];
//   // 	  a[d] = 1.-p[j].sum();
//
//   // 	  for(int l=0;l<=d;++l)
//   // 	    sa[l]+=h*pow(a[l],k);
//   // 	}
//
//   //     double se(1),see(1);
//   //     for (int i=1;i<=k;i++)
//   // 	se *= (double) i / (double) (i+d);
//   //     see=se;
//
//   //     for(int l=0;l<=d;++l)
//   // 	err = Max(err,Abs(se-sa[l]));
//
//   //     if (err>tol)
//   // 	{
//   // 	  cerr << " d= " << d << "T Ordre= " << k << " d!k!/(d+k)!= " <<
//   se
//   << " " ;
//   // 	  for(int l=0;l<=d;++l)
//   // 	    cerr << sa[l] << " ";
//   // 	  cerr << " err= " << err << endl;
//   // 	}
//   //   }
//
//   // if(err>tol)
//   //   {
//   //     cerr << "Erreur dans la formule d'integration d=" <<d  << " exact =
//   " << exact
//   // 	   << " Nb Point = " << n << endl;
//   //     assert(0);
//   //   }
//
// }
/*
  from:
http://www.cs.kuleuven.be/~nines/research/ecf/mtables.html
 */
typedef GQuadraturePoint<R3> PQP3;
typedef GQuadratureFormular<R3> PQF3;
// 2 4 (formule 2)
// A.H. Stroud, Approximate calculation of multiple integrals, Prentice-Hall,
// Englewood Cliffs, N.J., 1971. JavaScript:formule('t3-2-4b')

PQP3 QF_TET_1[] = {PQP3(R3(0.25, 0.25, 0.25), 1)};
PQF3 const QuadratureFormular_Tet_1(1, 1, QF_TET_1);

PQP3 QF_TET_2[] = {PQP3(R3(0.58541019662496845446137605030968, 0.138196601125010515179541316563436,
                           0.138196601125010515179541316563436),
                        0.25),
                   PQP3(R3(0.138196601125010515179541316563436, 0.58541019662496845446137605030968,
                           0.138196601125010515179541316563436),
                        0.25),
                   PQP3(R3(0.138196601125010515179541316563436, 0.138196601125010515179541316563436,
                           0.58541019662496845446137605030968),
                        0.25),
                   PQP3(R3(0.138196601125010515179541316563436, 0.138196601125010515179541316563436,
                           0.138196601125010515179541316563436),
                        0.25)};
PQF3 const QuadratureFormular_Tet_2(2, 4, QF_TET_2);
// 5  14  (formule 1)
/*
GM78
A. Grundmann and H.M. M�ller, Invariant integration formulas for the n-simplex
by combinatorial methods, SIAM J. Numer. Anal. 15 (1978), 282--290.
 */
PQP3 QF_TET_5[] = {PQP3(R3(0.7217942490673263207930282587889082, 0.0927352503108912264023239137370306,
                           0.0927352503108912264023239137370306),
                        0.0122488405193936582572850342477212 * 6.),
                   PQP3(R3(0.0927352503108912264023239137370306, 0.7217942490673263207930282587889082,
                           0.0927352503108912264023239137370306),
                        0.0122488405193936582572850342477212 * 6.),
                   PQP3(R3(0.0927352503108912264023239137370306, 0.0927352503108912264023239137370306,
                           0.7217942490673263207930282587889082),
                        0.0122488405193936582572850342477212 * 6.),
                   PQP3(R3(0.0927352503108912264023239137370306, 0.0927352503108912264023239137370306,
                           0.0927352503108912264023239137370306),
                        0.0122488405193936582572850342477212 * 6),
                   PQP3(R3(0.067342242210098170607962798709629, 0.310885919263300609797345733763457,
                           0.310885919263300609797345733763457),
                        0.0187813209530026417998642753888810 * 6.),
                   PQP3(R3(0.310885919263300609797345733763457, 0.067342242210098170607962798709629,
                           0.310885919263300609797345733763457),
                        0.0187813209530026417998642753888810 * 6.),
                   PQP3(R3(0.310885919263300609797345733763457, 0.310885919263300609797345733763457,
                           0.067342242210098170607962798709629),
                        0.0187813209530026417998642753888810 * 6.),
                   PQP3(R3(0.310885919263300609797345733763457, 0.310885919263300609797345733763457,
                           0.310885919263300609797345733763457),
                        0.0187813209530026417998642753888810 * 6.),
                   PQP3(R3(0.454496295874350350508119473720660, 0.454496295874350350508119473720660,
                           0.045503704125649649491880526279339),
                        7.09100346284691107301157135337624E-3 * 6.),
                   PQP3(R3(00.454496295874350350508119473720660, 0.045503704125649649491880526279339,
                           0.454496295874350350508119473720660),
                        7.09100346284691107301157135337624E-3 * 6.),
                   PQP3(R3(00.045503704125649649491880526279339, 0.454496295874350350508119473720660,
                           0.454496295874350350508119473720660),
                        7.09100346284691107301157135337624E-3 * 6.),
                   PQP3(R3(00.045503704125649649491880526279339, 0.045503704125649649491880526279339,
                           0.454496295874350350508119473720660),
                        7.09100346284691107301157135337624E-3 * 6.),
                   PQP3(R3(00.045503704125649649491880526279339, 0.454496295874350350508119473720660,
                           0.045503704125649649491880526279339),
                        7.09100346284691107301157135337624E-3 * 6.),
                   PQP3(R3(00.454496295874350350508119473720660, 0.045503704125649649491880526279339,
                           0.045503704125649649491880526279339),
                        7.09100346284691107301157135337624E-3 * 6.)};

PQF3 const QuadratureFormular_Tet_5(5, 14, QF_TET_5);

//-------------------------------------------------------------------------------------------------
PQP3 QF_TET_7[] = {PQP3(R3(0.5, 0.5, 0.0), 0.00388007050 * 1.5),
                   PQP3(R3(0.5, 0.0, 0.5), 0.00388007050 * 1.5),
                   PQP3(R3(0.0, 0.5, 0.5), 0.00388007050 * 1.5),
                   PQP3(R3(0.0, 0.0, 0.5), 0.00388007050 * 1.5),
                   PQP3(R3(0.0, 0.5, 0.0), 0.00388007050 * 1.5),
                   PQP3(R3(0.5, 0.0, 0.0), 0.00388007050 * 1.5),

                   PQP3(R3(0.25, 0.25, 0.25), 0.07305689385 * 1.5),

                   PQP3(R3(0.078213192350, 0.078213192350, 0.078213192350), 0.042399766050 * 1.5),
                   PQP3(R3(0.078213192350, 0.078213192350, 0.765360423000), 0.042399766050 * 1.5),
                   PQP3(R3(0.078213192350, 0.765360423000, 0.078213192350), 0.042399766050 * 1.5),
                   PQP3(R3(0.765360423000, 0.078213192350, 0.078213192350), 0.042399766050 * 1.5),

                   PQP3(R3(0.121843216700, 0.121843216700, 0.121843216700), -0.250070960450 * 1.5),
                   PQP3(R3(0.121843216700, 0.121843216700, 0.634470350000), -0.250070960450 * 1.5),
                   PQP3(R3(0.121843216700, 0.634470350000, 0.121843216700), -0.250070960450 * 1.5),
                   PQP3(R3(0.634470350000, 0.121843216700, 0.121843216700), -0.250070960450 * 1.5),

                   PQP3(R3(0.332539164450, 0.332539164450, 0.332539164450), 0.01956570105 * 1.5),
                   PQP3(R3(0.332539164450, 0.332539164450, 0.002382506700), 0.01956570105 * 1.5),
                   PQP3(R3(0.332539164450, 0.002382506700, 0.332539164450), 0.01956570105 * 1.5),
                   PQP3(R3(0.002382506700, 0.332539164450, 0.332539164450), 0.01956570105 * 1.5),

                   PQP3(R3(0.1, 0.1, 0.2), 0.110229276850 * 1.5),
                   PQP3(R3(0.1, 0.2, 0.1), 0.110229276850 * 1.5),
                   PQP3(R3(0.1, 0.1, 0.6), 0.110229276850 * 1.5),
                   PQP3(R3(0.1, 0.6, 0.1), 0.110229276850 * 1.5),
                   PQP3(R3(0.1, 0.2, 0.6), 0.110229276850 * 1.5),
                   PQP3(R3(0.1, 0.6, 0.2), 0.110229276850 * 1.5),

                   PQP3(R3(0.2, 0.1, 0.1), 0.110229276850 * 1.5),
                   PQP3(R3(0.2, 0.1, 0.6), 0.110229276850 * 1.5),
                   PQP3(R3(0.2, 0.6, 0.1), 0.110229276850 * 1.5),

                   PQP3(R3(0.6, 0.1, 0.2), 0.110229276850 * 1.5),
                   PQP3(R3(0.6, 0.1, 0.1), 0.110229276850 * 1.5),
                   PQP3(R3(0.6, 0.2, 0.1), 0.110229276850 * 1.5)

};

PQF3 const QuadratureFormular_Tet_7(7, 31, QF_TET_7);

template <> const GQuadratureFormular<R3> *QF_Simplex<R3>(int exact) {

    switch (exact) {
    case 1:
        return &QuadratureFormular_Tet_1;
    case 2:
        return &QuadratureFormular_Tet_2;
    case 3:
        return &QuadratureFormular_Tet_5;
    case 4:
        return &QuadratureFormular_Tet_5;
    case 5:
        return &QuadratureFormular_Tet_5;
    default:
        return &QuadratureFormular_Tet_7;
    }
}

template <> const GQuadratureFormular<R1> *GQuadratureFormular<R1>::Default = &QF_GaussLegendre3;
template <> const GQuadratureFormular<R2> *GQuadratureFormular<R2>::Default = &QuadratureFormular_T_9;
template <> const GQuadratureFormular<R3> *GQuadratureFormular<R3>::Default = &QuadratureFormular_Tet_7;




// Sebastian's stuff on Gauss-Lobatto

