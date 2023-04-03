/*
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
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

TEST_CASE("Test BDM2 elements", "[BDM2]") {

    Mesh2 mesh(2, 2, 0., 0., 1., 1.);
    FESpace2 Vh(mesh, DataFE<Mesh2>::BDM2);
    const auto &Kref(mesh[0]);
    const auto &BDM2(DataFE<Mesh2>::BDM2);
    KNM<R> J(2, 2), invJ(2, 2);

    SECTION("Test Piola mapping") {

        //  We define a triangle A(2,3), B(4,5), C(3,5)
        // Let check that the Piolat transformation preseve the normal
        Vertex2 vertices[3];
        Triangle2 K;
        R2 A(2, 3);
        R2 B(4, 5);
        R2 C(3, 5);
        (R2 &)vertices[0] = A;
        (R2 &)vertices[1] = B;
        (R2 &)vertices[2] = C;
        int iv[3]         = {0, 1, 2};
        K.set(vertices, iv, 0);
        jacobianLinearTransformation<2>(J, K);

        REQUIRE(isEqual(J(0, 0), B.x - A.x));
        REQUIRE(isEqual(J(0, 1), C.x - A.x));
        REQUIRE(isEqual(J(1, 0), B.y - A.y));
        REQUIRE(isEqual(J(1, 1), C.y - A.y));

        auto J0i = J(0, '.');
        auto J1i = J(1, '.');
        REQUIRE(isEqual(J0i(0), B.x - A.x));
        REQUIRE(isEqual(J0i(1), C.x - A.x));
        REQUIRE(isEqual(J1i(0), B.y - A.y));
        REQUIRE(isEqual(J1i(1), C.y - A.y));

        double inv_det_J = inverseDeterminant<2>(J);
        REQUIRE(isEqual(inv_det_J,
                        1. / std::fabs(J0i(0) * J1i(1) - J0i(1) * J1i(0))));

        //   inverseJacobian<2>(J, inv_det_J, invJ);
    }

    //  SECTION("Test bf on reference element") {

    //      const auto &K(mesh[0]);
    //      KNM<R> J(2, 2), invJ(2, 2);
    //      jacobianLinearTransformation<2>(J, K);
    //      double inv_det_J = inverseDeterminant<2>(J);
    //      inverseJacobian<2>(J, inv_det_J, invJ);

    //      for (int df = 0; df < 12; ++df) {
    //          FunFEM<Mesh2> fh(Vh, BDM2.referenceBasisFunction(df));

    //          for (int i = 0; i < 12; ++i) {

    //              REQUIRE(isEqual(fh(i), (df == i)));
    //          }
    //      }
    //  }

    SECTION(
        "Test that the Piola transform preserve normal : J phi_hat . n_e = 0") {
        Vertex2 vertices[3];
        Triangle2 K;

        R2 A(2, 3);
        R2 B(4, 5);
        R2 C(3, 5);

        (R2 &)vertices[0] = A;
        (R2 &)vertices[1] = B;
        (R2 &)vertices[2] = C;

        int iv[3] = {0, 1, 2};
        K.set(vertices, iv, 0);
        jacobianLinearTransformation<2>(J, K);

        auto J0i         = J(0, '.');
        auto J1i         = J(1, '.');
        double inv_det_J = inverseDeterminant<2>(J);

        // Normal to the edges of the reference element
        std::vector<R2> normal_ref{R2(1. / std::sqrt(2), 1. / std::sqrt(2)),
                                   R2(-1, 0), R2(0, -1)};
        // Normal to the edges of the element K
        R2 normal1(-1, 0.5);
        normal1 = 1. / normal1.norm() * normal1;
        std::vector<R2> normal_K{R2(0, 1), R2(normal1.x, normal1.y),
                                 R2(1. / std::sqrt(2), -1. / std::sqrt(2))};

        // A normal basis functions has to be orthogonal to normal of the
        // edges that they don't lie on
        for (int e = 0; e < 3; ++e) {
            for (int j = 0; j <= 3; ++j) {
                R2 x = Kref.mapToReferenceElement(j * 1. / 3, e);
                for (int i = 0; i < 9; ++i) {
                    if (i >= 3 * e && i < 3 * (e + 1))
                        continue;
                    auto phi = BDM2.referenceBasisFunction(i);
                    R2 phiX(phi(x, 0), phi(x, 1));
                    REQUIRE(isEqual((phiX, normal_ref[e]), 0.));
                }
            }
        }
        //  Need to check if the mapping preserve the orthogonality
        for (int e = 0; e < 3; ++e) {
            for (int j = 0; j <= 3; ++j) {
                R2 x = Kref.mapToReferenceElement(j * 1. / 3, e);
                for (int i = 0; i < 9; ++i) {
                    if (i >= 3 * e && i < 3 * (e + 1))
                        continue;
                    auto phiHat = BDM2.referenceBasisFunction(i);
                    std::vector<double> phi_ref{phiHat(x, 0), phiHat(x, 1)};
                    KN<double> phi_K(2);
                    phi_K(0) = piolatTransformation<2>(J0i, inv_det_J, phi_ref);
                    phi_K(1) = piolatTransformation<2>(J1i, inv_det_J, phi_ref);
                    R2 phiX(phi_K(0), phi_K(1));

                    REQUIRE(isEqual((phiX, normal_K[e]), 0.));
                }
            }
        }
    }
}