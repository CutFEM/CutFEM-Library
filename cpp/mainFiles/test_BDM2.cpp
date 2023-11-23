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
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tool.hpp"

typedef TestFunction<Mesh2> FunTest;
typedef FunFEM<Mesh2> Fun_h;
typedef Mesh2 Mesh;
typedef ActiveMeshT2 CutMesh;
using fespace_t = FESpace2;
typedef CutFESpaceT2 CutSpace;

namespace Data_Darcy_Two_Unfitted {
R d_x          = 1.;
R d_y          = 1.;
R shift        = 0.5;
R inRad        = 0.2;
R interfaceRad = 0.3;  // 0.350001
R outRad       = 0.4;  // 0.4901
R pie          = M_PI; // 3.14159265359;
R fun_levelSet_in(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - inRad;
}
R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}
R fun_levelSet_out(const R2 P, const int i) {
    return outRad - sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift));
}
R fun_test(const R2 P, const int i, int d) { return d; } // [Eriks Scotti example]
R rad2 = interfaceRad * interfaceRad;
R mu_G = 2 * interfaceRad / (4 * cos(rad2) + 3); // xi0*mu_G =
// 1/8*2/3*1/4

} // namespace Data_Darcy_Two_Unfitted
using namespace Data_Darcy_Two_Unfitted;
int main(int argc, char **argv) {
    typedef TestFunction<Mesh2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    int nx = 21; // 6
    int ny = 21; // 6

    Mesh Kh(nx, ny, 0., 0., d_x + globalVariable::Epsilon, d_y + globalVariable::Epsilon);
    const R h_i  = 1. / (nx - 1);
    const R invh = 1. / h_i;
    Space Lh(Kh, DataFE<Mesh2>::P1);

    Fun_h levelSet_out(Lh, fun_levelSet_out);
    InterfaceLevelSet<Mesh> boundary_out(Kh, levelSet_out);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);
    Fun_h levelSet_in(Lh, fun_levelSet_in);
    InterfaceLevelSet<Mesh> boundary_in(Kh, levelSet_in);

    // Cut mesh
    ActiveMesh<Mesh> Kh_i(Kh);
    Kh_i.truncate(boundary_in, -1);
    Kh_i.truncate(boundary_out, -1);
    // Kh_i.add(interface, -1);

    Paraview<Mesh> writer(Kh_i, "Kh_i.vtk");

    std::cout << " hey " << std::endl;
}

// int main(int argc, char **argv) {

//     // MPIcf cfMPI(argc, argv);

//     int nx = 11; // 6
//     int ny = 11; // 6

//     std::vector<double> uPrint, pPrint, divPrint, divPrintLoc, maxDivPrint, h, convuPr, convpPr, convdivPr,
//         convdivPrLoc, convmaxdivPr;

//     int iters = 4;
//     for (int i = 0; i < iters; ++i) {
//         Mesh Kh(nx, ny, 0., 0., 1, 1.);

//         fespace_t Vh1(Kh, DataFE<Mesh>::BDM2);
//         fespace_t Vh2(Kh, DataFE<Mesh>::RT1);
//         fespace_t Vh0(Kh, DataFE<Mesh>::RT1);

//         double h_i = 1. / (nx - 1);
//         Rn uh_data(Vh1.get_nb_dof(), 0.);
//         Fun_h uh(Vh1, uh_data);
//         Fun_h fv(Vh2, fun_exact_u);

//         Fun_h f0(Vh0, fun_exact_u);
//         projection(fv, uh);

//         R errU = L2norm(uh, fun_exact_u, 0, 2);
//         // [PLOTTING]
//         // {
//         //     Paraview<Mesh> writer(Kh, "BDM2_projection_" + std::to_string(i)
//         //     +
//         //                                   ".vtk");
//         //     writer.add(uh, "velocity", 0, 2);
//         //     writer.add(fv, "velocityBDM2", 0, 2);
//         //     writer.add(f0, "velocityRT1", 0, 2);
//         // }

//         uPrint.push_back(errU);
//         h.push_back(h_i);

//         if (i == 0) {
//             convuPr.push_back(0);
//         } else {
//             convuPr.push_back(log(uPrint[i] / uPrint[i - 1]) / log(h[i] / h[i - 1]));
//         }

//         nx = 2 * nx - 1;
//         ny = 2 * ny - 1;
//     }
//     std::cout << "\n"
//               << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setfill(' ') << "err u" <<
//               std::setw(15)
//               << std::setfill(' ') << "conv u" << std::setw(15) << "\n"
//               << std::endl;
//     for (int i = 0; i < uPrint.size(); ++i) {
//         std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill('
//         ')
//                   << uPrint[i] << std::setw(15) << std::setfill(' ') << convuPr[i] << std::endl;
//     }
// }
