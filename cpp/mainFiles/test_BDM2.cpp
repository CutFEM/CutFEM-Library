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

typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;
typedef Mesh2 Mesh;
typedef ActiveMeshT2 CutMesh;
typedef FESpace2 Space;
typedef CutFESpaceT2 CutSpace;

R fun_exact_u(double *P, int compInd) {
    return (compInd == 0) ? -P[0] : P[1] - 1;
}

int main(int argc, char **argv) {

    // MPIcf cfMPI(argc, argv);

    int nx = 11; // 6
    int ny = 11; // 6

    std::vector<double> uPrint, pPrint, divPrint, divPrintLoc, maxDivPrint, h,
        convuPr, convpPr, convdivPr, convdivPrLoc, convmaxdivPr;

    int iters = 4;
    for (int i = 0; i < iters; ++i) {
        Mesh Kh(nx, ny, 0., 0., 1, 1.);

        Space Vh1(Kh, DataFE<Mesh>::BDM2);
        Space Vh2(Kh, DataFE<Mesh>::RT1);
        Space Vh0(Kh, DataFE<Mesh>::RT1);

        double h_i = 1. / (nx - 1);
        Rn uh_data(Vh1.get_nb_dof(), 0.);
        Fun_h uh(Vh1, uh_data);
        Fun_h fv(Vh2, fun_exact_u);

        Fun_h f0(Vh0, fun_exact_u);
        projection(fv, uh);

        R errU = L2norm(uh, fun_exact_u, 0, 2);
        // [PLOTTING]
        // {
        //     Paraview<Mesh> writer(Kh, "BDM2_projection_" + std::to_string(i)
        //     +
        //                                   ".vtk");
        //     writer.add(uh, "velocity", 0, 2);
        //     writer.add(fv, "velocityBDM2", 0, 2);
        //     writer.add(f0, "velocityRT1", 0, 2);
        // }

        uPrint.push_back(errU);
        h.push_back(h_i);

        if (i == 0) {
            convuPr.push_back(0);
        } else {
            convuPr.push_back(log(uPrint[i] / uPrint[i - 1]) /
                              log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h"
              << std::setfill(' ') << "err u" << std::setw(15)
              << std::setfill(' ') << "conv u" << std::setw(15) << "\n"
              << std::endl;
    for (int i = 0; i < uPrint.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i]
                  << std::setw(15) << std::setfill(' ') << uPrint[i]
                  << std::setw(15) << std::setfill(' ') << convuPr[i]
                  << std::endl;
    }
}
