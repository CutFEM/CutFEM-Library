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

#include <fstream>
#include <iostream>
#include "Mesh2dn.hpp"
#include "RNM.hpp"
#include "libmesh5.h"

Mesh2::Mesh2(const char *filename) { // read the mesh

    int nt, nv, nbe;
    int ok = 1; // load(filename);
    if (ok) {
        std::ifstream f(filename);
        if (!f) {
            std::cerr << "Mesh2::Mesh2 Erreur openning " << filename << std::endl;
            exit(1);
        }
        if (verbosity)
            std::cout << " Read On file \"" << filename << "\"" << std::endl;
        f >> nv >> nt >> nbe;
        this->set(nv, nt, nbe);
        if (verbosity)
            std::cout << "  -- Nb of Vertex " << nv << " "
                      << " Nb of Triangles " << nt << " , Nb of border edges " << nbe << std::endl;
        assert(f.good() && nt && nv);
        for (int i = 0; i < nv; i++) {
            f >> this->vertices[i];
            assert(f.good());
        }
        mes = 0;
        for (int i = 0; i < nt; i++) {
            this->t(i).Read1(f, this->vertices, nv);
            mes += t(i).measure();
        }
        mesb = 0.;
        for (int i = 0; i < nbe; i++) {
            this->be(i).Read1(f, this->vertices, nv);
            mesb += be(i).measure();
        }
    }
    BuildBound();
    BuildAdj();

    if (verbosity)
        std::cout << "   - mesh mesure = " << mes << " border mesure: " << mesb << std::endl;
}

Mesh2::Mesh2(int nx, int ny, R orx, R ory, R lx, R ly) { this->init(nx, ny, orx, ory, lx, ly); }

void Mesh2::init(int nx, int ny, R orx, R ory, R lx, R ly) {

    // int idQ[2][3] = {{0,1,2},{1,3,2}};
    int idQ2[2][3] = {{0, 1, 2}, {3, 2, 1}};
    int idQ[2][3]  = {{0, 1, 3}, {2, 0, 3}};

    int mv     = nx * ny;
    int mt     = 2 * (nx - 1) * (ny - 1);
    int mbe    = 2 * ((nx - 1) + (ny - 1));
    const R hx = lx / (nx - 1);
    const R hy = ly / (ny - 1);
    this->set(mv, mt, mbe);

    KN<int> iv(4), indT(3);

    int jt = 0;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {

            int id = 0;
            for (int jj = j; jj < j + 2; ++jj) {
                for (int ii = i; ii < i + 2; ++ii) {

                    int ivl  = ii + jj * nx; // index
                    iv(id++) = ivl;

                    vertices[ivl].x = ii * hx + orx;
                    vertices[ivl].y = jj * hy + ory;
                }
            }

            for (int l = 0; l < 2; ++l) { // create 2 elements
                for (int e = 0; e < 3; ++e) {
                    if (jt == 0 || jt == 1 || jt == mt - 1 || jt == mt - 2)
                        indT(e) = iv(idQ2[l][e]);
                    else
                        indT(e) = iv(idQ2[l][e]);
                }

                elements[jt++].set(vertices, indT, 0);
            }
        }
    }

    // create the for borders
    int lab, k = 0;
    for (int i = 0; i < nx - 1; ++i) {
        indT(0) = i;
        indT(1) = i + 1;
        lab     = 1;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indT(0) = (i + 1) * nx - 1;
        indT(1) = indT(0) + nx;
        lab     = 2;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < nx - 1; ++i) {
        indT(0) = i + nx * (ny - 1);
        indT(1) = indT(0) + 1;
        lab     = 3;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indT(0) = i * nx;
        indT(1) = indT(0) + nx;
        lab     = 4;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    assert(k == nbe);

    BuildBound();
    BuildAdj();
}

MeshQuad2::MeshQuad2(int nx, int ny, R orx, R ory, R lx, R ly) {

    // int idQ[2][3] = {{0,1,2},{1,3,2}};
    // int idQ2[2][3] = {{0,1,2},{3,2,1}};
    int indQ[4] = {0, 1, 3, 2};

    int mv     = nx * ny;
    int mt     = (nx - 1) * (ny - 1);
    int mbe    = 2 * ((nx - 1) + (ny - 1));
    const R hx = lx / (nx - 1);
    const R hy = ly / (ny - 1);
    this->set(mv, mt, mbe);

    KN<int> iv(4), indT(4);

    int jt = 0;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {

            int id = 0;
            for (int jj = j; jj < j + 2; ++jj) {
                for (int ii = i; ii < i + 2; ++ii) {

                    int ivl  = ii + jj * nx; // index
                    iv(id++) = ivl;

                    vertices[ivl].x = ii * hx + orx;
                    vertices[ivl].y = jj * hy + ory;
                }
            }
            for (int e = 0; e < 4; ++e) {
                indT(e) = iv(indQ[e]);
            }
            elements[jt++].set(vertices, indT, 0);
        }
    }

    // create the for borders
    int lab, k = 0;
    for (int i = 0; i < nx - 1; ++i) {
        indT(0) = i;
        indT(1) = i + 1;
        lab     = 1;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indT(0) = (i + 1) * nx - 1;
        indT(1) = indT(0) + nx;
        lab     = 2;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < nx - 1; ++i) {
        indT(0) = i + nx * (ny - 1);
        indT(1) = indT(0) + 1;
        lab     = 3;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indT(0) = i * nx;
        indT(1) = indT(0) + nx;
        lab     = 4;
        for (int j = 0; j < 2; ++j)
            vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
        borderelements[k++].set(vertices, indT, lab);
    }
    assert(k == nbe);

    BuildBound();
    BuildAdj();
}
