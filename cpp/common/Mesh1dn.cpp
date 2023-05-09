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
#include "RNM.hpp"
#include "libmesh5.h"
#include "Mesh1dn.hpp"

long verbosity = 2;

Mesh1::Mesh1(const char *filename) { // read the mesh

    int nt, nv, nbe;
    int ok = 0; // load(filename);
    if (ok) {
        std::ifstream f(filename);
        if (!f) {
            std::cerr << "Mesh1::Mesh1 Erreur openning " << filename << std::endl;
            exit(1);
        }
        if (verbosity)
            std::cout << " Read On file \"" << filename << "\"" << std::endl;
        f >> nv >> nt >> nbe;
        this->set(nv, nt, nbe);
        if (verbosity)
            std::cout << "  -- Nb of Vertex " << nv << " "
                      << " Nb of Seg " << nt << " , Nb of border Vertex " << nbe << std::endl;
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

Mesh1::Mesh1(int nx, R orx, R lx) {

    int mv     = nx;
    int mt     = (nx - 1);
    int mbe    = 2;
    const R hx = lx / (nx - 1);

    this->set(mv, mt, mbe);

    for (int i = 0; i < nx; i++) {
        vertices[i].X() = i * hx + orx;
    }

    for (int k = 0; k < nx - 1; k++) {
        int iv[2] = {k, k + 1};
        elements[k].set(vertices, iv, 0);
    }

    int iv[1] = {0};
    be(0).set(vertices, iv, 0);
    iv[0] = nx - 1;
    be(1).set(vertices, iv, 1);

    BuildBound();
    BuildAdj();
}
