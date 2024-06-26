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
#include <cstring>
#include <algorithm>
#include "libmesh5.h"
#include "RNM.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "base_interface.hpp"

// Mesh3::Mesh3(const std::string filename) {
//     std::ifstream f(filename.c_str());
//     if (filename.rfind(".msh") == filename.length() - 4)
//         readmsh(f);
//     else {
//         std::cout << " not a good format" << std::endl;
//     }

//     BuildBound();
//     if (nt > 0) {
//         BuildAdj();
//     }

//     verbosity = 3;
//     if (verbosity > 2)
//         std::cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << std::endl;
//     if (verbosity)
//         std::cout << "  -- Mesh3 : " << filename << ", d " << 3 << ", n Tet " << nt << ", n Vtx " << nv << " n Bord "
//                   << nbe << std::endl;
//     ffassert(mes >= 0); // add F. Hecht sep 2009.
// }

Mesh3::Mesh3(const std::string filename, MeshFormat type_mesh) { // read the mesh

    std::ifstream f(filename);
    if (!f) {
        std::cerr << "Mesh3::Mesh3 Erreur openning " << filename << std::endl;
        exit(1);
    }
    if (verbosity)
        std::cout << " Read On file \"" << filename << "\"" << std::endl;

    // if (filename.rfind(".msh") == filename.length() - 4) {
    //     if (verbosity)
    //         std::cout << "  -- Read msh file " << std::endl;
    //     readMsh(f);
    //     // return;
    // } else if (filename.rfind(".mesh") == filename.length() - 5) {
    //     if (verbosity)
    //         std::cout << "  -- Read mesh file " << std::endl;
    //     readMesh(f);
    // }
    if (type_mesh == MeshFormat::mesh_gmsh)
        readMeshGmsh(f);
    else if (type_mesh == MeshFormat::mesh_freefem) {
        //readMeshFreefem(f);
        std::cerr << "No freefem implementation " << filename << std::endl;
        exit(1);
    }
    else {
        std::cerr << "No mesh format specified " << filename << std::endl;
        exit(1);
    }

    BuildBound();
    BuildAdj();

    // if (verbosity)
    std::cout << "   - mesh mesure = " << mes << " border mesure: " << mesb << std::endl;
}

void Mesh3::readmsh(std::ifstream &f) {
    f >> nv >> nt >> nbe;
    if (verbosity > 2)
        std::cout << " GRead : nv " << nv << " " << nt << " " << nbe << std::endl;
    this->vertices       = new Vertex[nv];
    this->elements       = new Element[nt];
    this->borderelements = new BorderElement[nbe];
    for (int k = 0; k < nv; k++) {
        Vertex &P = this->vertices[k];
        f >> P.x >> P.y >> P.z >> P.lab;
    }
    mes  = 0.;
    mesb = 0.;

    if (nt != 0) {
        for (int k = 0; k < nt; k++) {
            int i[4], lab;
            Element &K(this->elements[k]);
            f >> i[0] >> i[1] >> i[2] >> i[3] >> lab;
            K.set(this->vertices, i, lab);
            mes += K.measure();
        }
    }
    for (int k = 0; k < nbe; k++) {
        int i[4], lab;
        BorderElement &K(this->borderelements[k]);
        f >> i[0] >> i[1] >> i[2] >> lab;
        K.set(this->vertices, i, lab);
        mesb += K.measure();
    }
}

void Mesh3::readMeshGmsh(std::ifstream &f) {

    std::string field;
    while (std::getline(f, field)) {
        if (field.find("Vertices") != std::string::npos) {
            f >> nv;
            std::cout << "  -- Nb of Vertex " << nv << std::endl;
            assert(nv);
            vertices = new Vertex3[nv];

            R3 P;
            for (int i = 0; i < nv; i++) {
                f >> P >> vertices[i].lab;
                (R3 &)vertices[i] = P;
            }
        }

        if (field.find("Triangles") != std::string::npos) {
            f >> nbe;
            std::cout << "  -- Nb of border triangles " << nbe << std::endl;
            assert(nbe);
            borderelements = new Triangle3[nbe];
            mesb           = 0.;
            for (int i = 0; i < nbe; i++) {
                this->be(i).Read1(f, this->vertices, nv);
                mesb += be(i).measure();
            }
        }

        if (field.find("Tetrahedra") != std::string::npos) {
            f >> nt;
            std::cout << "  -- Nb of Tetrahedra " << nt << std::endl;
            assert(nt);
            elements = new Tet[nt];
            mes      = 0;
            for (int i = 0; i < nt; i++) {
                this->t(i).Read1(f, this->vertices, nv);
                mes += t(i).measure();
            }
        }
    }
}


Mesh3::Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb) {

    nv  = nnv;
    nt  = nnt;
    nbe = nnbe;

    vertices       = vv;
    elements       = tt;
    borderelements = bb;

    mes  = 0.;
    mesb = 0.;

    for (int i = 0; i < nt; i++)
        mes += this->elements[i].measure();

    for (int i = 0; i < nbe; i++)
        mesb += this->be(i).measure();

    //  Add FH to be consitant we all constructor ...  July 09
    BuildBound();
    if (nt > 0) {
        BuildAdj();
    }
    //  end add
    // if(verbosity>1)
    // 	  std::cout << "  -- End of read: mesure = " << mes << " border mesure "
    // << mesb << std::endl;

    //    if(verbosity>1)
    //      std::cout << "  -- End of read: mesure = " << mes << " border mesure
    //      " << mesb << std::endl;
    assert(mes >= 0.);
}

Mesh3::Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb) {

    nv  = nnv;
    nbe = nnbe;

    vertices       = vv;
    borderelements = bb;

    mes  = 0.;
    mesb = 0.;

    for (int i = 0; i < nbe; i++)
        mesb += this->be(i).measure();

    //  Add FH to be consitant we all constructor ...  July 09
    BuildBound();
    if (nt > 0) {
        BuildAdj();
    }
    //  end add

    if (verbosity > 1)
        std::cout << "  -- End of Construct  mesh3: mesure = " << mes << " border mesure " << mesb << std::endl;
    ffassert(mes >= 0); // add F. Hecht sep 2009.
}

int signe_permutation(int i0, int i1, int i2, int i3) {
    int p = 1;
    if (i0 > i1)
        std::swap(i0, i1), p = -p;
    if (i0 > i2)
        std::swap(i0, i2), p = -p;
    if (i0 > i3)
        std::swap(i0, i3), p = -p;
    if (i1 > i2)
        std::swap(i1, i2), p = -p;
    if (i1 > i3)
        std::swap(i1, i3), p = -p;
    if (i2 > i3)
        std::swap(i2, i3), p = -p;
    return p;
}

Mesh3::Mesh3(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz) {

    const int idQ[6][4] = {{1, 3, 2, 5}, {1, 2, 4, 5}, {0, 1, 2, 4}, {6, 3, 5, 2}, {6, 2, 5, 4}, {7, 3, 5, 6}};
    int mv              = nx * ny * nz;
    int mt              = 6 * (nx - 1) * (ny - 1) * (nz - 1);
    int mbe             = 4 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) + (nz - 1) * (ny - 1));
    const R hx          = lx / (nx - 1);
    const R hy          = ly / (ny - 1);
    const R hz          = lz / (nz - 1);

    this->set(mv, mt, mbe);

    KN<int> iv(8), indT(4);

    int jt = 0;
    for (int k = 0; k < nz - 1; k++) {
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {

                int id = 0;
                for (int kk = k; kk < k + 2; ++kk) {
                    for (int jj = j; jj < j + 2; ++jj) {
                        for (int ii = i; ii < i + 2; ++ii) {

                            int ivl  = ii + jj * nx + kk * nx * ny; // index
                            iv(id++) = ivl;

                            vertices[ivl].x = ii * hx + orx;
                            vertices[ivl].y = jj * hy + ory;
                            vertices[ivl].z = kk * hz + orz;
                        }
                    }
                }

                for (int l = 0; l < 6; ++l) { // create 2 elements
                    for (int e = 0; e < 4; ++e) {
                        indT(e) = iv(idQ[l][e]);
                    }
                    elements[jt++].set(vertices, indT, 0);
                }
            }
        }
    }

    assert(jt == nt);

    const R xmin = orx, ymin = ory, zmin = orz;
    const R xmax = orx + lx, ymax = ory + ly, zmax = orz + lz;

    // create the for borders
    int kf = 0;
    int indV[3], indKV[3];
    const R step = 1e-12;
    for (int jt = 0; jt < nt; ++jt) { // loop over element
        const Element &K((*this)[jt]);
        for (int k = 0; k < 4; ++k) { // loop over faces
            for (int i = 0; i < 3; ++i) {
                indV[i]  = Element::nvface[k][i];
                indKV[i] = (*this)(K[indV[i]]);
            }

            const int iv1 = Element::nvface[k][0];
            const int iv2 = Element::nvface[k][1];
            const int iv3 = Element::nvface[k][2];

            //  Faces on the plane x=orx
            //-----------------------------------------------------
            if ((fabs(K[iv1].x - xmin) < step) && (fabs(K[iv2].x - xmin) < step) && (fabs(K[iv3].x - xmin) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 1);
                borderelements[kf++].set(vertices, indKV, 1);
            }

            //  Faces on the plane x=orx + glx
            //-----------------------------------------------------
            else if ((fabs(K[iv1].x - xmax) < step) && (fabs(K[iv2].x - xmax) < step) &&
                     (fabs(K[iv3].x - xmax) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 2);
                borderelements[kf++].set(vertices, indKV, 2);
            }

            //  Faces on the plane y=ory
            //-----------------------------------------------------
            else if ((fabs(K[iv1].y - ymin) < step) && (fabs(K[iv2].y - ymin) < step) &&
                     (fabs(K[iv3].y - ymin) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 3);
                borderelements[kf++].set(vertices, indKV, 3);
            }

            //  Faces on the plane y=ory + gly
            //-----------------------------------------------------
            else if ((fabs(K[iv1].y - ymax) < step) && (fabs(K[iv2].y - ymax) < step) &&
                     (fabs(K[iv3].y - ymax) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 4);
                borderelements[kf++].set(vertices, indKV, 4);
            }
            //  Faces on the plane z=orz
            //-----------------------------------------------------
            else if ((fabs(K[iv1].z - zmin) < step) && (fabs(K[iv2].z - zmin) < step) &&
                     (fabs(K[iv3].z - zmin) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 5);
                borderelements[kf++].set(vertices, indKV, 5);
            }

            //  Faces on the plane z=orz + glz
            //-----------------------------------------------------
            else if ((fabs(K[iv1].z - zmax) < step) && (fabs(K[iv2].z - zmax) < step) &&
                     (fabs(K[iv3].z - zmax) < step)) {

                for (int j = 0; j < 3; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 6);
                borderelements[kf++].set(vertices, indKV, 6);
            }
        } // end loop over faces
    }

    assert(kf == nbe);

    BuildBound();
    if (nt > 0) {
        BuildAdj();
    }

    // std::cout << "  -- End of read: mesure = " << mes << " border mesure " <<
    // mesb << std::endl; std::cout << "  -- Mesh3 : "<< " d "<< 3  << ", n Tet "
    // << nt << ", n Vtx "
    // 	    << nv << " n Bord " << nbe << std::endl;
    ffassert(mes >= 0);
}

MeshHexa::MeshHexa(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz) {

    const int idQ[8] = {0, 1, 3, 2, 4, 5, 7, 6};
    int mv           = nx * ny * nz;
    int mt           = (nx - 1) * (ny - 1) * (nz - 1);
    int mbe          = 2 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) + (nz - 1) * (ny - 1));
    const R hx       = lx / (nx - 1);
    const R hy       = ly / (ny - 1);
    const R hz       = lz / (nz - 1);

    this->set(mv, mt, mbe);

    KN<int> iv(8), indT(8);

    int jt = 0;
    for (int k = 0; k < nz - 1; k++) {
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {

                int id = 0;
                for (int kk = k; kk < k + 2; ++kk) {
                    for (int jj = j; jj < j + 2; ++jj) {
                        for (int ii = i; ii < i + 2; ++ii) {

                            int ivl  = ii + jj * nx + kk * nx * ny; // index
                            iv(id++) = ivl;

                            vertices[ivl].x = ii * hx + orx;
                            vertices[ivl].y = jj * hy + ory;
                            vertices[ivl].z = kk * hz + orz;
                        }
                    }
                }
                for (int e = 0; e < 8; ++e) {
                    indT(e) = iv(idQ[e]);
                }
                elements[jt++].set(vertices, indT, 0);
            }
        }
    }

    assert(jt == nt);

    const R xmin = orx, ymin = ory, zmin = orz;
    const R xmax = orx + lx, ymax = ory + ly, zmax = orz + lz;

    // create the for borders
    int kf = 0;
    int indV[4], indKV[4];
    const R step = 1e-12;
    for (int jt = 0; jt < nt; ++jt) { // loop over element
        const Element &K((*this)[jt]);
        for (int k = 0; k < 6; ++k) { // loop over faces
            for (int i = 0; i < 4; ++i) {
                indV[i]  = Element::nvface[k][i];
                indKV[i] = (*this)(K[indV[i]]);
            }

            const int iv1 = Element::nvface[k][0];
            const int iv2 = Element::nvface[k][1];
            const int iv3 = Element::nvface[k][2];
            const int iv4 = Element::nvface[k][3];

            //  Faces on the plane x=orx
            //-----------------------------------------------------
            if ((fabs(K[iv1].x - xmin) < step) && (fabs(K[iv2].x - xmin) < step) && (fabs(K[iv3].x - xmin) < step) &&
                (fabs(K[iv4].x - xmin) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 1);
                borderelements[kf++].set(vertices, indKV, 1);
            }

            //  Faces on the plane x=orx + glx
            //-----------------------------------------------------
            else if ((fabs(K[iv1].x - xmax) < step) && (fabs(K[iv2].x - xmax) < step) &&
                     (fabs(K[iv3].x - xmax) < step) && (fabs(K[iv4].x - xmax) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 2);
                borderelements[kf++].set(vertices, indKV, 2);
            }

            //  Faces on the plane y=ory
            //-----------------------------------------------------
            else if ((fabs(K[iv1].y - ymin) < step) && (fabs(K[iv2].y - ymin) < step) &&
                     (fabs(K[iv3].y - ymin) < step) && (fabs(K[iv4].y - ymin) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 3);
                borderelements[kf++].set(vertices, indKV, 3);
            }

            //  Faces on the plane y=ory + gly
            //-----------------------------------------------------
            else if ((fabs(K[iv1].y - ymax) < step) && (fabs(K[iv2].y - ymax) < step) &&
                     (fabs(K[iv3].y - ymax) < step) && (fabs(K[iv4].y - ymax) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 4);
                borderelements[kf++].set(vertices, indKV, 4);
            }
            //  Faces on the plane z=orz
            //-----------------------------------------------------
            else if ((fabs(K[iv1].z - zmin) < step) && (fabs(K[iv2].z - zmin) < step) &&
                     (fabs(K[iv3].z - zmin) < step) && (fabs(K[iv4].z - zmin) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 5);
                borderelements[kf++].set(vertices, indKV, 5);
            }

            //  Faces on the plane z=orz + glz
            //-----------------------------------------------------
            else if ((fabs(K[iv1].z - zmax) < step) && (fabs(K[iv2].z - zmax) < step) &&
                     (fabs(K[iv3].z - zmax) < step) && (fabs(K[iv4].z - zmax) < step)) {

                for (int j = 0; j < 4; ++j)
                    vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 6);
                borderelements[kf++].set(vertices, indKV, 6);
            }
        } // end loop over faces
    }

    assert(kf == nbe);

    BuildBound();
    if (nt > 0) {
        BuildAdj();
    }

    ffassert(mes >= 0);
}
