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

Mesh2::Mesh2(const std::string filename, MeshFormat type_mesh) { // read the mesh

    std::ifstream f(filename);
    if (!f) {
        std::cerr << "Mesh2::Mesh2 Erreur openning " << filename << std::endl;
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
    else if (type_mesh == MeshFormat::mesh_freefem)
        readMeshFreefem(f);
    else {
        std::cerr << "Mesh2::Mesh2 Erreur openning " << filename << std::endl;
        exit(1);
    }

    BuildBound();
    BuildAdj();

    // if (verbosity)
    std::cout << "   - mesh mesure = " << mes << " border mesure: " << mesb << std::endl;
}

void Mesh2::readMeshGmsh(std::ifstream &f) {

    std::string field;
    while (std::getline(f, field)) {
        if (field.find("Vertices") != std::string::npos) {
            f >> nv;
            std::cout << "  -- Nb of Vertex " << nv << std::endl;
            assert(nv);
            vertices = new Vertex2[nv];
            double Pz;
            R2 P;
            for (int i = 0; i < nv; i++) {
                f >> P >> Pz >> vertices[i].lab;
                (R2 &)vertices[i] = P;
            }
        }

        if (field.find("Edges") != std::string::npos) {
            f >> nbe;
            std::cout << "  -- Nb of border edges " << nbe << std::endl;
            assert(nbe);
            borderelements = new BoundaryEdge2[nbe];
            mesb           = 0.;
            for (int i = 0; i < nbe; i++) {
                this->be(i).Read1(f, this->vertices, nv);
                mesb += be(i).measure();
            }
        }

        if (field.find("Triangles") != std::string::npos) {
            f >> nt;
            std::cout << "  -- Nb of Triangles " << nt << std::endl;
            assert(nt);
            elements = new Triangle2[nt];
            mes      = 0;
            for (int i = 0; i < nt; i++) {
                this->t(i).Read1(f, this->vertices, nv);
                mes += t(i).measure();
            }
        }
    }
}

void Mesh2::readMeshFreefem(std::ifstream &f) {
    int mv, mt, mbe;
    f >> mv >> mt >> mbe;
    if (verbosity)
        std::cout << "  -- Nb of Vertex " << mv << " " << " Nb of Triangles " << mt << " , Nb of border edges " << mbe
                  << std::endl;
    this->set(mv, mt, mbe);

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

Mesh2 refine(const Mesh2 &Th) {
    using Element       = typename Mesh2::Element;
    using BorderElement = typename Mesh2::BorderElement;
    using Key           = std::pair<size_t, size_t>;
    int idx_ref[4][3]   = {{0, 5, 4}, {1, 3, 5}, {2, 4, 3}, {5, 3, 4}};

    Mesh2 Th2;

    // Euler characteristic V-E+F = 1
    size_t E = Th.nv + Th.nt - 1;

    int mt  = 4 * Th.nt;
    int mv  = Th.nv + E;
    int mbe = Th.nbe * 2;
    Th2.set(mv, mt, mbe);

    std::cout << "  -- Nb of Vertex " << mv << " " << " Nb of Triangles " << mt << " , Nb of border edges " << mbe
              << std::endl;

    size_t tinfty = -1;
    std::map<Key, size_t> new_nodes;

    //   int kt = 0;
    size_t iv = 0;
    for (int k = 0; k < Th.nt; ++k) {
        const Element &K(Th[k]); // the element
        std::array<int, 6> idxK;
        int ivl = 0;
        // loop over the nodes
        for (int i = 0; i < 3; ++i) {
            int iv_glob = Th(K[i]);
            Key ki(iv_glob, tinfty);

            auto pk = new_nodes.find(ki);
            if (pk == new_nodes.end()) {
                idxK[ivl++]        = iv;
                new_nodes[ki]      = iv;
                R2 P               = K[i];
                Th2.vertices[iv].x = P.x;
                Th2.vertices[iv].y = P.y;
                iv++;
            } else {
                idxK[ivl++] = pk->second;
            }
        }

        for (int i = 0; i < Element::ne; ++i) {
            int i0  = Element::nvedge[i][0];
            int i1  = Element::nvedge[i][1];
            Key ki  = (Th(K[i0]) < Th(K[i1])) ? Key(Th(K[i0]), Th(K[i1])) : Key(Th(K[i1]), Th(K[i0]));
            auto pk = new_nodes.find(ki);
            if (pk == new_nodes.end()) {
                idxK[ivl++]        = iv;
                new_nodes[ki]      = iv;
                R2 P               = 0.5 * (K[i0] + K[i1]);
                Th2.vertices[iv].x = P.x;
                Th2.vertices[iv].y = P.y;
                iv++;
            } else {
                idxK[ivl++] = pk->second;
            }
        }

        for (int i = 0; i < 4; ++i) {
            std::array<int, 3> indT = {idxK[idx_ref[i][0]], idxK[idx_ref[i][1]], idxK[idx_ref[i][2]]};
            Th2.elements[4 * k + i].set(Th2.vertices, indT.data(), 0);
        }
    }

    int n = 0;
    for (int ke = 0; ke < Th.nbe; ++ke) { // loop over boundary element
        int kf;
        int k = Th.BoundaryElement(ke, kf);
        const Element &K(Th[k]);
        const BorderElement &Be(Th.be(ke));

        std::array<int, 3> idxBE;

        Key k0(Th(K[Element::nvedge[kf][0]]), tinfty);
        auto pk = new_nodes.find(k0);
        if (pk == new_nodes.end()) {
            std::cout << "Error: node not found" << std::endl;
            exit(1);
        }
        idxBE[0] = pk->second;

        Key k1(Th(K[Element::nvedge[kf][1]]), tinfty);
        pk = new_nodes.find(k1);
        if (pk == new_nodes.end()) {
            std::cout << "Error: node not found" << std::endl;
            exit(1);
        }
        idxBE[2] = pk->second;

        int iv0 = Th(Be[0]);
        int iv1 = Th(Be[1]);
        Key ki  = (iv0 < iv1) ? Key(iv0, iv1) : Key(iv1, iv0);
        pk      = new_nodes.find(ki);
        if (pk == new_nodes.end()) {
            std::cout << "Error: node not found" << std::endl;
            exit(1);
        }
        idxBE[1] = pk->second;

        Th2.borderelements[n++].set(Th2.vertices, idxBE.data(), Be.lab);
        Th2.borderelements[n++].set(Th2.vertices, idxBE.data() + 1, Be.lab);
    }

    Th2.BuildBound();
    Th2.BuildAdj();

    return Th2;
}

Mesh2 refine_barycentric(const Mesh2 &Th) {
    using Element       = typename Mesh2::Element;
    using BorderElement = typename Mesh2::BorderElement;
    using Key           = std::pair<size_t, size_t>;
    int idx_ref[3][3]   = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}};

    Mesh2 Th2;
    size_t E = Th.nt;
    int mt   = 3 * Th.nt;
    int mv   = Th.nv + E;
    int mbe  = Th.nbe;
    Th2.set(mv, mt, mbe);

    std::cout << "  -- Nb of Vertex " << mv << " " << " Nb of Triangles " << mt << " , Nb of border edges " << mbe
              << std::endl;

    size_t tinfty = -1;
    std::map<Key, size_t> new_nodes;

    //   int kt = 0;
    size_t iv  = 0;
    size_t if0 = Th.nv + 1;
    for (int k = 0; k < Th.nt; ++k) {
        const Element &K(Th[k]); // the element
        std::array<int, 6> idxK;
        int ivl = 0;
        int iv_glob;
        // loop over the nodes
        for (int i = 0; i < 3; ++i) {
            iv_glob = Th(K[i]);
            Key ki(iv_glob, tinfty);

            auto pk = new_nodes.find(ki);
            if (pk == new_nodes.end()) {
                idxK[ivl++]        = iv;
                new_nodes[ki]      = iv;
                R2 P               = K[i];
                Th2.vertices[iv].x = P.x;
                Th2.vertices[iv].y = P.y;
                iv++;
            } else {
                idxK[ivl++] = pk->second;
            }
        }

        iv_glob = k + if0;
        Key ki(iv_glob, tinfty);
        auto pk = new_nodes.find(ki);
        if (pk == new_nodes.end()) {
            idxK[ivl++]        = iv;
            new_nodes[ki]      = iv;
            R2 P               = 1. / 3 * (K[0] + K[1] + K[2]);
            Th2.vertices[iv].x = P.x;
            Th2.vertices[iv].y = P.y;
            iv++;
        } else {
            idxK[ivl++] = pk->second;
        }

        for (int i = 0; i < 3; ++i) {
            std::array<int, 3> indT = {idxK[idx_ref[i][0]], idxK[idx_ref[i][1]], idxK[idx_ref[i][2]]};
            Th2.elements[3 * k + i].set(Th2.vertices, indT.data(), 0);
        }
    }

    int n = 0;
    for (int ke = 0; ke < Th.nbe; ++ke) { // loop over boundary element
        int kf;
        int k = Th.BoundaryElement(ke, kf);
        const Element &K(Th[k]);
        const BorderElement &Be(Th.be(ke));

        std::array<int, 2> idxBE;

        Key k0(Th(K[Element::nvedge[kf][0]]), tinfty);
        auto pk = new_nodes.find(k0);
        if (pk == new_nodes.end()) {
            std::cout << "Error: node not found" << std::endl;
            exit(1);
        }
        idxBE[0] = pk->second;
        Key k1(Th(K[Element::nvedge[kf][1]]), tinfty);
        pk = new_nodes.find(k1);
        if (pk == new_nodes.end()) {
            std::cout << "Error: node not found" << std::endl;
            exit(1);
        }
        idxBE[1] = pk->second;

        Th2.borderelements[n++].set(Th2.vertices, idxBE.data(), Be.lab);
    }

    Th2.BuildBound();
    Th2.BuildAdj();

    return Th2;
}
