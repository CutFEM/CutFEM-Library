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

#ifndef COMMON_MESH3DN_HPP_
#define COMMON_MESH3DN_HPP_

#include <cstdlib>
#include "dataStruct3D.hpp"
#include "GenericMesh.hpp"

class Mesh3 : public GenericMesh<Tet, Triangle3, Vertex3> {

  public:
    static const int D = 3;

    Mesh3() {}
    //Mesh3(const std::string);
    Mesh3(const std::string filename, MeshFormat type_mesh);
    Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb);
    Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb); // surface mesh
    Mesh3(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz);

  private:
    void readmsh(std::ifstream &f);
    void readMeshGmsh(std::ifstream &f);

    Mesh3(const Mesh3 &);          // pas de construction par copie
    void operator=(const Mesh3 &); // pas affectation par copy
};

class MeshHexa : public GenericMesh<Hexa, Quad3, Vertex3> {
  public:
    static const int D = 3;

    MeshHexa(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly,
             R lz); // build structured mesh
  private:
    MeshHexa(const MeshHexa &);       // no copy constructor
    void operator=(const MeshHexa &); // no copy allowed
};

template <typename Mesh> void setLabelBorder(Mesh &Th, R3 normal, int newLabel) {
    typedef typename Mesh::Rd Rd;
    for (int ifac = 0; ifac < Th.nbe; ifac += 1) {
        // for( int ifac = Th.first_element(); ifac < Th.last_boundary_element();
        // ifac+=Th.next_element()) {
        typename Mesh::BorderElement &face(Th.be(ifac));
        int ifaceK; // index of face of triangle corresp to edge (0,1,2)
        const int k = Th.BoundaryElement(ifac,
                                         ifaceK); // index of element (triangle), ifaceK gets modified inside

        int ib = ifaceK;
        if (Th.ElementAdj(k, ib) != -1)
            continue; // not on the boundary. Because mesh built with buildlayers
                      // in freefem++

        const typename Mesh::Element &T(Th[k]);
        Rd normal_ext = T.N(ifaceK);

        Rd diff_normal(normal - normal_ext);
        if (diff_normal.norm() < 1e-15) {
            face.lab = newLabel;
        }
    }
}

typedef Mesh3 MeshT3;
typedef MeshHexa MeshQ3;

#endif
