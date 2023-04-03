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

#ifndef COMMON_MESH2DN_HPP_
#define COMMON_MESH2DN_HPP_

#include "dataStruct2D.hpp"
#include "GenericMesh.hpp"
#include <cstdlib>

class Mesh2 : public GenericMesh<Triangle2, BoundaryEdge2, Vertex2> {
  public:
    static const int D = 2;

    Mesh2() : GenericMesh<Triangle2, BoundaryEdge2, Vertex2>() {}
    Mesh2(const char *);                                 // build from mesh generated by Freefem
    Mesh2(int nx, int ny, R orx, R ory, R lx, R ly);     // build structured mesh
    void init(int nx, int ny, R orx, R ory, R lx, R ly); // build structured mesh

  private:
    Mesh2(const Mesh2 &);          // no copy constructor
    void operator=(const Mesh2 &); // no copy allowed
};

class MeshQuad2 : public GenericMesh<Quad2, BoundaryEdge2, Vertex2> {
  public:
    
    static const int D = 2;
    MeshQuad2(int nx, int ny, R orx, R ory, R lx, R ly); // build structured mesh

  private:
    MeshQuad2(const MeshQuad2 &);      // no copy constructor
    void operator=(const MeshQuad2 &); // no copy allowed
};

#endif
