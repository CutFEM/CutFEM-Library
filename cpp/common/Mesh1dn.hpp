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
#ifndef COMMON_MESH1DN_HPP_
#define COMMON_MESH1DN_HPP_

#include "GenericMesh.hpp"
#include <cstdlib>

typedef GenericVertex<R1> Vertex1;

// class SignPatternTrait1;
class RefPatch2;
// class RefPartition1;
// class Partition1;

class Mesh1 : public GenericMesh<Seg1, BoundaryPoint1, Vertex1> {
 public:
   // typedef SignPatternTrait1 SignPattern;
   typedef RefPatch2 RefPatch;
   // typedef RefPartition1 RefPartition;
   // typedef Partition1 Partition;

   Mesh1(int nx, R orx, R lx); // build structured mesh
   Mesh1(const char *);        //
   const Element *Find(R1 P, R1 &Phat, bool &outside,
                       const Element *tstart) const;

 private:
   int load(const std::string &filename);
   Mesh1(const Mesh1 &);          // pas de construction par copie
   void operator=(const Mesh1 &); // pas affectation par copy
};

#endif
