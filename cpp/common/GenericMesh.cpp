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

#include "GenericMesh.hpp"

template <typename T, typename B, typename V> void GenericMesh<T, B, V>::BuildAdj() {
    if (TheAdjacencesLink != 0)
        return; //  already build

    BuildAdjacencyOfMesh<GenericMesh<T, B, V>> a(*this);
}

template <typename GMesh> BuildAdjacencyOfMesh<GMesh>::BuildAdjacencyOfMesh(GMesh &m) : mesh(m) {
    initializeArray();

    findAdjacencyElement();

    findBoundaryElement();
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::initializeArray() {
    mesh.TheAdjacencesLink       = new int[mesh.nea * mesh.nt];
    mesh.BoundaryElementHeadLink = new int[mesh.nbe];
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::findAdjacencyElement() {
    nk = 0, nba = 0;
    ne = 0;
    for (int k = 0; k < mesh.nt; ++k) {
        for (int i = 0; i < mesh.nea; ++i) {
            addFace(k, i);
        }
    }
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::addFace(const int k, const int i) {
    SArray a(mesh.itemadj(k, i)); // sort nodes on the adj item
    auto p = h.find(a);
    if (p == h.end()) {
        ne++;
        addFaceFirstTime(a);
    } else {
        addFaceAlreadySeen(p, a);
    }
    ++nk;
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::addFaceFirstTime(const SArray &a) {
    h[a]                       = nk;
    mesh.TheAdjacencesLink[nk] = -1;
    nba++;
}

template <typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::addFaceAlreadySeen(typename std::map<SArray, int>::iterator p, const SArray &a) {
    assert(p->second >= 0);
    mesh.TheAdjacencesLink[nk]        = p->second;
    mesh.TheAdjacencesLink[p->second] = nk;

    p->second = -1 - nk;
    nba--;
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::findBoundaryElement() {
    for (int k = 0; k < mesh.nbe; ++k) {
        addBoundary(k);
    }
}

template <typename GMesh> void BuildAdjacencyOfMesh<GMesh>::addBoundary(const int k) {
    int err = 0;
    SArray a(mesh.itembe(k));
    auto p = h.find(a);
    if (p == h.end()) {
        err++;
        if (err == 1)
            std::cerr << "Err  Border element not in mesh \n";
        if (err < 10)
            std::cerr << " \t " << k << " " << a << std::endl;
    } else {
        mesh.BoundaryElementHeadLink[k] = p->second < 0 ? -p->second - 1 : p->second;
    }
}

template <typename T, typename B, typename V> void GenericMesh<T, B, V>::BuildBound() {
    mes  = 0.;
    mesb = 0.;

    for (int i = 0; i < nt; i++)
        mes += this->elements[i].measure();

    for (int i = 0; i < nbe; i++)
        mesb += this->be(i).measure();
}

template class GenericMesh<Seg1, BoundaryPoint1, Vertex1>;
template class GenericMesh<Triangle2, BoundaryEdge2, Vertex2>;
template class GenericMesh<Tet, Triangle3, Vertex3>;
template class GenericMesh<Quad2, BoundaryEdge2, Vertex2>;
template class GenericMesh<Hexa, Quad3, Vertex3>;

template class BuildAdjacencyOfMesh<GenericMesh<Seg1, BoundaryPoint1, Vertex1>>;
template class BuildAdjacencyOfMesh<GenericMesh<Triangle2, BoundaryEdge2, Vertex2>>;
template class BuildAdjacencyOfMesh<GenericMesh<Tet, Triangle3, Vertex3>>;
template class BuildAdjacencyOfMesh<GenericMesh<Quad2, BoundaryEdge2, Vertex2>>;
template class BuildAdjacencyOfMesh<GenericMesh<Hexa, Quad3, Vertex3>>;
