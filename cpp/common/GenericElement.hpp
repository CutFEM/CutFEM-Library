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

#ifndef COMMON_GENERIC_ELEMENT_HPP
#define COMMON_GENERIC_ELEMENT_HPP
#include <cassert>
#include <vector>
#include "global.hpp"
#include "GenericVertex.hpp"

inline R1 ExtNormal(const std::array<GenericVertex<R1> *, 2> &v, const std::vector<int> &f) {
    return (f[0] == 0) ? R1(-1) : R1(1);
}

inline R2 ExtNormal(const std::array<GenericVertex<R2> *, 3> &v, const std::vector<int> &f) {
    return R2(*v[f.at(1)], *v[f.at(0)]).perp();
}
inline R2 ExtNormal(const std::array<GenericVertex<R2> *, 4> &v, const std::vector<int> &f) {
    return R2(*v[f.at(1)], *v[f.at(0)]).perp();
}

inline R3 ExtNormal(const std::array<GenericVertex<R3> *, 4> &v, const std::vector<int> &f) {
    return R3(*v[f.at(0)], *v[f.at(2)]) ^ R3(*v[f.at(0)], *v[f.at(1)]);
}

template <typename Data> class GenericElement : public Label {
  public:
    typedef typename Data::V Vertex;
    typedef typename Data::Face Face;
    typedef typename Data::V::Rd Rd;
    typedef typename Data::RdHat RdHat;
    typedef typename Data::RdHatBord RdHatBord;
    typedef typename Rd::R R;

    static const int nv              = Data::NbOfVertices;
    static const int ne              = Data::NbOfEdges;
    static const int nf              = Data::NbOfFaces;
    static const int nt              = Data::NbOfTet;
    static const int nitem           = nv + ne + nf + nt;
    static const int nva             = Data::NbOfVertexOnHyperFace;
    static const int nea             = Data::NbOfAdjElem;
    static const int d               = Rd::d;
    static const int nvOnFace        = Data::NvOnFace;
    static const int ParaviewNumCell = Data::ParaviewNumCell;
    static const int nb_ntcut        = Data::NbNtCut;
    static const int nvc             = Data::NbOfVerticesCut;
    static const int nb_sign_pattern = Data::NbSignPattern;

    static const std::vector<std::vector<int>> nvedge;
    static const std::vector<std::vector<int>> nvface;
    static const std::vector<std::vector<int>> edgeOfFace;
    static const std::vector<std::vector<int>> faceOfEdge;
    static const std::vector<std::vector<int>> nvhyperFace;
    static const std::vector<std::vector<int>> commonVertOfEdges;
    static std::vector<int> itemTopology() { return {nv, ne, nf, nt}; }
    static int oppVertOfEdge(int edge, int vert) { return vert == nvedge[edge][0] ? nvedge[edge][1] : nvedge[edge][0]; }

  public:
    std::array<Vertex *, nv> vertices;
    double mes;

  public:
    GenericElement() {}

    const Vertex &operator[](int i) const { return *vertices.at(i); }
    Vertex &operator[](int i) { return *vertices.at(i); }
    const Vertex &at(int i) const { return *vertices.at(i); }
    Vertex &at(int i) { return *vertices.at(i); }

    GenericElement &set(Vertex *v0, int *iv, int r, double mss = globalVariable::UnSetMesure) {
        for (int i = 0; i < nv; ++i)
            vertices[i] = v0 + iv[i];
        mes = (mss != globalVariable::UnSetMesure) ? mss : Data::mesure(vertices.data());
        lab = r;
        if (mss != globalVariable::UnSetMesure || mes <= 0) {
            for (int i = 0; i < nv; ++i)
                std::cout << *vertices[i] << std::endl;
            std::cout << mss << "\t" << mes << std::endl;
            getchar();
        }
        assert(mss == globalVariable::UnSetMesure && mes > 0);
        return *this;
    }

    void set_face(int ifac, Face &face) const {
        for (int i = 0; i < nva; ++i) {
            face.vertices[i] = vertices[nvface[ifac][i]];
        }
        face.mes = Data::mesure(face.vertices);
        face.lab = 0;
    }

    std::istream &Read1(std::istream &f, Vertex *v0, int n) {
        int iv[nv], ir, err = 0;
        for (int i = 0; i < nv; ++i) {
            f >> iv[i];
            iv[i]--;
            if (!(iv[i] >= 0 && iv[i] < n))
                err++;
        }
        f >> ir;
        if (err || !f.good()) {
            std::cerr << " Erreur GenericElement::Read1 " << nv << " " << n << "  : ";
            for (int j = 0; j < nv; ++j)
                std::cerr << iv[j] << " ";
            std::cerr << " , " << ir << std::endl;
            std::abort();
        }

        set(v0, iv, ir);
        return f;
    }

    /**
     * @brief Get edge i of element
     *
     * @param i
     * @return Rd Coordinates
     */
    Rd Edge(int i) const {
        assert(i >= 0 && i < ne);

        // nvedge[i] = numbering of edge i (e.g. (0,1), (1,2) or (2,0))
        // at(nvedge[i][0]) = Rd coordinates of vertex 0 of edge i, (x_0^i, y_0^i)
        // at(nvedge[i][1]) = Rd coordinates of vertex 1 of edge i, (x_1^i, y_1^i)

        return Rd(at(nvedge[i][0]), at(nvedge[i][1]));
        // (x_1^i-x_0^i, y_1^i-y_0^i) physical vector from vertex 0 to vertex 1 of edge i
    } // opposite edge vertex i

    Rd N(int i) const { return ExtNormal(vertices, nvhyperFace[i]) / (ExtNormal(vertices, nvhyperFace[i]).norm()); }
    Rd N_notNormalized(int i) const { return ExtNormal(vertices, nvhyperFace[i]); }

    Rd PBord(int i, RdHatBord P) const { return Data::PBord(nvhyperFace[i], P); }

    // THIS IS DIFFERENT FOR RECTANGLES
    virtual Rd operator()(const RdHat &Phat) const = 0;

    Rd barycenter() const {
        Rd Q;
        for (int i = 1; i < nv; ++i)
            Q += *(Rd *)vertices[i];
        return 1. / nv * Q;
    }

    int EdgeOrientation(int i) const { return 2 * (&at(nvedge[i][0]) < &at(nvedge[i][1])) - 1; }

    int edgePermutation(int i) const { return &at(nvedge[i][1]) < &at(nvedge[i][0]); } // 0 : no permutation

    R lenEdge(int i) const {
        assert(i >= 0 && i < ne);
        Rd E = Edge(i);
        return sqrt((E, E));
    }

    R hElement() const {
        double h = 0;
        for (int i = 0; i < ne; ++i)
            h += lenEdge(i);
        return h / ne;
    }
    R hMax() const {
        double h = 0;
        for (int i = 0; i < ne; ++i)
            h = max(h, lenEdge(i));
        return h;
    }
    R measure() const { return mes; }
    R get_h() const {
        double h = 0;
        for (int i = 0; i < ne; ++i)
            h += lenEdge(i);
        return h / ne;
    }
    Rd mapToPhysicalElement(const RdHat &Phat) const { return (*this)(Phat); }

    const auto begin() const { return vertices.begin(); }
    const auto end() const { return vertices.end(); }
    auto begin() { return vertices.begin(); }
    auto end() { return vertices.end(); }

  private:
    // pas de copie
    GenericElement(const GenericElement &);
    GenericElement &operator=(const GenericElement &);
};

#endif
