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

#ifndef COMMON_GENERIC_VERTEX_HPP
#define COMMON_GENERIC_VERTEX_HPP

#include "Label.hpp"
#include "point.hpp"
/*
 *    Class of the generic vertex
 *
 */
template <typename Rn> class GenericVertex : public Rn, public Label {
    template <typename T, typename B, typename V> friend class GenericMesh;

    friend inline std::ostream &operator<<(std::ostream &f, const GenericVertex &v) {
        f << (const Rn &)v << ' ' << (const Label &)v;
        return f;
    }
    friend inline std::istream &operator>>(std::istream &f, GenericVertex &v) {
        f >> (Rn &)v >> (Label &)v;
        return f;
    }

  public:
    using value_type = typename Rn::value_type;

    typedef Rn Rd;
    static const int d = Rd::d;

    GenericVertex() : Rd(), Label(){};                         //,normal(0) {};
    GenericVertex(const Rd &P, int r = 0) : Rd(P), Label(r){}; //,normal(0){}

  private: // pas de copie pour ne pas prendre l'adresse
    GenericVertex(const GenericVertex &);
    void operator=(const GenericVertex &);
};

#endif
