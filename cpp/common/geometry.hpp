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
#ifndef COMMON_GEOMETRY_HPP
#define COMMON_GEOMETRY_HPP

#include "Mesh2dn.hpp"

namespace geometry {

template <typename Rd> struct Segment {
   Rd A;
   Rd B;

   Segment(Rd a, Rd b) : A(a), B(b) {}
   bool is_between(const Rd C) const;
};

// equation cartesienne d'une droite
typedef struct coefDroite {
   double pente;
   double ord_or;
} Droite;

struct Interval {
   double v_min, v_max;
   Interval(double a, double b) : v_min(a), v_max(b) {}
   Interval operator*(double x) const { return Interval(x * v_min, x * v_max); }
   Interval operator+(double x) const { return Interval(x + v_min, x + v_max); }
};

Droite equation(const R2 a, const R2 b);
bool p_boncote(const R2 a, const R2 b, const R2 c, const R2 p);
bool p_dans_triangle(const typename Mesh2::Element &K, const R2 P);
bool p_dans_triangle(const R2 P, R2 a, R2 b, R2 c);
int find_triangle_contenant_p(const Mesh2 &Th, const R2 P, int k_init = 0);

R3 map_point_to_simplex(const R3 N[4], const R3 P);
R2 map_point_to_simplex(const R2 N[3], const R2 P);
// maps to face simplex
R3 map_point_to_simplex(const R3 N[3], const R2 P);
R2 map_point_to_simplex(const R2 N[2], const R1 P);

R measure_hyper_simplex(R2 N[2]);
R measure_hyper_simplex(R3 N[3]);

// change to template maybe
template <int d> R mesure_simplex(R2 N[d + 1]);
template <int d> R mesure_simplex(R3 N[d + 1]);
} // namespace geometry

#endif
