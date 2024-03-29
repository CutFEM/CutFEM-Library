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
#include "geometry.hpp"

double geometry::measure_hyper_simplex(R2 N[2]) {
   R2 u(N[0], N[1]);
   return u.norm();
}
double geometry::measure_hyper_simplex(R3 N[3]) {

   R3 u(N[0], N[1]);
   R3 v(N[0], N[2]);
   return Norme2(0.5 * (u ^ v));
}

template <> double geometry::mesure_simplex<1>(R2 N[2]) {
   return measure_hyper_simplex(N);
}
template <> double geometry::mesure_simplex<2>(R2 N[3]) {
   return std::fabs(det(N[0], N[1], N[2])) * 0.5;
}

template <> double geometry::mesure_simplex<1>(R3 N[2]) {
   R3 u(N[0], N[1]);
   return u.norm();
}
template <> double geometry::mesure_simplex<2>(R3 N[3]) {
   return measure_hyper_simplex(N);
}
template <> double geometry::mesure_simplex<3>(R3 N[4]) {
   R3 AB(N[0], N[1]);
   R3 AC(N[0], N[2]);
   R3 AD(N[0], N[3]);
   return std::fabs(det(AB, AC, AD)) / 6.;
}

R2 geometry::map_point_to_simplex(const R2 N[3], const R2 Phat) {
   R2 P = (1 - Phat.sum()) * N[0] + Phat[0] * N[1] + Phat[1] * N[2];
   return P;
}
R3 geometry::map_point_to_simplex(const R3 N[4], const R3 Phat) {
   R3 P = (1 - Phat.sum()) * N[0] + Phat[0] * N[1] + Phat[1] * N[2] +
          Phat[2] * N[3];
   return P;
}

R2 geometry::map_point_to_simplex(const R2 N[2], const R1 Phat) {
   R2 P = (1 - Phat[0]) * N[0] + Phat[0] * N[1];
   return P;
}
R3 geometry::map_point_to_simplex(const R3 N[3], const R2 Phat) {

   R3 P = (1 - Phat.sum()) * N[0] + Phat[0] * N[1] + Phat[1] * N[2];
   return P;
}

template <> bool geometry::Segment<R2>::is_between(const R2 C) const {
   assert(0);
   return true;
}

geometry::Droite geometry::equation(const R2 a, const R2 b) {
   Droite d;
   if (!(b.x - a.x))
      d.pente = b.x;
   else
      d.pente = (b.y - a.y) / (b.x - a.x);
   d.ord_or = a.y - a.x * d.pente;
   return d;
}

bool geometry::p_boncote(const R2 a, const R2 b, const R2 c, const R2 p) {
   // teste si la droite est perpendiculaire à l'axe des abcisses
   // et retourne le resultat en conséquence

   if (a.x == b.x) {
      return ((p.x - a.x) * (c.x - a.x) >= -globalVariable::Epsilon);
   }

   Droite d;
   d = equation(a, b);
   return ((d.pente * c.x + d.ord_or - c.y) *
               (d.pente * p.x + d.ord_or - p.y) >=
           -globalVariable::Epsilon);
}
bool geometry::p_dans_triangle(const typename Mesh2::Element &K, const R2 P) {
   return ((p_boncote(K[0], K[1], K[2], P)) &&
           (p_boncote(K[2], K[1], K[0], P)) &&
           (p_boncote(K[0], K[2], K[1], P)));
}
bool geometry::p_dans_triangle(const R2 P, R2 a, R2 b, R2 c) {
   return ((p_boncote(a, b, c, P)) && (p_boncote(c, b, a, P)) &&
           (p_boncote(a, c, b, P)));
}
int geometry::find_triangle_contenant_p(const Mesh2 &Th, const R2 P,
                                        int k_init) {

   int idx_elt = k_init;
   int count   = 0;

   // std::ofstream plot;
   // plot.open("node.dat",  std::ofstream::out);
   // plot<<  P << std::endl;
   // plot.close();
   // plot.open("searchElement.dat", std::ofstream::out);
   while (!p_dans_triangle(Th[idx_elt], P)) {

      const typename Mesh2::Element &K(Th[idx_elt]);
      // passe au triangle voisin qui est dans la "direction" de p
      // première arête

      int iface;
      if (!p_boncote(K[0], K[1], K[2], P)) {
         iface = 2;
      } else {
         if (!p_boncote(K[1], K[2], K[0], P)) {
            iface = 0;
         } else {
            iface = 1;
         }
      }
      idx_elt = Th.ElementAdj(idx_elt, iface);
   }
   return idx_elt;
}
