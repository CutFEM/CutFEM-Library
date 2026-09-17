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
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
 */

#include "cut_method.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

template <> int SignElement<Hexa>::nb_node_positif() const {
   int s = util::fsign(sum_);
   return 0.5 * (8 + sum_);
}

template <> RefPatch<Edge2>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;

   // byte ls[2];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
   //     if ( ls[0] == 0 && ls[1] == 0) continue;
   //     RefPatch<Edge2>::instance( ls);
   //
   //   }
   // }
}
template <> RefPatch<Triangle2>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;

   byte ls[3];
   for (ls[0] = -1; ls[0] < 2; ++ls[0]) {
      for (ls[1] = -1; ls[1] < 2; ++ls[1]) {
         for (ls[2] = -1; ls[2] < 2; ++ls[2]) {
            if (ls[0] == 0 && ls[1] == 0 && ls[2] == 0)
               continue;
            RefPatch<Triangle2>::instance(ls);
         }
      }
   }
}
template <> RefPatch<Triangle3>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;

   // byte ls[3];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
   //     for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
   //       if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0) continue;
   //       RefPatch<Triangle3>::instance( ls);
   //     }
   //   }
   // }
}
template <> RefPatch<Quad2>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;

   byte ls[4];
   for (ls[0] = -1; ls[0] < 2; ++ls[0]) {
      for (ls[1] = -1; ls[1] < 2; ++ls[1]) {
         for (ls[2] = -1; ls[2] < 2; ++ls[2]) {
            for (ls[3] = -1; ls[3] < 2; ++ls[3]) {
               if (ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
                  continue;
               RefPatch<Quad2>::instance(ls);
            }
         }
      }
   }
}
template <> RefPatch<Tet>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;

   byte ls[4];
   for (ls[0] = -1; ls[0] < 2; ++ls[0]) {
      for (ls[1] = -1; ls[1] < 2; ++ls[1]) {
         for (ls[2] = -1; ls[2] < 2; ++ls[2]) {
            for (ls[3] = -1; ls[3] < 2; ++ls[3]) {
               if (ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
                  continue;
               RefPatch<Tet>::instance(ls);
            }
         }
      }
   }
}
template <> RefPatch<Hexa>::InitializerCL::InitializerCL() {
   init_count_++;
   if (init_count_++ > 0)
      return;
}

// RefPatch
template <> bool RefPatch<Edge2>::assign(const SignPattern<Edge2> &cut) {
   for (size_ = 0; size_ < num_elements(cut); ++size_) {
      face_[size_] = FaceIdx(cut(size_));
   }
   return empty();
}
template <>
bool RefPatch<Triangle2>::assign(const SignPattern<Triangle2> &cut) {

   for (size_ = 0; size_ < num_elements(cut); ++size_) {
      face_[size_] = FaceIdx(cut(size_), cut(size_ + 1));
   }
   return empty();
}
template <>
bool RefPatch<Triangle3>::assign(const SignPattern<Triangle3> &cut) {

   for (size_ = 0; size_ < num_elements(cut); ++size_) {
      face_[size_] = FaceIdx(cut(size_), cut(size_ + 1));
   }
   return empty();
}
template <> bool RefPatch<Quad2>::assign(const SignPattern<Quad2> &cut) {

   for (size_ = 0; size_ < num_elements(cut); ++size_) {
      face_[size_] = FaceIdx(cut(size_), cut(size_ + 1));
   }
   return empty();
}
template <> bool RefPatch<Tet>::assign(const SignPattern<Tet> &cut) {
   for (size_ = 0; size_ < num_elements(cut); ++size_) {
      face_[size_] = FaceIdx(cut(size_), cut(size_ + 1), cut(size_ + 2));
   }
   return empty();
}
template <>
Ubyte RefPatch<Hexa>::num_elements(const SignPattern<Hexa> &cut) const {
   switch (cut.num_cut_simplexes()) {
   case 3:
      return 1;
   case 4:
      return 2;
   case 5:
      return 3;
   case 6:
      return 4;
   default:
      return 0;
   }
}
template <> bool RefPatch<Hexa>::assign(const SignPattern<Hexa> &cut) {
   switch (static_cast<int>(cut.num_cut_simplexes())) {
   case 3:
      face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
      break;
   case 4:
      face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
      face_[size_++] = FaceIdx(cut(1), cut(2), cut(3));
      break;
   case 5:
      face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
      face_[size_++] = FaceIdx(cut(0), cut(2), cut(3));
      face_[size_++] = FaceIdx(cut(0), cut(3), cut(4));
      break;
   case 6:
      face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
      face_[size_++] = FaceIdx(cut(0), cut(2), cut(3));
      face_[size_++] = FaceIdx(cut(0), cut(3), cut(5));
      face_[size_++] = FaceIdx(cut(3), cut(4), cut(5));
      break;
   }
   return empty();
}

// RefPartition
template <> RefPartition<Edge2>::InitializerCL::InitializerCL() {
   // std::cout << " initialize instances Triangle" << std::endl;
   if (init_count_++ > 0)
      return;
   //
   // byte ls[2];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
   //     if ( ls[0] == 0 && ls[1] == 0 ) continue;
   //     RefPartition<Edge2>::instance( ls);
   //
   //   }
   // }
}
template <> RefPartition<Triangle2>::InitializerCL::InitializerCL() {
   // std::cout << " initialize instances Triangle" << std::endl;

   if (init_count_++ > 0)
      return;

   byte ls[3];
   for (ls[0] = -1; ls[0] < 2; ++ls[0]) {
      for (ls[1] = -1; ls[1] < 2; ++ls[1]) {
         for (ls[2] = -1; ls[2] < 2; ++ls[2]) {
            if (ls[0] == 0 && ls[1] == 0 && ls[2] == 0)
               continue;
            RefPartition<Triangle2>::instance(ls);
         }
      }
   }
}
template <> RefPartition<Triangle3>::InitializerCL::InitializerCL() {
   return;
   // if (init_count_++ > 0) return;
   // byte ls[3];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
   //     for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
   //       if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0) continue;
   //       RefPartition<Triangle2>::instance( ls);
   //     }
   //   }
   // }
}
template <> RefPartition<Quad2>::InitializerCL::InitializerCL() {
   // std::cout << " initialize instances Quad" << std::endl;
   if (init_count_++ > 0)
      return;
   // if (init_count_++ > 0) return;
   //
   // byte ls[4];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]){
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]){
   //     for (ls[2]= -1; ls[2] < 2; ++ls[2]){
   //       for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
   //
   //         if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
   //         continue; RefPartition<Quad2>::instance( ls);
   //       }
   //     }
   //   }
   // }
}
template <> RefPartition<Tet>::InitializerCL::InitializerCL() {

   if (init_count_++ > 0)
      return;
   //
   // byte ls[4];
   // for (ls[0]= -1; ls[0] < 2; ++ls[0]){
   //   for (ls[1]= -1; ls[1] < 2; ++ls[1]){
   //     for (ls[2]= -1; ls[2] < 2; ++ls[2]){
   //       for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
   //         if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
   //         continue; RefPartition<Tet>::instance( ls);
   //       }
   //     }
   //   }
   // }
}
template <> RefPartition<Hexa>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;
}

// Gustaf
template <> RefPartition<Quad3>::InitializerCL::InitializerCL() {
   if (init_count_++ > 0)
      return;
}

template <> bool RefPartition<Edge2>::assign(const SignPattern<Edge2> &cut) {
   end_ = begin_ = elements_ + start_array;

   if (cut.empty()) { // Most common case: no cut.
      Ubyte list_v[] = {0, 1};
      AddElement(list_v, cut.sign(0));
   } else if (cut.no_zero_vertex()) { // next common case: cuts without vertices
                                      // on the zero level
      is_cut_        = true;
      Ubyte list_v[] = {0, cut(0)};
      AddElement(list_v, cut.sign(0));
      Ubyte list_v1[] = {1, cut(0)};
      AddElement(list_v1, cut.sign(1));
   }
   return is_uncut();
}
template <>
bool RefPartition<Triangle2>::assign(const SignPattern<Triangle2> &cut) {
   end_ = begin_ = elements_ + start_array;

   if (cut.empty()) { // Most common case: no cut.
      Ubyte list_v[] = {0, 1, 2};
      // AddTriangle( 0, 1, 2, cut.sign( 0));
      AddElement(list_v, cut.sign(0));

   } else if (cut.no_zero_vertex()) { // next common case: cuts without vertices
                                      // on the zero level
      is_cut_        = true;
      const Ubyte v  = Triangle2::commonVertOfEdges[cut[0]][cut[1]];
      Ubyte list_v[] = {v, cut(0), cut(1)};
      AddElement(list_v, cut.sign(v));
      // AddTriangle( v, cut(0), cut(1), cut.sign( v));
      AddQuadrilateral(Triangle2::oppVertOfEdge(cut[0], v), cut(0),
                       Triangle2::oppVertOfEdge(cut[1], v), cut(1),
                       -cut.sign(v));
   }
   // else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next
   // common case: there are cut edges, and also 1 or 2 vertices of the tetra
   // with value 0 (the latter as we are in the else-part of
   // cut.no_zero_vertex())
   //     if (cut.num_zero_vertexes() == 1) { // triangular cut through a
   //     vertex: a tetra and a remaining pyramid with quadrilateral base
   //         const Ubyte e= cut[1], f= cut[2];
   //         const Ubyte v= VertByEdge( e, f);
   //         AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
   //         const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1)
   //         : VertOfEdge( e, 0); const Ubyte opp_v_in_f= v == VertOfEdge( f,
   //         0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
   //         // the pyramid
   //         AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
   //         AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
   //     }
   //     else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2
   //     vertexes: two tetras
   //         const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
   //         const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
   //         AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
   //         AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
   //     }
   // }
   // else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
   //     AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
   return is_uncut();
}
template <>
bool RefPartition<Triangle3>::assign(const SignPattern<Triangle3> &cut) {
   end_ = begin_ = elements_ + start_array;

   if (cut.empty()) { // Most common case: no cut.
      Ubyte list_v[] = {0, 1, 2};
      // AddTriangle( 0, 1, 2, cut.sign( 0));
      AddElement(list_v, cut.sign(0));

   } else if (cut.no_zero_vertex()) { // next common case: cuts without vertices
                                      // on the zero level
      is_cut_        = true;
      const Ubyte v  = Triangle3::commonVertOfEdges[cut[0]][cut[1]];
      Ubyte list_v[] = {v, cut(0), cut(1)};
      AddElement(list_v, cut.sign(v));
      // AddTriangle( v, cut(0), cut(1), cut.sign( v));
      AddQuadrilateral(Triangle3::oppVertOfEdge(cut[0], v), cut(0),
                       Triangle3::oppVertOfEdge(cut[1], v), cut(1),
                       -cut.sign(v));
   }
   // else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next
   // common case: there are cut edges, and also 1 or 2 vertices of the tetra
   // with value 0 (the latter as we are in the else-part of
   // cut.no_zero_vertex())
   //     if (cut.num_zero_vertexes() == 1) { // triangular cut through a
   //     vertex: a tetra and a remaining pyramid with quadrilateral base
   //         const Ubyte e= cut[1], f= cut[2];
   //         const Ubyte v= VertByEdge( e, f);
   //         AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
   //         const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1)
   //         : VertOfEdge( e, 0); const Ubyte opp_v_in_f= v == VertOfEdge( f,
   //         0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
   //         // the pyramid
   //         AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
   //         AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
   //     }
   //     else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2
   //     vertexes: two tetras
   //         const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
   //         const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
   //         AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
   //         AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
   //     }
   // }
   // else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
   //     AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
   return is_uncut();
}
template <> bool RefPartition<Quad2>::assign(const SignPattern<Quad2> &cut) {
   end_ = begin_ = elements_ + start_array;
   if (cut.empty()) { // Most common case: no cut.
      AddQuadrilateral(0, 1, 2, 3, cut.sign(0));
   } else if (cut.no_zero_vertex()) { // next common case: cuts without vertices
                                      // on the zero level
      // case when interface cut two adjacent edges
      is_cut_ = true;

      if (cut.sign_element.sum() != 0) {
         const Ubyte v  = Quad2::commonVertOfEdges[cut[0]][cut[1]];
         Ubyte list_v[] = {v, cut(0), cut(1)};
         AddElement(list_v, cut.sign(v));
         // AddTriangle (v, cut(0), cut(1), cut.sign( v));
         const Ubyte v_op = (v + 2) % 4;
         AddPentagone(cut(0), Quad2::oppVertOfEdge(cut[0], v), cut(1),
                      Quad2::oppVertOfEdge(cut[1], v), v_op, -cut.sign(v));
      } else {
         const Ubyte v0 = Quad2::nvedge[cut[0]][0];
         const Ubyte u0 = Quad2::nvedge[cut[0]][1];
         const Ubyte v1 = (u0 + 2) % 4;
         const Ubyte u1 = (v0 + 2) % 4;
         AddQuadrilateral(v0, cut(0), v1, cut(1), cut.sign(v0));
         AddQuadrilateral(u1, cut(1), u0, cut(0), -cut.sign(v0));
      }
   } else {
   }

   return is_uncut();
}
template <> bool RefPartition<Tet>::assign(const SignPattern<Tet> &cut) {
   end_ = begin_ = elements_ + start_array;

   if (cut.empty()) { // Most common case: no cut.
      Ubyte list_v[] = {0, 1, 2, 3};
      AddElement(list_v, cut.sign(0));
      // AddTetra( 0, 1, 2, 3, cut.sign( 0));
   } else if (cut.no_zero_vertex()) { // next common case: cuts without vertices
                                      // on the zero level
      is_cut_ = true;

      if (cut.num_cut_simplexes() ==
          3) { // triangular cut: a tetra and a remaining prism

         const Ubyte v  = Tet::commonVertOfEdges[cut[0]][cut[1]];
         Ubyte list_v[] = {v, cut(0), cut(1), cut(2)};
         AddElement(list_v, cut.sign(v));
         // AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));

         AddPrism(Tet::oppVertOfEdge(cut[0], v), cut(0),
                  Tet::oppVertOfEdge(cut[1], v), cut(1),
                  Tet::oppVertOfEdge(cut[2], v), cut(2), -cut.sign(v));
      } else if (cut.num_cut_simplexes() ==
                 4) { // quadrilateral cut: two prisms

         const Ubyte e = first_uncut_edge(cut);
         const Ubyte f = Tet::oppEdgeOfEdge[e];
         AddPrism(Tet::nvedge[e][0], Tet::nvedge[e][1], cut(0), cut(2), cut(1),
                  cut(3), cut.sign(Tet::nvedge[e][0]));
         AddPrism(Tet::nvedge[f][0], Tet::nvedge[f][1], cut(0), cut(1), cut(2),
                  cut(3), cut.sign(Tet::nvedge[f][0]));
      }
   } else {
   }
   // else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next
   // common case: there are cut edges, and also 1 or 2 vertices of the tetra
   // with value 0 (the latter as we are in the else-part of
   // cut.no_zero_vertex())
   //     if (cut.num_zero_vertexes() == 1) { // triangular cut through a
   //     vertex: a tetra and a remaining pyramid with quadrilateral base
   //         const Ubyte e= cut[1], f= cut[2];
   //         const Ubyte v= VertByEdge( e, f);
   //         AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
   //         const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1)
   //         : VertOfEdge( e, 0); const Ubyte opp_v_in_f= v == VertOfEdge( f,
   //         0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
   //         // the pyramid
   //         AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
   //         AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
   //     }
   //     else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2
   //     vertexes: two tetras
   //         const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
   //         const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
   //         AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
   //         AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
   //     }
   // }
   // else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
   //     AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
   return is_uncut();
}
namespace {

// ---------------------------------------------------------------------------
//  Validation for the hexahedral reference partition.
//
//  RefPartition<Hexa>::assign implements only single-sheet cut patterns: 3, 4,
//  5 or 6 cut edges forming one loop, with a supported distribution of positive
//  vertices, and no vertex exactly on the zero set.  Everything outside that
//  set used to be caught by assert() alone, which production builds remove with
//  -DNDEBUG.  An unsupported pattern then read Hexa::commonVertOfEdges[..][..]
//  == -1 as a vertex index and wrote through it, which is how both 3D active
//  surface production jobs died with SIGSEGV on 2026-09-12/13 while writing a
//  ParaView file of a fragmented interface.
//
//  These checks are always compiled.  They cost a handful of comparisons per
//  cut cell and turn a silent buffer overrun into a diagnosable exception that
//  names the pattern, which a caller such as a visualization writer can catch
//  and skip.  Supported patterns behave exactly as before.
// ---------------------------------------------------------------------------
std::string describe(const SignPattern<Hexa> &cut) {
   std::ostringstream description;
   description << "cut edges=" << static_cast<int>(cut.num_cut_simplexes())
               << ", zero vertices=" << static_cast<int>(cut.num_zero_vertexes())
               << ", positive vertices=" << cut.sign_element.nb_node_positif()
               << ", vertex signs=[";
   for (int v = 0; v < Hexa::nv; ++v)
      description << (v == 0 ? "" : " ") << static_cast<int>(cut.sign(v));
   description << "], cut edge ids=[";
   for (int i = 0; i < static_cast<int>(cut.num_cut_simplexes()); ++i)
      description << (i == 0 ? "" : " ") << static_cast<int>(cut[i]);
   description << "]";
   return description.str();
}

[[noreturn]] void reject(const SignPattern<Hexa> &cut, const std::string &what) {
   throw std::runtime_error(
       "RefPartition<Hexa>: unsupported cut configuration (" + what + "); " +
       describe(cut) +
       ". The trilinear reference partition only tessellates single-sheet cuts; "
       "a cell whose interface has two sheets, a vertex on the zero set, or a "
       "cut-edge loop outside the supported table cannot be represented. The "
       "Algoim MultiPoly quadrature is unaffected: this path is used by the "
       "ParaView writers and by linear-interface assembly, not by the cut rules.");
}

// Connectivity lookups return -1 when the two edges share no vertex, which
// happens exactly for the patterns the table below cannot express.
int require_vertex(int vertex, const SignPattern<Hexa> &cut, const char *what) {
   if (vertex < 0 || vertex >= Hexa::nv)
      reject(cut, what);
   return vertex;
}

} // namespace

template <> bool RefPartition<Hexa>::assign(const SignPattern<Hexa> &cut) {
   end_ = begin_ = elements_ + start_array;
   if (cut.empty()) { // Most common case: no cut.
      AddHexa(0, 1, 3, 2, 4, 5, 7, 6, cut.sign(0));
   } else {
      is_cut_ = true;

      // Always-on preconditions; see the note above the helpers.
      if (!cut.no_zero_vertex())
         reject(cut, "a vertex lies exactly on the zero set");
      int nb_pos = cut.sign_element.nb_node_positif();

      switch (static_cast<int>(cut.num_cut_simplexes())) {
      case 3: {
         // std::cout << " case \t 3" <<std::endl;
         // only one option, 1 tetra + 1 prisme + complet with 4 tetra
         // same way we cut a hexa into tetra (1tetra+prism => one fat tetra)
         int v0 = require_vertex(Hexa::commonVertOfEdges[cut[0]][cut[1]], cut,
                                 "the three cut edges do not meet in one vertex");
         Ubyte list_v[] = {static_cast<Ubyte>(v0), cut(0), cut(1), cut(2)};
         AddElement(list_v, cut.sign(v0));
         AddPrism(cut(0), Hexa::oppVertOfEdge(cut[0], v0), cut(1),
                  Hexa::oppVertOfEdge(cut[1], v0), cut(2),
                  Hexa::oppVertOfEdge(cut[2], v0), -cut.sign(v0));

         int e_op0 = Hexa::oppEdgeOfEdge[cut[0]];
         int e_op1 = Hexa::oppEdgeOfEdge[cut[1]];
         int e_op2 = Hexa::oppEdgeOfEdge[cut[2]];
         Ubyte u1  = static_cast<Ubyte>(require_vertex(
             Hexa::commonVertOfEdges[e_op0][e_op1], cut,
             "the opposite edges of a 3-edge cut do not meet in one vertex"));

         // sommet des 3 autres tetra de coin
         Ubyte s0 = Hexa::oppVertOfEdge(e_op0, u1);
         Ubyte s1 = Hexa::oppVertOfEdge(e_op1, u1);
         Ubyte s2 = Hexa::oppVertOfEdge(e_op2, u1);

         // center tetra
         Ubyte list_v5[] = {
             static_cast<Ubyte>(u1),
             static_cast<Ubyte>(Hexa::oppVertOfEdge(cut[0], v0)),
             static_cast<Ubyte>(Hexa::oppVertOfEdge(cut[1], v0)),
             static_cast<Ubyte>(Hexa::oppVertOfEdge(cut[2], v0))};
         AddElement(list_v5, cut.sign(s0));

         Ubyte list_v1[] = {static_cast<Ubyte>(s0),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s0][0]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s0][1]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s0][2])};
         AddElement(list_v1, cut.sign(s0));
         Ubyte list_v2[] = {static_cast<Ubyte>(s1),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s1][0]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s1][1]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s1][2])};
         AddElement(list_v2, cut.sign(s0));
         Ubyte list_v3[] = {static_cast<Ubyte>(s2),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s2][0]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s2][1]),
                            static_cast<Ubyte>(Hexa::nodeConnectivity[s2][2])};
         AddElement(list_v3, cut.sign(s0));

      } break;
      case 4:
         // std::cout << " case \t 4" <<std::endl;
         // add 1 prisme in part with 2 nodes
         // add 1 prisme && 1 hexa in part with 6 nodes
         if (nb_pos == 2 || nb_pos == 6) {
            // Prisme => add 3 edges connecting the 2 triangular faces
            int e0 = 0;
            int e1 = (Hexa::commonVertOfEdges[cut[0]][cut[1]] == -1) ? 1 : 3;
            int f0 = (Hexa::commonVertOfEdges[cut[0]][cut[1]] == -1) ? 3 : 1;
            int f1 = 2;
            int g0 = require_vertex(Hexa::commonVertOfEdges[cut[e0]][cut[f0]], cut,
                                    "4-edge cut: edges e0 and f0 do not meet");
            int g1 = require_vertex(Hexa::commonVertOfEdges[cut[e1]][cut[f1]], cut,
                                    "4-edge cut: edges e1 and f1 do not meet");
            AddPrism(cut(e0), cut(e1), cut(f0), cut(f1), g0, g1, cut.sign(g0));

            // add Hexa
            // send a cartesian index way
            int h0 = Hexa::oppVertOfEdge(cut[e0], g0);
            int h1 = Hexa::oppVertOfEdge(cut[e1], g1);
            int l0 = Hexa::oppVertOfEdge(cut[f0], g0);
            int l1 = Hexa::oppVertOfEdge(cut[f1], g1);
            AddHexa(cut(e0), cut(e1), cut(f0), cut(f1), h0, h1, l0, l1,
                    -cut.sign(g0));

            int ed_op00 = Hexa::oppEdgeOfEdge[cut[e1]];
            int ed_op01 = Hexa::oppEdgeOfEdge[cut[f1]];
            int ed_op10 = Hexa::oppEdgeOfEdge[cut[e0]];
            int ed_op11 = Hexa::oppEdgeOfEdge[cut[f0]];
            int k0      = require_vertex(Hexa::commonVertOfEdges[ed_op00][ed_op01],
                                         cut, "4-edge cut: opposite edges do not meet (k0)");
            int k1      = require_vertex(Hexa::commonVertOfEdges[ed_op10][ed_op11],
                                         cut, "4-edge cut: opposite edges do not meet (k1)");
            AddPrism(h0, h1, l0, l1, k0, k1, -cut.sign(g0));
         }
         // add 1 hexa on each part
         else if (nb_pos == 4) {
            Rn n_pos(4), n_neg(4);
            for (int i = 0; i < 4; ++i) {
               int v0   = Hexa::nvedge[cut[i]][0];
               int v1   = Hexa::nvedge[cut[i]][1];
               int s0   = (int)cut.sign(v0);
               n_neg(i) = (s0 < 0) ? v0 : v1;
               n_pos(i) = (s0 > 0) ? v0 : v1;
            }
            AddHexa(cut(0), cut(1), cut(3), cut(2), n_pos(0), n_pos(1),
                    n_pos(3), n_pos(2), 1);
            AddHexa(cut(0), cut(1), cut(3), cut(2), n_neg(0), n_neg(1),
                    n_neg(3), n_neg(2), -1);
         } else {
            reject(cut, "4-edge cut with an unsupported number of positive "
                        "vertices (expected 2, 4 or 6)");
         }
         break;
      case 5:
         // std::cout << " case \t 5" <<std::endl;
         {
            int nb_pos = cut.sign_element.nb_node_positif();
            if (nb_pos != 3 && nb_pos != 5)
               reject(cut, "5-edge cut with an unsupported number of positive "
                           "vertices (expected 3 or 5)");
            Rn v(8);
            v = -1;
            Rn e(8);
            e      = -1;
            // find starting cut node
            int ii = 0;
            for (int i = 0; i < 5; ++i) {
               // std::cout << (int) cut[i] << std::endl;
               int i_next = (i + 1) % 5;
               int i_prev = (i == 0) ? 4 : i - 1;
               if (Hexa::commonVertOfEdges[cut[i]][cut[i_next]] == -1 &&
                   Hexa::commonVertOfEdges[cut[i]][cut[i_prev]] == -1) {
                  e(0) = i;
                  e(1) = i_next;
                  e(2) = i_prev;
               } else if (Hexa::commonVertOfEdges[cut[i]][cut[i_next]] != -1 &&
                          Hexa::commonVertOfEdges[cut[i]][cut[i_prev]] != -1) {
                  e(3 + ii) = i;
                  ii += 1;
               }
            }
            // std::cout << std::endl;
            // for(int i=0;i<5;++i) std::cout << (int) e[i] << "\t";
            // std::cout << std::endl;

            if (e(0) == -1 || e(4) == -1)
               reject(cut, "5-edge cut whose edges do not form the single "
                           "supported loop");
            if (Hexa::commonVertOfEdges[cut[e(2)]][cut[e(3)]] != -1) {
               int eee = e(3);
               e(3)    = e(4);
               e(4)    = eee;
            }
            // find 3 nodes in pos part
            {
               int vv0 = Hexa::nvedge[cut[e(0)]][0];
               int ss  = (int)cut.sign(vv0);
               v(0)    = (ss == -1) ? Hexa::nvedge[cut[e(0)]][0]
                                    : Hexa::nvedge[cut[e(0)]][1];
               v(4)    = (ss == +1) ? Hexa::nvedge[cut[e(0)]][0]
                                    : Hexa::nvedge[cut[e(0)]][1];
               if (nb_pos == 5)
                  std::swap(v(0), v(4));
            }
            {
               int vv0 = Hexa::nvedge[cut[e(1)]][0];
               int ss  = (int)cut.sign(vv0);
               v(1)    = (ss == -1) ? Hexa::nvedge[cut[e(1)]][0]
                                    : Hexa::nvedge[cut[e(1)]][1];
               v(5)    = (ss == +1) ? Hexa::nvedge[cut[e(1)]][0]
                                    : Hexa::nvedge[cut[e(1)]][1];
               if (nb_pos == 5)
                  std::swap(v(1), v(5));
            }
            {
               int vv0 = Hexa::nvedge[cut[e(2)]][0];
               int ss  = (int)cut.sign(vv0);
               v(3)    = (ss == -1) ? Hexa::nvedge[cut[e(2)]][0]
                                    : Hexa::nvedge[cut[e(2)]][1];
               v(7)    = (ss == +1) ? Hexa::nvedge[cut[e(2)]][0]
                                    : Hexa::nvedge[cut[e(2)]][1];
               if (nb_pos == 5)
                  std::swap(v(3), v(7));
            }
            {
               v(6)     = require_vertex(
                   Hexa::commonVertOfEdges[cut[static_cast<int>(e(3))]]
                                          [cut[static_cast<int>(e(4))]],
                   cut, "5-edge cut: edges e3 and e4 do not meet");
               int e_op = Hexa::oppEdgeOfEdge[cut[e(0)]];
               v(2)     = Hexa::oppVertOfEdge(e_op, v(6));
            }

            // add prisme 1
            AddPrism(cut(e(0)), v(4), cut(e(1)), v(5), cut(e(2)), v(7),
                     cut.sign(v(7)));
            AddPrism(v(5), v(7), cut(e(1)), cut(e(2)), cut(e(3)), cut(e(4)),
                     cut.sign(v(5)));
            AddPrism(cut(e(0)), v(0), cut(e(1)), v(1), cut(e(2)), v(3),
                     cut.sign(v(0)));
            AddPrism(cut(e(4)), cut(e(2)), cut(e(3)), cut(e(1)), v(6), v(2),
                     cut.sign(v(2)));

            Ubyte list_v2[] = {static_cast<Ubyte>(v(1)), cut(e(1)), cut(e(2)),
                               static_cast<Ubyte>(v(2))};
            AddElement(list_v2, cut.sign(v(2)));
            Ubyte list_v3[] = {static_cast<Ubyte>(v(3)),
                               static_cast<Ubyte>(v(1)), cut(e(2)),
                               static_cast<Ubyte>(v(2))};
            AddElement(list_v3, cut.sign(v(2)));
         }

         break;
      case 6:
         // std::cout << " case \t 6" <<std::endl;
         {
            KN<int> v(8);
            v    = -1;
            // A 6-edge cut is only tessellable when the six edges form one
            // hexagonal ring.  Two disconnected triangles (opposite corners cut
            // off) also give six cut edges and land here; every lookup below
            // then returns -1.  This is the configuration a fragmenting
            // interface produces.
            v(0) = require_vertex(Hexa::commonVertOfEdges[cut[1]][cut[2]], cut,
                                  "6-edge cut: edges 1 and 2 do not meet");
            v(1) = require_vertex(Hexa::commonVertOfEdges[cut[2]][cut[3]], cut,
                                  "6-edge cut: edges 2 and 3 do not meet");
            v(3) = require_vertex(Hexa::commonVertOfEdges[cut[0]][cut[1]], cut,
                                  "6-edge cut: edges 0 and 1 do not meet");
            v(5) = require_vertex(Hexa::commonVertOfEdges[cut[3]][cut[4]], cut,
                                  "6-edge cut: edges 3 and 4 do not meet");
            v(6) = require_vertex(Hexa::commonVertOfEdges[cut[4]][cut[5]], cut,
                                  "6-edge cut: edges 4 and 5 do not meet");
            v(7) = require_vertex(Hexa::commonVertOfEdges[cut[0]][cut[5]], cut,
                                  "6-edge cut: edges 0 and 5 do not meet");
            AddPrism(cut(0), cut(3), cut(5), cut(4), v(7), v(5),
                     cut.sign(v(5)));
            AddPrism(cut(1), cut(2), cut(0), cut(3), v(3), v(1),
                     cut.sign(v(3)));

            { // get the two last nodes
               int s0 = cut.sign(v(0));
               for (int i = 0; i < 3; ++i) {
                  int t0 = Hexa::nodeConnectivity[v(0)][i];
                  if (cut.sign(t0) == s0) {
                     v(4) = t0;
                     break;
                  }
               }
               int s6 = cut.sign(v(6));
               for (int i = 0; i < 3; ++i) {
                  int t0 = Hexa::nodeConnectivity[v(6)][i];
                  if (cut.sign(t0) == s6) {
                     v(2) = t0;
                     break;
                  }
               }
            }
            AddPrism(cut(1), cut(0), cut(2), cut(3), v(0), v(4),
                     cut.sign(v(0)));
            AddPrism(cut(4), cut(3), cut(5), cut(0), v(6), v(2),
                     cut.sign(v(2)));
            Ubyte list_v0[] = {static_cast<Ubyte>(v(5)), cut(3), cut(0),
                               static_cast<Ubyte>(v(4))};
            AddElement(list_v0, cut.sign(v(4)));
            Ubyte list_v1[] = {static_cast<Ubyte>(v(5)),
                               static_cast<Ubyte>(v(7)), cut(0),
                               static_cast<Ubyte>(v(4))};
            AddElement(list_v1, cut.sign(v(4)));

            Ubyte list_v2[] = {static_cast<Ubyte>(v(1)), cut(3), cut(0),
                               static_cast<Ubyte>(v(2))};
            AddElement(list_v2, cut.sign(v(2)));
            Ubyte list_v3[] = {static_cast<Ubyte>(v(3)),
                               static_cast<Ubyte>(v(1)), cut(0),
                               static_cast<Ubyte>(v(2))};
            AddElement(list_v3, cut.sign(v(2)));
         }
         break;
      default:
         reject(cut, "unsupported number of cut edges (the table covers 3 to 6)");
      };
   }
   return is_uncut();
}
// Gustaf
template <>
bool RefPartition<Quad3>::assign(const SignPattern<Quad3> &cut) {
   end_ = begin_ = elements_ + start_array;

   if (cut.empty()) {
      // Most common case: no cut.
      AddQuadrilateral(0, 1, 2, 3, cut.sign(0));

   } else if (cut.no_zero_vertex()) {
      is_cut_ = true;

      if (cut.sign_element.sum() != 0) {
         // Interface cuts two adjacent edges:
         // one triangle on one side and one pentagon on the other.
         const Ubyte v  = Quad3::commonVertOfEdges[cut[0]][cut[1]];
         Ubyte list_v[] = {v, cut(0), cut(1)};
         AddElement(list_v, cut.sign(v));

         const Ubyte v_op = (v + 2) % 4;
         AddPentagone(cut(0), Quad3::oppVertOfEdge(cut[0], v),
                      cut(1), Quad3::oppVertOfEdge(cut[1], v),
                      v_op, -cut.sign(v));
      } else {
         // Interface cuts opposite edges:
         // split into two quadrilaterals.
         const Ubyte v0 = Quad3::nvedge[cut[0]][0];
         const Ubyte u0 = Quad3::nvedge[cut[0]][1];
         const Ubyte v1 = (u0 + 2) % 4;
         const Ubyte u1 = (v0 + 2) % 4;

         AddQuadrilateral(v0, cut(0), v1, cut(1), cut.sign(v0));
         AddQuadrilateral(u1, cut(1), u0, cut(0), -cut.sign(v0));
      }
   }

   return is_uncut();
}

template <>
Ubyte RefPartition<Hexa>::first_uncut_edge(const SignPattern<Hexa> &cut) const {
   assert(0);
   return static_cast<Ubyte>(0);
}
template <>
void Partition<Edge2>::get_list_node(std::vector<R2> &node, int s) const {
   assert(0);
}
template <>
void Partition<Triangle2>::get_list_node(std::vector<R2> &node, int s) const {
   node.resize(0);
   const SignPattern<Triangle2> cut(ls);
   if (cut.empty() || s == 0) {
      for (int i = 0; i < 3; ++i)
         node.push_back(T[i]);
   } else {
      const Ubyte v = Triangle2::commonVertOfEdges[cut[0]][cut[1]];
      if (cut.sign(v) == s) { // add triangle v, cut(0), cut(1)
         node.push_back(T[v]);
         node.push_back(get_vertex(cut(0)));
         node.push_back(get_vertex(cut(1)));
      } else {
         node.push_back(get_vertex(Triangle2::oppVertOfEdge(cut[0], v)));
         node.push_back(get_vertex(cut(0)));
         node.push_back(get_vertex(cut(1)));
         node.push_back(get_vertex(Triangle2::oppVertOfEdge(cut[1], v)));
      }
   }
}
template <>
void Partition<Quad2>::get_list_node(std::vector<R2> &node, int s) const {
   node.resize(0);
   const SignPattern<Quad2> cut(ls);
   if (cut.empty() || s == 0) {
      for (int i = 0; i < 4; ++i)
         node.push_back(T[i]);
   } else {
      if (cut.sign_element.sum() != 0) {
         const Ubyte v = Quad2::commonVertOfEdges[cut[0]][cut[1]];
         if (cut.sign(v) == s) { // add triangle v, cut(0), cut(1)
            node.push_back(T[v]);
            node.push_back(get_vertex(cut(0)));
            node.push_back(get_vertex(cut(1)));
         } else { // add pentagone
            const Ubyte v_op = (v + 2) % 4;
            node.push_back(get_vertex(Quad2::oppVertOfEdge(cut[0], v)));
            node.push_back(get_vertex(cut(0)));
            node.push_back(get_vertex(cut(1)));
            node.push_back(get_vertex(Quad2::oppVertOfEdge(cut[1], v)));
            node.push_back(get_vertex(v_op));
         }
      } else {

         const Ubyte v0 = Quad2::nvedge[cut[0]][0];
         const Ubyte u0 = Quad2::nvedge[cut[0]][1];
         const Ubyte v1 = (u0 + 2) % 4;
         const Ubyte u1 = (v0 + 2) % 4;
         if (cut.sign(v0) == s) {
            node.push_back(get_vertex(v0));
            node.push_back(get_vertex(cut(0)));
            node.push_back(get_vertex(cut(1)));
            node.push_back(get_vertex(v1));
         } else {
            node.push_back(get_vertex(u0));
            node.push_back(get_vertex(cut(0)));
            node.push_back(get_vertex(cut(1)));
            node.push_back(get_vertex(u1));
         }
      }
   }
}
template <>
void Partition<Hexa>::get_list_node(std::vector<R3> &node, int s) const {
   assert(0);
}
template <>
void Partition<Tet>::get_list_node(std::vector<R3> &node, int s) const {
   assert(0);
}
template <>
void Partition<Triangle3>::get_list_node(std::vector<R3> &node, int s) const {
   assert(0);
}

// Gustaf
template <>
void Partition<Quad3>::get_list_node(std::vector<R3> &node, int s) const {
   node.resize(Quad3::nv + Quad3::ne);

   for (int i = 0; i < Quad3::nv; ++i) {
      node[i] = T[i];
   }

   for (int e = 0; e < Quad3::ne; ++e) {
      const Ubyte v0 = Quad3::nvedge[e][0];
      const Ubyte v1 = Quad3::nvedge[e][1];

      if (ls[v0] == ls[v1]) {
         node[Quad3::nv + e] = T[v0];
      } else {
         const R t = -ls[v0] / (ls[v1] - ls[v0]);
         node[Quad3::nv + e] = (1. - t) * ((R3)T[v0]) + t * ((R3)T[v1]);
      }
   }
}
