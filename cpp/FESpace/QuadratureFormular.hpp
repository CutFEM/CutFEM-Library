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

#ifndef _QuadratureFormular_h
#define _QuadratureFormular_h

#include <iostream>
#include <cassert>
#include "../common/point.hpp"

struct QuadratureWeight {
   double a;
   QuadratureWeight(double aa) : a(aa) {}
};

template <class Rd>
class GQuadraturePoint : public QuadratureWeight, public Rd {
   typedef double R;

 public:
   R getWeight() const { return this->a; }

   typedef GQuadraturePoint QP;
   GQuadraturePoint() : QuadratureWeight(0), Rd() {}
   GQuadraturePoint(R aa, const Rd &xx) : QuadratureWeight(aa), Rd(xx) {}
   GQuadraturePoint(const Rd &xx, R aa) : QuadratureWeight(aa), Rd(xx) {}
   operator R() const { return a; }
   GQuadraturePoint(R aa, R xx) : QuadratureWeight(aa), Rd(xx) {}
   GQuadraturePoint(R aa, R x, R y) : QuadratureWeight(aa), Rd(x, y) {}
   GQuadraturePoint(R aa, R x, R y, R z) : QuadratureWeight(aa), Rd(x, y, z) {}
};

template <class Rdd> class GQuadratureFormular {
 public:
   typedef Rdd Rd;
   typedef GQuadraturePoint<Rd> QuadraturePoint;
   typedef GQuadraturePoint<Rd> QP;
   typedef double R;
   const int exact; // exact
   const int n;     // nombre de point d'integration
 private:
   QP *p; // les point d'integration
   const bool clean;

 public:
   R getNbrOfQuads() const { return this->n; }

   // -- les fonctions ------------------
   void Verification(){}; // for verification
   GQuadratureFormular(int e, int NbOfNodes, QuadraturePoint *pp,
                       bool c = false)
       : exact(e), n(NbOfNodes), p(pp), clean(c) {
      Verification();
   }
   GQuadratureFormular(int e, int NbOfNodes, const QuadraturePoint *pp,
                       bool c = false)
       : exact(e), n(NbOfNodes), p(pp), clean(c) {
      Verification();
   }
   GQuadratureFormular(int ex, QP p0, QP p1, QP p2, QP p3, QP p4)
       : exact(ex), n(5), p(new QP[5]), clean(true) {
      p[0] = p0;
      p[1] = p1;
      p[2] = p2;
      p[3] = p3;
      p[4] = p4;
      Verification();
   }
   GQuadratureFormular(int ex, QP p0, QP p1, QP p2, QP p3)
       : exact(ex), n(4), p(new QP[4]), clean(true) {
      p[0] = p0, p[1] = p1, p[2] = p2;
      p[3] = p3;
      Verification();
   }
   GQuadratureFormular(int ex, QP p0, QP p1, QP p2)
       : exact(ex), n(3), p(new QP[3]), clean(true) {
      p[0] = p0, p[1] = p1, p[2] = p2;
      Verification();
   }
   GQuadratureFormular(int ex, QP p0, QP p1)
       : exact(ex), n(2), p(new QP[2]), clean(true) {
      p[0] = p0, p[1] = p1;
      Verification();
   }
   GQuadratureFormular(int ex, QP p0)
       : exact(ex), n(1), p(new QP[1]), clean(true) {
      p[0] = p0;
      Verification();
   }

   const QP &operator[](int i) const { return p[i]; }
   const QP &operator()(int i) const { return p[i]; }
   const QP &at(int i) const { return p[i]; }
   ~GQuadratureFormular() {
      if (clean)
         delete[] p;
   }

 private:
   GQuadratureFormular(const GQuadratureFormular &) : exact(0), n(0), p(0) {
      assert(0);
   }
   void operator=(const GQuadratureFormular &) { assert(0); }
   GQuadratureFormular() : exact(0), n(0), p(0) { assert(0); }
   static const GQuadratureFormular *Default;
};

template <class Rd>
std::ostream &operator<<(std::ostream &f, const GQuadraturePoint<Rd> &p) {
   f << '{' << (const double)p << '\t' << (const Rd &)p << '}';
   return f;
}

template <class Rd>
std::ostream &operator<<(std::ostream &f, const GQuadratureFormular<Rd> &fi) {
   f << "nb de point integration " << fi.n << ", adr = " << &f << std::endl;
   for (int i = 0; i < fi.n; i++)
      f << '\t' << fi[i] << std::endl;
   return f;
}

template <class Rd>
std::ostream &operator<<(std::ostream &, const GQuadratureFormular<Rd> &);
template <class Rd>
std::ostream &operator<<(std::ostream &, GQuadraturePoint<Rd> &);

typedef GQuadratureFormular<R1> QuadratureFormular1d;
typedef GQuadratureFormular<R2> QuadratureFormular2d;
typedef GQuadratureFormular<R3> QuadratureFormular3d;

typedef GQuadraturePoint<R1> QuadraturePoint1d;
typedef GQuadraturePoint<R2> QuadraturePoint2d;
typedef GQuadraturePoint<R3> QuadraturePoint3d;

extern const QuadratureFormular1d QF_GaussLegendre1;
extern const QuadratureFormular1d QF_GaussLegendre2;
extern const QuadratureFormular1d QF_GaussLegendre3;
extern const QuadratureFormular1d QF_GaussLegendre4;
extern const QuadratureFormular1d QF_GaussLegendre5;
extern const QuadratureFormular1d QF_LumpP1_1D;
extern const QuadratureFormular1d QF_Lobatto1;
extern const QuadratureFormular1d QF_Lobatto3;
extern const QuadratureFormular1d QF_Lobatto4;
extern const QuadratureFormular1d QF_Lobatto5;
extern const QuadratureFormular1d QF_Lobatto6;
extern const QuadratureFormular1d QF_Lobatto7;
extern const QuadratureFormular1d QF_Lobatto10;
extern const QuadratureFormular1d QF_Euler;

extern const GQuadratureFormular<R2> QF_NODE_TRIANGLE;
extern const GQuadratureFormular<R2> QuadratureFormular_T_1;
extern const GQuadratureFormular<R2> QuadratureFormular_T_1lump;
extern const GQuadratureFormular<R2> QuadratureFormular_T_2;
extern const GQuadratureFormular<R2> QuadratureFormular_T_5;
extern const GQuadratureFormular<R2> QuadratureFormular_T_2_4P1;
extern const GQuadratureFormular<R2> QuadratureFormular_T_7;
extern const GQuadratureFormular<R2> QuadratureFormular_T_9;

extern const GQuadratureFormular<R3> QuadratureFormular_Tet_1;
extern const GQuadratureFormular<R3> QuadratureFormular_Tet_1lump;
extern const GQuadratureFormular<R3> QuadratureFormular_Tet_2;
extern const GQuadratureFormular<R3> QuadratureFormular_Tet_5;
extern const GQuadratureFormular<R3> QuadratureFormular_Tet_7;

template <class Rd> const GQuadratureFormular<Rd> *QF_Simplex(int exact);
const GQuadratureFormular<R2> *QF_Quad(int exact);
const GQuadratureFormular<R3> *QF_Hexa(int exact);

GQuadraturePoint<R1> *GaussLegendre(int nn);
GQuadraturePoint<R2> *GaussLegendre2D(int nn);
GQuadraturePoint<R3> *GaussLegendre3D(int nn);

const GQuadratureFormular<R1> *Lobatto(int nn);
int exactLobatto_nPt(int n);

#endif
