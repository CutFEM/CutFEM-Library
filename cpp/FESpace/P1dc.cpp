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

#include "FESpace.hpp"

// P1 discontinuous 2D
class TypeOfFE_P1dcLagrange2d : public GTypeOfFE<Mesh2> {

   typedef Mesh2 Mesh;
   typedef typename Mesh::Element E;
   static const int nbNodeOnItem[4];

 public:
   static const int k   = 1;
   static const int ndf = (k + 2) * (k + 1) / 2;
   static int Data[];
   static double alpha_Pi_h[];

   TypeOfFE_P1dcLagrange2d() : GTypeOfFE<Mesh2>(3, 1, Data, 3, 3, alpha_Pi_h) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P1dc;
      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R2 Pt[3] = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};

      for (int i = 0; i < ndf; ++i) {
         Pt_Pi_h[i]  = Pt[i];
         ipj_Pi_h[i] = IPJ(i, i, 0);
      }
   }

   // void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
   //   for (int i = 0; i < 3; ++i) v[i] = 1;
   // }

   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};
const int TypeOfFE_P1dcLagrange2d::nbNodeOnItem[4] = {1, 0, 0, 0};
int TypeOfFE_P1dcLagrange2d::Data[]                = {
    6, 6, 6,    // we use the face because we want discontinuous element
    0, 1, 2,    // the number of the df on  the node
    0, 0, 0,    // the node of the df
    0, 1, 2,    // which are de df on sub FE
    0, 0, 1, 0, // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, // begin_dfcomp
    3  // end_dfcomp
};
double TypeOfFE_P1dcLagrange2d::alpha_Pi_h[] = {1., 1., 1.};
void TypeOfFE_P1dcLagrange2d::FB(const What_d whatd, const Element &K,
                                 const R2 &P, RNMK_ &val) const {
   R l[] = {1. - P.sum(), P.x, P.y};

   assert(val.N() >= Element::nv);
   assert(val.M() == 1);

   val = 0;
   RN_ f0(val('.', 0, op_id));

   if (whatd & Fop_D0) {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
   }

   if (whatd & Fop_D1) {
      R2 Dl[3];
      K.Gradlambda(Dl);
      if (whatd & Fop_dx) {
         RN_ f0x(val('.', 0, op_dx));
         f0x[0] = Dl[0].x;
         f0x[1] = Dl[1].x;
         f0x[2] = Dl[2].x;
      }

      if (whatd & Fop_dy) {
         RN_ f0y(val('.', 0, op_dy));
         f0y[0] = Dl[0].y;
         f0y[1] = Dl[1].y;
         f0y[2] = Dl[2].y;
      }
   }
}
static TypeOfFE_P1dcLagrange2d P1dc_2d;
GTypeOfFE<Mesh2> &P1dcLagrange2d(P1dc_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P1dc = P1dc_2d;

//---------------------------------------------------------

/**
 * @brief P1dc basis functions for 2D quadrilateral elements
 *
 */
class TypeOfFE_P1QdcLagrange2d : public GTypeOfFE<MeshQuad2> {

    typedef MeshQuad2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 1;                 /// P1
    static const int ndf = (k + 1) * (k + 1); /// number of degrees of freedom
    static int Data[];
    static double alpha_Pi_h[];

    TypeOfFE_P1QdcLagrange2d() : GTypeOfFE<MeshQuad2>(4, 1, Data, 4, 4, alpha_Pi_h) {

        static const R2 Pt[4] = {R2(0., 0.), R2(1., 0.), R2(1., 1.), R2(0., 1.)}; /// vertices of the reference element

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i]  = Pt[i];
            ipj_Pi_h[i] = IPJ(i, i, 0);
        }
    }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P1QdcLagrange2d::nbNodeOnItem[4] = {1, 0, 0, 0};

int TypeOfFE_P1QdcLagrange2d::Data[] = {
    8, 8, 8, 8, // the support number  of the node of the df
    0, 1, 2, 3, // the number of the df on the node
    0, 0, 0, 0, // the node of the df
    0, 1, 2, 3, // which are de df on sub FE
    0, 0, 1, 0, // nb node on what
    0,          // for each component $j=0,N-1$ it give the sub FE associated
    0,          // begin_dfcomp
    4           // end_dfcomp
};

double TypeOfFE_P1QdcLagrange2d::alpha_Pi_h[] = {1., 1., 1., 1.};

void TypeOfFE_P1QdcLagrange2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
    R lx[] = {1. - P.x, P.x};
    R ly[] = {1. - P.y, P.y};

    R2 edge_x(K[0], K[1]); // [hx, 0]
    R2 edge_y(K[0], K[3]); // [0, hy]

    R x0 = K[0][0], y0 = K[0][1];
    R x1 = K[2][0], y1 = K[2][1];

    // P.x, P.y are the coordinates in the reference element
    // If we want to express in terms of the basis functions
    // on the physical element, we must map these points to
    // the physical element
    R2 p      = K(R2(P.x, P.y));
    R phi_x[] = {x1 - p[0], p[0] - x0};
    R phi_y[] = {y1 - p[1], p[1] - y0};

    // std::cout << "(x0, P.x, x1) = (" << x0 << ", " << P.x << ", " << x1 << ")" << std::endl;
    // std::cout << "(y0, P.y, y1) = (" << y0 << ", " << P.y << ", " << y1 << ")" << std::endl;
    // std::cout << "K[0] = " << K[0] << "\n";
    // std::cout << "K[1] = " << K[1] << "\n";
    // std::cout << "K[2] = " << K[2] << "\n";
    // std::cout << "K[3] = " << K[3] << "\n";
    // std::cout << "edge_x(K[0], K[1]) = " << edge_x << "\n";
    // std::cout << "edge_y(K[0], K[3]) = " << edge_y << "\n";
    // getchar();
    // Linear_Transformation<Element> map(K);

    assert(val.N() >= Element::nv);
    assert(val.M() == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    // Basis functions
    // Note: Since the range of the basis functions coincide on the reference element
    // and the physical element, it suffices to give the functions on the reference element here.
    // I.e.: we only care about evaluating the function in nodes, and on the points, the values of
    // the basis function coincides

    if (whatd & Fop_D0) {
        f0[0] = lx[0] * ly[0]; // On physical element: (x1-x)*(y1-y)/(hx*hy)
        f0[1] = lx[1] * ly[0]; // On physical element: (x-x0)*(y1-y)/(hx*hy)
        f0[2] = lx[1] * ly[1]; // On physical element: (x-x0)*(y-y0)/(hx*hy)
        f0[3] = lx[0] * ly[1]; // On physical element: (x1-x)*(y-y0)/(hx*hyt)
    }

    // Derivatives of basis functions
    // Note: here we need to care about the scaling, since we are in fact taking derivatives of
    // the basis functions on the physical element, not the reference element

    if (whatd & Fop_D1) {

        RN_ f0x(val('.', 0, op_dx)); // derivatives in x direction
        RN_ f0y(val('.', 0, op_dy)); // derivatives in y direction
        R Dlx[] = {-1. / Norme2(edge_x), 1. / Norme2(edge_x)};
        R Dly[] = {-1. / Norme2(edge_y), 1. / Norme2(edge_y)};

        R hx = Norme2(edge_x), hy = Norme2(edge_y);

        // yhat = (y-y0)/hy, 1-yhat = (y1 - y)/hy where y = p[1], x=p[0], yhat = P.y, xhat = P.x

        // Thomas'
        f0x[0] = Dlx[0] * ly[0]; // d(phi_0)/dx = -1*(1-y)/hx
        f0x[1] = Dlx[1] * ly[0]; // d(phi_1)/dx = 1*(1-y)/hx
        f0x[2] = Dlx[1] * ly[1]; // ...
        f0x[3] = Dlx[0] * ly[1];
        f0y[0] = Dly[0] * lx[0];
        f0y[1] = Dly[0] * lx[1];
        f0y[2] = Dly[1] * lx[1];
        f0y[3] = Dly[1] * lx[0];

        // Sebastian's
        //std::cout << "phi_x[0]/hx = " << phi_x[0] / hx << std::endl;

        // // Assert that the basis functions on the physical element take values between 0 and 1
        // assert((0. - 1e-12 <= phi_x[0] / hx) && (phi_x[0] / hx <= 1. + 1e-12));
        // assert((0. - 1e-12 <= phi_x[1] / hx) && (phi_x[1] / hx <= 1. + 1e-12));
        // assert((0. - 1e-12 <= phi_y[0] / hy) && (phi_y[0] / hy <= 1. + 1e-12));
        // assert((0. - 1e-12 <= phi_y[1] / hy) && (phi_y[1] / hy <= 1. + 1e-12));

        // // Assert that the basis functions on the reference element take values between 0 and 1
        // assert((0. - 1e-12 <= ly[0]) && (ly[0] <= 1. + 1e-12));
        // assert((0. - 1e-12 <= ly[1]) && (ly[1] <= 1. + 1e-12));
        // assert((0. - 1e-12 <= lx[0]) && (lx[0] <= 1. + 1e-12));
        // assert((0. - 1e-12 <= lx[1]) && (lx[1] <= 1. + 1e-12));

        // // Assert that points are on the reference element
        // assert((0. - 1e-12 <= P.x) && (P.x <= 1. + 1e-12));
        // assert((0. - 1e-12 <= P.y) && (P.y <= 1. + 1e-12));

        // Below works as well, but is written entirely in terms of the physical basis functions
        // // dx
        // f0x[0] = -1.0*phi_y[0]/(hx*hy);			// d(phi_0)/dx = -1*(y1-y)/(hx*hy)
        // f0x[1] = phi_y[0]/(hx*hy);             // d(phi_1)/dx = 1*(y1-y)/(hx*hy)
        // f0x[2] = phi_y[1]/(hx*hy);             // d(phi_2)/dx = 1*(y-y0)/(hx*hy)
        // f0x[3] = -1.0*phi_y[1]/(hx*hy);        // d(phi_3)/dx = -1*(y-y0)/(hx*hy)

        // // dy
        // f0y[0] = -1.0*phi_x[0]/(hx*hy);        // d(phi_0)/dy = -1*(x1-x)/(hx*hy)
        // f0y[1] = -1.0*phi_x[1]/(hx*hy);        // d(phi_1)/dy = -1*(x-x0)/(hx*hy)
        // f0y[2] = phi_x[1]/(hx*hy);             // d(phi_2)/dy = (x-x0)/(hx*hy)
        // f0y[3] = phi_x[0]/(hx*hy);             // d(phi_3)/dy = (x1-x)/(hx*hy)
    }
}

static TypeOfFE_P1QdcLagrange2d P1Qdc_2d;
GTypeOfFE<MeshQuad2> &P1QdcLagrange2d(P1Qdc_2d);
template <> GTypeOfFE<MeshQuad2> &DataFE<MeshQuad2>::P1dc = P1Qdc_2d;

//---------------------------------------------------------
//---------------------------------------------------------

// P1 discontinuous scaled 2D
class TypeOfFE_P1dcScLagrange2d : public GTypeOfFE<Mesh2> {

   typedef Mesh2 Mesh;
   typedef typename Mesh::Element E;
   static const int nbNodeOnItem[4];

 public:
   static const int k   = 1;
   static const int ndf = (k + 2) * (k + 1) / 2;
   static int Data[];

   TypeOfFE_P1dcScLagrange2d() : GTypeOfFE<Mesh2>(3, 1, Data, 3, 3) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P1dc;
      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R2 Pt[3] = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};

      for (int i = 0; i < ndf; ++i) {
         Pt_Pi_h[i]  = Pt[i];
         ipj_Pi_h[i] = IPJ(i, i, 0);
      }
   }

   // void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
   //   for (int i = 0; i < 3; ++i) v[i] = 1;
   // }
   void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const {
      const Element &T = K.T;
      double s         = T.measure();
      v[0]             = s;
      v[1]             = s;
      v[2]             = s;
   }

   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};
const int TypeOfFE_P1dcScLagrange2d::nbNodeOnItem[4] = {1, 0, 0, 0};
int TypeOfFE_P1dcScLagrange2d::Data[]                = {
    6, 6, 6, // we use the face because we want discontinuous element
    0, 1, 2,    // the number of the df on  the node
    0, 0, 0,    // the node of the df
    0, 1, 2,    // which are de df on sub FE
    0, 0, 1, 0, // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, // begin_dfcomp
    3  // end_dfcomp
};
void TypeOfFE_P1dcScLagrange2d::FB(const What_d whatd, const Element &K,
                                   const R2 &P, RNMK_ &val) const {
   R l[] = {1. - P.sum(), P.x, P.y};

   assert(val.N() >= Element::nv);
   assert(val.M() == 1);
   double scaling = 1. / K.measure();
   val            = 0;
   RN_ f0(val('.', 0, op_id));

   if (whatd & Fop_D0) {
      f0[0] = scaling * l[0];
      f0[1] = scaling * l[1];
      f0[2] = scaling * l[2];
   }

   if (whatd & Fop_D1) {
      R2 Dl[3];
      K.Gradlambda(Dl);
      if (whatd & Fop_dx) {
         RN_ f0x(val('.', 0, op_dx));
         f0x[0] = scaling * Dl[0].x;
         f0x[1] = scaling * Dl[1].x;
         f0x[2] = scaling * Dl[2].x;
      }

      if (whatd & Fop_dy) {
         RN_ f0y(val('.', 0, op_dy));
         f0y[0] = scaling * Dl[0].y;
         f0y[1] = scaling * Dl[1].y;
         f0y[2] = scaling * Dl[2].y;
      }
   }
}

static TypeOfFE_P1dcScLagrange2d P1dcsc_2d;
GTypeOfFE<Mesh2> &P1dcScLagrange2d(P1dcsc_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P1dcsc = P1dcsc_2d;

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
// P1 3D
class TypeOfFE_P1dcLagrange3d : public GTypeOfFE<Mesh3> {

   typedef Mesh3 Mesh;
   typedef typename Mesh::Element E;
   static const int nbNodeOnItem[4];

 public:
   static const int k   = 1;
   static const int ndf = 4;
   static int Data[];
   static double alpha_Pi_h[];

   TypeOfFE_P1dcLagrange3d() : GTypeOfFE<Mesh3>(4, 1, Data, 4, 4, alpha_Pi_h) {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::P1dc;
      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R3 Pt[4] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
                               R3(0., 0., 1.)};

      for (int i = 0; i < ndf; ++i) {
         Pt_Pi_h[i]  = Pt[i];
         ipj_Pi_h[i] = IPJ(i, i, 0);
      }
   }

   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};
const int TypeOfFE_P1dcLagrange3d::nbNodeOnItem[4] = {1, 0, 0, 0};
int TypeOfFE_P1dcLagrange3d::Data[]                = {
    14, 14, 14, 14, // we use the volume because we want discontinuous element
    0,  1,  2,  3, // the number of the df on  the node
    0,  0,  0,  0, // the node of the df
    0,  1,  2,  3, // which are de df on sub FE
    0,  0,  0,  1, // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, // begin_dfcomp
    4  // end_dfcomp
};
double TypeOfFE_P1dcLagrange3d::alpha_Pi_h[] = {1., 1., 1., 1.};
void TypeOfFE_P1dcLagrange3d::FB(const What_d whatd, const Element &K,
                                 const R3 &P, RNMK_ &val) const {
   R l[] = {1. - P.sum(), P.x, P.y, P.z};

   assert(val.N() >= Element::nv);
   assert(val.M() == 1);

   val = 0;
   RN_ f0(val('.', 0, op_id));
   if (whatd & Fop_D0) {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
      f0[3] = l[3];
   }
   if (whatd & Fop_D1) {
      R3 Dl[4];
      K.Gradlambda(Dl);

      if (whatd & Fop_dx) {
         RN_ f0x(val('.', 0, op_dx));
         f0x[0] = Dl[0].x;
         f0x[1] = Dl[1].x;
         f0x[2] = Dl[2].x;
         f0x[3] = Dl[3].x;
      }

      if (whatd & Fop_dy) {
         RN_ f0y(val('.', 0, op_dy));
         f0y[0] = Dl[0].y;
         f0y[1] = Dl[1].y;
         f0y[2] = Dl[2].y;
         f0y[3] = Dl[3].y;
      }

      if (whatd & Fop_dz) {
         RN_ f0z(val('.', 0, op_dz));
         f0z[0] = Dl[0].z;
         f0z[1] = Dl[1].z;
         f0z[2] = Dl[2].z;
         f0z[3] = Dl[3].z;
      }
   }
}

static TypeOfFE_P1dcLagrange3d P1dc_3d;
GTypeOfFE<Mesh3> &P1dcLagrange3d(P1dc_3d);
template <> GTypeOfFE<Mesh3> &DataFE<Mesh3>::P1dc = P1dc_3d;

// P1
class TypeOfFE_P1dcTaylor2d : public GTypeOfFE<Mesh2> {

   typedef Mesh2 Mesh;
   typedef typename Mesh::Element E;
   static const int nbNodeOnItem[4];

 public:
   static const int k   = 1;
   static const int ndf = (k + 2) * (k + 1) / 2;
   static int Data[];
   static double alpha_Pi_h[];
   static const double h_eps;
   TypeOfFE_P1dcTaylor2d() : GTypeOfFE<Mesh2>(3, 1, Data, 5, 5, alpha_Pi_h) {

      GTypeOfFE<Mesh>::polynomialOrder = k;

      static const R2 Pt[5] = {
          R2(1. / 3., 1. / 3.), R2(1. / 3. + h_eps, 1. / 3.),
          R2(1. / 3. - h_eps, 1. / 3.), R2(1. / 3., 1. / 3. + h_eps),
          R2(1. / 3., 1. / 3. - h_eps)};

      // u_h(x_G)
      Pt_Pi_h[0]  = Pt[0];
      ipj_Pi_h[0] = IPJ(0, 0, 0);
      int ii      = 1;
      for (int i = 1; i < 3; ++i) {
         Pt_Pi_h[ii]  = Pt[ii];
         ipj_Pi_h[ii] = IPJ(i, ii, 0);
         ii++;
         Pt_Pi_h[ii]  = Pt[ii];
         ipj_Pi_h[ii] = IPJ(i, ii, 0);
         ii++;
      }
   }

   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

const int TypeOfFE_P1dcTaylor2d::nbNodeOnItem[4] = {1, 0, 0, 0};
const double TypeOfFE_P1dcTaylor2d::h_eps        = 1e-10;
// tre trick is to use the face because we want discontinuous element
// the the node on the face will have 3 df, one for each ndfonVertex
// like P3 on edge
int TypeOfFE_P1dcTaylor2d::Data[]                = {
    6, 6, 6, // we use the face because we want discontinuous element
    0, 1, 2,    // the number of the df on  the node
    0, 0, 0,    // the node of the df
    0, 1, 2,    // which are de df on sub FE
    0, 0, 1, 0, // nb node on what
    0, // for each compontant $j=0,N-1$ it give the sub FE associated
    0, // begin_dfcomp
    3  // end_dfcomp
};

double TypeOfFE_P1dcTaylor2d::alpha_Pi_h[] = {
    1., 1. / (2. * h_eps), -1. / (2. * h_eps), 1. / (2. * h_eps),
    -1. / (2. * h_eps)};

void TypeOfFE_P1dcTaylor2d::FB(const What_d whatd, const Element &K,
                               const R2 &Phat, RNMK_ &val) const {

   R2 X = K(Phat);
   R2 G = K.centroid();
   // R meas = K.mesure();

   assert(val.N() >= Element::nv);
   assert(val.M() == 1);

   val = 0;
   RN_ f0(val('.', 0, op_id));

   if (whatd & Fop_D0) {
      f0[0] = 1.;
      f0[1] = (Phat.x - 1. / 3);
      f0[2] = (Phat.y - 1. / 3);
   }

   if (whatd & Fop_D1) {
      if (whatd & Fop_dx) {
         RN_ f0x(val('.', 0, op_dx));
         f0x[0] = 0.;
         f0x[1] = 1.;
         f0x[2] = 0.;
      }

      if (whatd & Fop_dy) {
         RN_ f0y(val('.', 0, op_dy));
         f0y[0] = 0.;
         f0y[1] = 0.;
         f0y[2] = 1.;
      }
   }
}

static TypeOfFE_P1dcTaylor2d P1dcTaylor_2d;
GTypeOfFE<Mesh2> &P1dcTaylor2d(P1dcTaylor_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::P1dcTaylor = P1dcTaylor_2d;
