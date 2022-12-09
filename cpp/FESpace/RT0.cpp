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
#include "transformation.hpp"
// P0 2D
class TypeOfFE_RT0_2d : public GTypeOfFE<Mesh2> {
   typedef Mesh2 Mesh;               // Define 2D mesh as Mesh
   typedef typename Mesh::Element E; // Define mesh element as E
   static const int ifGeomOnItem[4]; // Nodes, Edges, Faces, Volume

 public:
   static const int k = 0; // Degree, RT_0
   static const int nbDof =
       (k + 1) * (k + 3); // #dof, and dimension of FE space

   // Initialised below class definition
   static int Data[];
   // static double alpha_Pi_h[];

   // dont understand this constructor aside from it being for 2d
   TypeOfFE_RT0_2d()
       : GTypeOfFE<Mesh2>(
             3,    // #dof on the element
             2,    // dim of basis fun
             Data, // geometry data defined below
             // 1, // #subdivisions for plotting [??]
             // 1, // # sub finite elements (generally 1) [??]
             6, // =kPi=m number of coeff alpha_k
             3  // =npPi=n_p number of geometry locations (points) for dofs
                // alpha_Pi_h // array to store coeff alpha_k
             )  // Inputs to GTypeOfFE
   {
      GTypeOfFE<Mesh>::basisFctType    = BasisFctType::RT0;
      GTypeOfFE<Mesh>::polynomialOrder = 0;

      static const R2 Pt[3] = {R2(1. / 2, 1. / 2), R2(0., 1. / 2),
                               R2(1. / 2, 0)}; // Reference coordinate of dof

      // for loop syntax, its ok to define several ints as below inside the
      // scope of the for loop and only incr one
      for (int l = 0, k = 0; l < nbDof;
           l++) { // Fortin interpolator (defines p_l=coorDof[l], j_k, l_k (or
                  // p_k), and i_k)
         Pt_Pi_h[l] = Pt[l]; // ::: REPLACE WITH ::: coorDof[l] = Pt_Pi_h[l]; //
                             // coordinate of the dofs
         // ipj_Pi_h some kind of interpolation structure, i=dof, p=pointOfDof,
         // j=componentOfBasisFunction Gives for each dof, get what node
         // corresponds to it See GTypeOfFE.hpp for IPJ struct
         for (int j = 0; j < 2; j++) {  // 6 = nbDof*dim entries for ipj_Pi_h
            ipj_Pi_h[k] = IPJ(l, l, j); // defines in order: i_k, l_k, j_k
            k           = k + 1;        // k goes from 0 to 5
         }
      }
   }

   // Defined below class definition
   // freefem: const bool *What_d, const Mesh &Th, const Triangle &K, const R2
   // &PHat, RNMK_ &val
   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const; // [???]
   // const What_d whatd, const Element & K, const R2 & P,RNMK_ & bfMat) const

   // Initialises the coefficients in front of the basis functions
   void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const;
};

const int TypeOfFE_RT0_2d::ifGeomOnItem[4] = {
    0, 1, 0, 0}; // which geometry=vertex,edge,face,volume contains dofs
int TypeOfFE_RT0_2d::Data[]{
    // geometry=vertex,edge,face,volume; 0,1,2 vertices, 3,4,5 edges, 6 face
    3, 4, 5, // on what hard coded geometry the dofs are located
    0, 0, 0, // the (number of the dof - 1) on corresponding geometry/node (here
             // each edge has 1 dof)
    0, 1, 2, // the geometry(vertex) of the basis fcn corresponding to the dof
    0, 1, 2, // the dof of the sub FE (default 0,1,2) [???]
    0, 1, 0, 0, // which geometry=vertex,edge,face,volume contains dofs
    0, // for each component $j=0,N-1$ it give the sub FE associated [???]
    0, // begin_dfcomp
    3  // end_dfcomp (=total number of dof in one element)
};

void TypeOfFE_RT0_2d::FB(const What_d whatd, const Element &K, const R2 &PHat,
                         RNMK_ &bfMat) const {
   assert(bfMat.N() >= 3); // if statement inside not true -- STOP code, in this
                           // case bfMat.N() = #rows of matrix = 3
   bfMat = 0.; // set to zeros Basis function matrix: rows dofs; cols id, dx, dy
   R2 P(K(PHat));
   R2 a(K[0]), b(K[1]), c(K[2]);
   R scaling = 1. / (2 * K.mesure());
   double s  = sqrt(K.measure());
   R const0  = scaling * K.EdgeOrientation(0) * s;
   R const1  = scaling * K.EdgeOrientation(1) * s;
   R const2  = scaling * K.EdgeOrientation(2) * s;

   // whatd = 0,1,2
   if (whatd &
       Fop_D0) { // checks whether whatd = 0, ie function and no derivative ?
      // RN_ baseFuns0(bfMat('.',0,op_id)); // Pointer to (all rows, 1st col) of
      // bfMat RN_ baseFuns1(bfMat('.',1,op_id)); // Pointer to (all rows, 2nd
      // col) of bfMat

      bfMat(0, 0, op_id) =
          const0 * (P.x - a.x); // first component, first basis fun
      bfMat(0, 1, op_id) = const0 * (P.y - a.y); // second comp, first basis fun

      bfMat(1, 0, op_id) = const1 * (P.x - b.x); // first comp, second basis fun
      bfMat(1, 1, op_id) =
          const1 *
          (P.y - b.y); // [???] there is a MINUS sign on these for freefem!!

      bfMat(2, 0, op_id) = const2 * (P.x - c.x);
      bfMat(2, 1, op_id) = const2 * (P.y - c.y);
   }

   if (whatd & Fop_dx) { // here first comp gets differentiated away to 0
      bfMat(0, 0, op_dx) = const0; // first basis fun, first component
      bfMat(1, 0, op_dx) = const1;
      bfMat(2, 0, op_dx) = const2;
   }

   if (whatd & Fop_dy) { // here second comp gets differentiated away to 0
      bfMat(0, 1, op_dy) = const0; // first basis fun, second component
      bfMat(1, 1, op_dy) = const1;
      bfMat(2, 1, op_dy) = const2;
   }

   // some weird assertion, probably not necessary
}

// Initialises v as a vector containing at index k the coeff alpha_k
void TypeOfFE_RT0_2d::get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                                    KN_<double> &v) const {
   const Element &T = K.T;
   double s         = 1. / sqrt(T.measure());
   for (int i = 0, k = 0; i < 3; i++) {
      R2 E(T.Edge(i));
      R sgn = T.EdgeOrientation(i);

      v[k++] = sgn * E.y * s;  // on first run, k=0 and then incremented
      v[k++] = -sgn * E.x * s; // -.-, k=1 and then incremented to 2
   }
}

static TypeOfFE_RT0_2d myRT0_2d;
GTypeOfFE<Mesh2> &RT0_2d(myRT0_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::RT0 = myRT0_2d;

// P0 2D using Piola Contravariant
class TypeOfFE_RT0m_2d : public GTypeOfFE<Mesh2> {
   typedef Mesh2 Mesh;               // Define 2D mesh as Mesh
   typedef typename Mesh::Element E; // Define mesh element as E
   static const int ifGeomOnItem[4]; // Nodes, Edges, Faces, Volume

 public:
   static const int k = 0; // Degree, RT_0
   static const int nbDof =
       (k + 1) * (k + 3); // #dof, and dimension of FE space

   // Initialised below class definition
   static int Data[];
   // static double alpha_Pi_h[];
   PiolaContravariant<E> piola;

   // dont understand this constructor aside from it being for 2d
   TypeOfFE_RT0m_2d()
       : GTypeOfFE<Mesh2>(
             3,    // #dof on the element
             2,    // dim of basis fun
             Data, // geometry data defined below
             // 1, // #subdivisions for plotting [??]
             // 1, // # sub finite elements (generally 1) [??]
             6, // =kPi=m number of coeff alpha_k
             3  // =npPi=n_p number of geometry locations (points) for dofs
                // alpha_Pi_h // array to store coeff alpha_k
             )  // Inputs to GTypeOfFE
   {
      static const R2 Pt[3] = {R2(1. / 2, 1. / 2), R2(0., 1. / 2),
                               R2(1. / 2, 0)}; // Reference coordinate of dof

      // for loop syntax, its ok to define several ints as below inside the
      // scope of the for loop and only incr one
      for (int l = 0, k = 0; l < nbDof;
           l++) { // Fortin interpolator (defines p_l=coorDof[l], j_k, l_k (or
                  // p_k), and i_k)
         Pt_Pi_h[l] = Pt[l]; // ::: REPLACE WITH ::: coorDof[l] = Pt_Pi_h[l]; //
                             // coordinate of the dofs
         // ipj_Pi_h some kind of interpolation structure, i=dof, p=pointOfDof,
         // j=componentOfBasisFunction Gives for each dof, get what node
         // corresponds to it See GTypeOfFE.hpp for IPJ struct
         for (int j = 0; j < 2; j++) {  // 6 = nbDof*dim entries for ipj_Pi_h
            ipj_Pi_h[k] = IPJ(l, l, j); // defines in order: i_k, l_k, j_k
            k           = k + 1;        // k goes from 0 to 5
         }
      }
   }

   // Defined below class definition
   // freefem: const bool *What_d, const Mesh &Th, const Triangle &K, const R2
   // &PHat, RNMK_ &val
   void FB(const What_d, const Element &, const Rd &, RNMK_ &) const; // [???]
   // const What_d whatd, const Element & K, const R2 & P,RNMK_ & bfMat) const

   // Initialises the coefficients in front of the basis functions
   void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const;
};

const int TypeOfFE_RT0m_2d::ifGeomOnItem[4] = {
    0, 1, 0, 0}; // which geometry=vertex,edge,face,volume contains dofs
int TypeOfFE_RT0m_2d::Data[]{
    // geometry=vertex,edge,face,volume; 0,1,2 vertices, 3,4,5 edges, 6 face
    3, 4, 5, // on what hard coded geometry the dofs are located
    0, 0, 0, // the (number of the dof - 1) on corresponding geometry/node (here
             // each edge has 1 dof)
    0, 1, 2, // the geometry(vertex) of the basis fcn corresponding to the dof
    0, 1, 2, // the dof of the sub FE (default 0,1,2) [???]
    0, 1, 0, 0, // which geometry=vertex,edge,face,volume contains dofs
    0, // for each component $j=0,N-1$ it give the sub FE associated [???]
    0, // begin_dfcomp
    3  // end_dfcomp (=total number of dof in one element)
};

void TypeOfFE_RT0m_2d::FB(const What_d whatd, const Element &K, const R2 &PHat,
                          RNMK_ &bfMat) const {
   assert(bfMat.N() >= 3); // if statement inside not true -- STOP code, in this
                           // case bfMat.N() = #rows of matrix = 3
   bfMat = 0.; // set to zeros Basis function matrix: rows dofs; cols id, dx, dy
   R2 P(K(PHat));
   R2 a(K[0]), b(K[1]), c(K[2]);
   R scaling = 1. / (2 * K.mesure());
   double s  = sqrt(K.measure());
   R const0  = scaling * K.EdgeOrientation(0) * s;
   R const1  = scaling * K.EdgeOrientation(1) * s;
   R const2  = scaling * K.EdgeOrientation(2) * s;

   // whatd = 0,1,2
   if (whatd & Fop_D0) {
      bfMat(0, 0, op_id) = -PHat.x; // first component, first basis fun
      bfMat(0, 1, op_id) = -PHat.y; // second comp, first basis fun

      bfMat(1, 0, op_id) = (PHat.x - 1); // first comp, second basis fun
      bfMat(1, 1, op_id) = PHat.y;

      bfMat(2, 0, op_id) = -PHat.x;
      bfMat(2, 1, op_id) = 1 - PHat.y;
   }
   // piola.transform_phi(K, bfMat);

   if (whatd & Fop_dx) {       // here first comp gets differentiated away to 0
      bfMat(0, 0, op_dx) = -1; // first basis fun, first component
      bfMat(1, 0, op_dx) = 1;
      bfMat(2, 0, op_dx) = -1;
   }

   if (whatd & Fop_dy) {       // here second comp gets differentiated away to 0
      bfMat(0, 1, op_dy) = -1; // first basis fun, second component
      bfMat(1, 1, op_dy) = 1;
      bfMat(2, 1, op_dy) = -1;
   }
}

// Initialises v as a vector containing at index k the coeff alpha_k
void TypeOfFE_RT0m_2d::get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                                     KN_<double> &v) const {
   const Element &T = K.T;
   double s         = 1. / sqrt(T.measure());
   for (int i = 0, k = 0; i < 3; i++) {
      R2 E(T.Edge(i));
      R sgn = T.EdgeOrientation(i);

      v[k++] = sgn * E.y * s;  // on first run, k=0 and then incremented
      v[k++] = -sgn * E.x * s; // -.-, k=1 and then incremented to 2
   }
}

static TypeOfFE_RT0m_2d myRT0m_2d;
GTypeOfFE<Mesh2> &RT0m_2d(myRT0m_2d);
template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::RT0m = myRT0m_2d;

class TypeOfFE_RT0_3d : public GTypeOfFE<Mesh3> {
 public:
   typedef Mesh3 Mesh;
   typedef typename Mesh::Element Element;
   // typedef GFElement<Mesh3> FElement;

   static int dfon[];
   static const int d = Mesh::Rd::d;
   TypeOfFE_RT0_3d();
   // int edgeface[4][3] ;
   static int Data[];

   void FB(const What_d whatd, const Element &K, const R3 &PHat,
           RNMK_ &bfMat) const;
   void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const;
};

int TypeOfFE_RT0_3d::dfon[] = {0, 0, 1, 0};
int TypeOfFE_RT0_3d::Data[]{
    // geometry=vertex,edge,face,volume; 0,1,2 vertices, 3,4,5 edges, 6 face
    10, 11, 12,
    13, // on what hard coded geometry the dofs are located
    0, 0, 0,
    0, // the (number of the dof - 1) on corresponding geometry/node (here each
    // edge has 1 dof)
    0, 1, 2,
    3, // the geometry(vertex) of the basis fcn corresponding to the dof
    0, 1, 2,
    3, // the dof of the sub FE (default 0,1,2) [???]
    0, 0, 1,
    0, // which geometry=vertex,edge,face,volume contains dofs
    0, // for each component $j=0,N-1$ it give the sub FE associated [???]
    0, // begin_dfcomp
    4  // end_dfcomp (=total number of dof in one element)
};

TypeOfFE_RT0_3d::TypeOfFE_RT0_3d() : GTypeOfFE<Mesh3>(4, 3, Data, 12, 4) {

   GTypeOfFE<Mesh>::basisFctType    = BasisFctType::RT0;
   GTypeOfFE<Mesh>::polynomialOrder = 0;

   R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};
   for (int i = 0; i < Element::nf; ++i) {
      this->Pt_Pi_h[i] =
          (Pt[Element::nvface[i][0]] + Pt[Element::nvface[i][1]] +
           Pt[Element::nvface[i][2]]) /
          3.;
   }
   int i = 0;
   for (int f = 0; f < 4; f++) {
      for (int c = 0; c < 3; c++, i++) {
         // this->pInterpolation[i]=e;
         // this->cInterpolation[i]=c;
         // this->dofInterpolation[i]=f;
         // this->coefInterpolation[i]=0.;
         ipj_Pi_h[i] = IPJ(f, f, c); // defines in order: i_k, l_k, j_k
      }
   }
}

void TypeOfFE_RT0_3d::get_Coef_Pi_h(const GbaseFElement<Mesh> &K,
                                    KN_<double> &v) const {
   const Element &T = K.T;
   for (int f = 0, k = 0; f < 4; f++) {
      R3 N = T.N_notNormalized(f); //  exterior and  ||N|| = 2* area f
      N *= T.faceOrient(f) / 2.;

      v[k++] = N.x;
      v[k++] = N.y;
      v[k++] = N.z;
   }
}
void TypeOfFE_RT0_3d::FB(const What_d whatd, const Element &K, const R3 &PHat,
                         RNMK_ &bfMat) const {
   assert(bfMat.N() >= 4);
   assert(bfMat.M() == 3);
   // wi = signe * (x - qi)/ (volume*d)
   bfMat   = 0;
   // cout << " TypeOfFE_RT0_3d "<< Th(K) << " " << Th.nt << " /
   // "<<K.faceOrient(0)<< " " << K.faceOrient(1) << " " << K.faceOrient(2) << "
   // " << K.faceOrient(3) <<endl;
   R cc    = 1. / (d * K.mesure());
   R ci[4] = {cc * K.faceOrient(0), cc * K.faceOrient(1), cc * K.faceOrient(2),
              cc * K.faceOrient(3)};

   if (whatd & Fop_D0) {
      R3 X  = K(PHat);
      int k = 0;
      for (int i = 0; i < 4; ++i) {
         R3 wi              = (X - K[i]) * ci[i];
         bfMat(i, 0, op_id) = wi.x;
         bfMat(i, 1, op_id) = wi.y;
         bfMat(i, 2, op_id) = wi.z;
      }
   }

   if (whatd & Fop_D1) {
      RN_ Ci(ci, 4);
      if (whatd & Fop_dx)
         bfMat('.', 0, op_dx) = Ci;
      if (whatd & Fop_dy)
         bfMat('.', 1, op_dy) = Ci;
      if (whatd & Fop_dz)
         bfMat('.', 2, op_dz) = Ci;
   }
}

static TypeOfFE_RT0_3d RT0_3d;
GTypeOfFE<Mesh3> &RT03d(RT0_3d);
template <> GTypeOfFE<Mesh3> &DataFE<Mesh3>::RT0 = RT03d;

//
