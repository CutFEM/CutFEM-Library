#ifndef COMMON_DATA_STRUCT_3D_HPP
#define COMMON_DATA_STRUCT_3D_HPP

#include "GenericElement.hpp"
#include <cassert>

typedef double R;
typedef GenericVertex<R3> Vertex3;

struct DataPoint3 {
   static const int NbOfVertices          = 1;
   static const int NbOfEdges             = 0;
   static const int NbOfFaces             = 0;
   static const int NbOfTet               = 0;
   static const int NbOfAdjElem           = 1;
   static const int NbOfVertexOnHyperFace = 1;
   static const int ParaviewNumCell       = 3;
   static const int NbOfRef               = 0;
   static const int NbOfVerticesCut       = 0;
   static const int nva                   = 0;
   static const int NvOnFace              = 1;
   static const int NbSignPattern         = 3;
   static const int NbNtCut               = 1;
   static const int NbNtPatch             = 1;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   typedef R0 Face;
   static R mesure(V *pv[NbOfVertices]) { return 1.; }
   // static R mesure() {
   //   return 1.;
   // }
   typedef R0 RdHatBord;
   typedef R0 RdHat;
   static RdHat PBord(const int *nvb, const RdHatBord &P) { return R0(); }
   // static RdHat PBord()  { return R0() ;}
};
class Node3 : public GenericElement<DataPoint3> {
 public:
   Node3(){}; // constructor empty for array
   Rd operator()(const RdHat &Phat) const {
      Rd r = (*(Rd *)vertices[0]);
      return r;
   }
   // std::array<int, 2> index_adjacent_element_;
   // void set_adjacent_element(int k1, int k2) {
   //   index_adjacent_element_ = {k1, k2};
   // }
   // int get_indes_adjacent_element(int i) const {
   //   return index_adjacent_element_[i];
   // }
};

struct DataSeg3 {
   static const int NbOfVertices          = 2;
   static const int NbOfEdges             = 1;
   static const int NbOfFaces             = 0;
   static const int NbOfTet               = 0;
   static const int NbOfAdjElem           = NbOfVertices;
   static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
   static const int NbOfRef               = 2;
   static const int NbOfVerticesCut       = 1;
   static const int ParaviewNumCell       = 3;
   static const int nva                   = 1;
   static const int NvOnFace              = 1;
   static const int NbSignPattern         = 9;
   static const int NbNtCut               = 1;
   static const int NbNtPatch             = 1;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   static R mesure(V *pv[NbOfVertices]) { return R3(*pv[0], *pv[1]).norme(); }
   typedef R1 RdHat;
   typedef R0 RdHatBord;
   typedef Node3 Face;
   static RdHat PBord(const int *nvb, const RdHatBord &P) {
      return RdHat(*nvb);
   }
   // static RdHat PBord(const int * nvb)  { return RdHat(*nvb) ;}

   // static const int (* const nvface)[3];// = nvfaceSeg ;
   // static const int (* const nvedge)[2];//  = nvedgeSeg;
};
class Edge3 : public GenericElement<DataSeg3> {
 public:
   Edge3(){}; // constructor empty for array
   Rd operator()(const RdHat &Phat) const {
      Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
      for (int i = 1; i < nv; ++i)
         r += Phat[i - 1] * (*(Rd *)vertices[i]);
      return r;
   }
   // std::array<int, 2> index_adjacent_element_;
   // void set_adjacent_element(int k1, int k2) {
   //   index_adjacent_element_ = {k1, k2};
   // }
   // int get_index_adjacent_element(int i) const {
   //   return index_adjacent_element_[i];
   // }
};

struct DataTriangle3 {
   static const int NbOfVertices          = 3;
   static const int NbOfVerticesCut       = 3;
   static const int NbOfEdges             = 3;
   static const int NbOfFaces             = 1;
   static const int NbOfTet               = 0;
   static const int NbOfAdjElem           = NbOfVertices;
   static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
   static const int NbOfRef               = 4;
   static const int NvOnFace              = 3;
   static const int ParaviewNumCell       = 5;
   static const int NbSignPattern         = 27;
   static const int NbNtCut               = 1;
   static const int NbNtPatch             = 1;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   typedef Edge3 Face;
   typedef R2 RdHat;
   typedef R1 RdHatBord;
   static RdHat PBord(const int *nvb, const RdHatBord &P) {
      return RdHat::KHat[nvb[0]] * (1 - P.X()) + R2::KHat[nvb[1]] * (P.X());
   }

   static R mesure(V *pv[NbOfVertices]) {
      return (R3(*pv[0], *pv[1]) ^ R3(*pv[0], *pv[2])).norme() * 0.5;
   }
};
class Triangle3 : public GenericElement<DataTriangle3> {
 public:
   typedef Triangle3 TypeCutElement;
   typedef Edge3 Face;
   Triangle3(){}; // constructor empty for array

   R3 Edge(int i) const {
      assert(i >= 0 && i < 3);
      return Rd(this->at((i + 1) % 3), this->at((i + 2) % 3));
   }

   Rd operator()(const RdHat &Phat) const {
      Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
      for (int i = 1; i < nv; ++i)
         r += Phat[i - 1] * (*(Rd *)vertices[i]);
      return r;
   }
   // std::array<int, 2> index_adjacent_element_;
   // void set_adjacent_element(int k1, int k2) {
   //   index_adjacent_element_ = {k1, k2};
   // }
   // int get_index_adjacent_element(int i) const {
   //   return index_adjacent_element_[i];
   // }
};

struct DataQuad3 {
   static const int NbOfVertices          = 4;
   static const int NbOfFaces             = 1;
   static const int NbOfEdges             = 4;
   static const int NbOfTet               = 0;
   static const int NbOfAdjElem           = 4;
   static const int NbOfVertexOnHyperFace = 2;
   static const int NvOnFace              = 4;
   static const int NbOfRef               = 4;
   static const int ParaviewNumCell       = 9;
   static const int NbSignPattern         = 81;
   static const int NbNtCut               = 2;
   static const int NbNtPatch             = 2;
   static const int NbOfVerticesCut       = 2;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   typedef Edge3 Face;

   static R mesure(V *pv[NbOfVertices]) { return det(*pv[0], *pv[1], *pv[2]); }
   typedef R2 RdHat;
   typedef R1 RdHatBord;
   static RdHat PBord(const int *nvb, const RdHatBord &P) {
      return RdHat::KHat[nvb[0]] * (1 - P.X()) + R2::KHat[nvb[1]] * (P.X());
   }
};
class Quad3 : public GenericElement<DataQuad3> {
 public:
   typedef Edge3 Face;
   typedef Triangle3 TypeCutElement;
   Quad3(){}; // constructor empty for array

   Rd operator()(const RdHat &Phat) const {
      Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]) +
             Phat[0] * (*(Rd *)vertices[1]) + Phat[1] * (*(Rd *)vertices[3]);
      return r;
   }
   // std::array<int, 2> index_adjacent_element_;
   // void set_adjacent_element(int k1, int k2) {
   //   index_adjacent_element_ = {k1, k2};
   // }
   // int get_index_adjacent_element(int i) const {
   //   return index_adjacent_element_[i];
   // }
};

struct DataTet {
   static const int NbOfVertices          = 4;
   static const int NbOfEdges             = 6;
   static const int NbOfFaces             = 4;
   static const int NbOfTet               = 1;
   static const int NbOfAdjElem           = NbOfVertices;
   static const int NbOfVertexOnHyperFace = NbOfVertices - 1;
   static const int NbOfRef               = 8;

   static const int ParaviewNumCell = 10;
   static const int NvOnFace        = 3;
   static const int NbSignPattern   = 81;
   static const int NbNtCut         = 1;
   static const int NbNtPatch       = 1;
   static const int NbOfVerticesCut = 4;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   typedef Triangle3 Face;

   static R mesure(V *pv[NbOfVertices]) {
      R3 AB(*pv[0], *pv[1]);
      R3 AC(*pv[0], *pv[2]);
      R3 AD(*pv[0], *pv[3]);
      return det(AB, AC, AD) / 6.;
   }
   static const int (*const nvface)[3]; // = nvfaceTet;
   static const int (*const nvedge)[2]; //  = nvedgeTet;
   typedef R3 RdHat;
   typedef R2 RdHatBord;
   static RdHat PBord(const int *nvb, const RdHatBord &P) {
      return RdHat::KHat[nvb[0]] * (1 - P.x - P.y) +
             RdHat::KHat[nvb[1]] * (P.x) + RdHat::KHat[nvb[2]] * (P.y);
   }
};
class Tet : public GenericElement<DataTet> {
 public:
   typedef Triangle3 Face;
   typedef Tet TypeCutElement;

   static const std::vector<int> oppEdgeOfEdge; //

   Tet(){}; // constructor empty for array

   Rd operator()(const RdHat &Phat) const {
      Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]);
      for (int i = 1; i < nv; ++i)
         r += Phat[i - 1] * (*(Rd *)vertices[i]);
      return r;
   }

   R3 H(int i) const;
   // { assert(i>=0 && i <4);
   //     R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
   //     R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
   //     return AB^AC/(6.*this->mesure());} // heigth

   R3 n(int i) const;
   //    { assert(i>=0 && i <4);
   // 	R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
   // 	R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
   // 	R3 N=AB^AC;
   //     return N/N.norme();} //  exterior normal

   void Gradlambda(R3 *GradL) const {
      R3 V1(at(0), at(1));
      R3 V2(at(0), at(2));
      R3 V3(at(0), at(3));
      R det1   = 1. / (6. * mesure());
      GradL[1] = (V2 ^ V3) * det1;
      GradL[2] = (V3 ^ V1) * det1;
      GradL[3] = (V1 ^ V2) * det1;
      GradL[0] = -GradL[1] - GradL[2] - GradL[3];
   }

   int faceOrient(int i) const;

   R3 toKref(const R3 &P) const {
      R l[4];
      const R3 &A = at(0);
      const R3 &B = at(1);
      const R3 &C = at(2);
      const R3 &D = at(3);
      R3 PA(P, A), PB(P, B), PC(P, C), PD(P, D);
      l[0] = fabs(det(PB, PC, PD) / (6 * mes));
      l[1] = fabs(det(PC, PA, PD) / (6 * mes));
      l[2] = fabs(det(PB, PA, PD) / (6 * mes));
      l[3] = 1 - l[0] - l[1] - l[2];
      return l[0] * R3::KHat[0] + l[1] * R3::KHat[1] + l[2] * R3::KHat[2] +
             l[3] * R3::KHat[3];
   }
   R3 mapToReferenceElement(const R3 &P) const {
      R l[4];
      const R3 &A = at(0);
      const R3 &B = at(1);
      const R3 &C = at(2);
      const R3 &D = at(3);
      R3 PA(P, A), PB(P, B), PC(P, C), PD(P, D);
      l[0] = fabs(det(PB, PC, PD) / (6 * mes));
      l[1] = fabs(det(PC, PA, PD) / (6 * mes));
      l[2] = fabs(det(PB, PA, PD) / (6 * mes));
      l[3] = 1 - l[0] - l[1] - l[2];
      return l[0] * R3::KHat[0] + l[1] * R3::KHat[1] + l[2] * R3::KHat[2] +
             l[3] * R3::KHat[3];
   }

   R3 toKref(const R2 &P, int i) const;
   R3 mapToReferenceElement(const R2 &P, int i) const;
   R mesureBord(int i) const;
};

struct DataHexa {
   static const int NbOfVertices          = 8;
   static const int NbOfEdges             = 12;
   static const int NbOfFaces             = 6;
   static const int NbOfTet               = 1;
   static const int NbOfAdjElem           = 6;
   static const int NbOfVertexOnHyperFace = 4;
   static const int NbOfRef               = 8;
   static const int ParaviewNumCell       = 12;
   static const int NvOnFace              = 4;
   static const int NbSignPattern         = 6561; // 3^8
   static const int NbNtCut         = 5; // nb of tetra for an non cut element
   static const int NbNtPatch       = 4; // maximum 4 triangles in a patch
   static const int NbOfVerticesCut = 6;

   typedef Vertex3 V;
   typedef V::Rd Rd;
   typedef Quad3 Face;

   static R mesure(V *pv[NbOfVertices]) {
      R3 AB(*pv[0], *pv[1]);
      R3 AC(*pv[0], *pv[3]);
      R3 AD(*pv[0], *pv[4]);
      return det(AB, AC, AD);
   }

   typedef R3 RdHat;
   typedef R2 RdHatBord;
   // static RdHat PBord(const int * nvb,const RdHatBord& P)  {
   //   return RdHat::KHat[nvb[0]]*(1-P.x-P.y)+
   //     RdHat::KHat[nvb[1]]*(P.x)+RdHat::KHat[nvb[2]]*(P.y) ;}
};
class Hexa : public GenericElement<DataHexa> {
 public:
   typedef Quad3 Face;
   typedef Tet TypeCutElement;

   static const std::vector<int> oppEdgeOfEdge;
   static const std::vector<std::vector<int>> nodeConnectivity;

   Hexa(){}; // constructor empty for array
   Hexa(Vertex *v0, int *iv, int r = 0,
        double mss = globalVariable::UnSetMesure) {
      this->set(v0, iv, r, mss);
   }; // constructor empty for array

   Rd operator()(const RdHat &Phat) const {
      Rd r = (1. - Phat.sum()) * (*(Rd *)vertices[0]) +
             Phat[0] * (*(Rd *)vertices[1]) + Phat[1] * (*(Rd *)vertices[3]) +
             Phat[2] * (*(Rd *)vertices[4]);
      return r;
   }

   void Gradlambda(R3 *GradL) const {
      assert(0);
      // R3 V1(at(0),at(1));
      // R3 V2(at(0),at(2));
      // R3 V3(at(0),at(3));
      // R det1=1./(6.*mesure());
      // GradL[1]= (V2^V3)*det1;
      // GradL[2]= (V3^V1)*det1;
      // GradL[3]= (V1^V2)*det1;
      // GradL[0]=-GradL[1]-GradL[2]-GradL[3];
   }

   R3 toKref(const R3 &P) const {
      const R3 &A = *vertices[0];
      const R3 &B = *vertices[1];
      const R3 &C = *vertices[3];
      const R3 &D = *vertices[4];

      R3 AB(A, B), AC(A, C), AD(A, D);
      R3 AP(A, P);
      double c1 = (AP, AB) / AB.norme2();
      double c2 = (AP, AC) / AC.norme2();
      double c3 = (AP, AD) / AD.norme2();

      return R3(c1, c2, c3);
   }

   R3 mapToReferenceElement(const R3 &P) const {
      const R3 &A = *vertices[0];
      const R3 &B = *vertices[1];
      const R3 &C = *vertices[3];
      const R3 &D = *vertices[4];

      R3 AB(A, B), AC(A, C), AD(A, D);
      R3 AP(A, P);
      double c1 = (AP, AB) / AB.norme2();
      double c2 = (AP, AC) / AC.norme2();
      double c3 = (AP, AD) / AD.norme2();

      return R3(c1, c2, c3);
   }

   // R2 centroid() const {
   //   return 1./3*((*vertices[0])+(*vertices[1])+(*vertices[2]));
   // }

   // R2 toKref(const R1& P, int i) const;
   // R mesureBord(int i) const ;
};

#endif
