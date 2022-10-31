#ifndef _DATA_STRUCT_1D_HPP
#define _DATA_STRUCT_1D_HPP

#include "R1.hpp"
#include "GenericVertex.hpp"
#include "GenericElement.hpp"
#include<array>

typedef double R;
typedef GenericVertex<R1> Vertex1;

struct DataPoint1  {
  static const int NbOfVertices =1;
  static const int NbOfEdges =0;
  static const int NbOfFaces =0;
  static const int NbOfTet =0;
  static const int NbOfAdjElem =1;
  static const int NbOfVertexOnHyperFace =1;
  static const int NbOfRef = 0;
  static const int NbOfVerticesCut = 0;
  static const int ParaviewNumCell = 1;
  static const int nva = 0;
  static const int NvOnFace = 1;
  static const int NbSignPattern = 3;
  static const int NbNtCut = 1;
  static const int NbNtPatch = 1;

  typedef Vertex1 V;
  typedef  V::Rd Rd;
  typedef  R0 Face;
  static R mesure(  V * pv[NbOfVertices]  ) {
    return 1.;
  }
  // static R mesure() {
  //   return 1.;
  // }
  typedef R0 RdHatBord;
  typedef R0 RdHat;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  { return R0() ;}
  // static RdHat PBord()  { return R0() ;}

};
class Node1: public GenericElement<DataPoint1> {
public:
  Node1() {}; // constructor empty for array
  Rd operator()(const RdHat & Phat) const {
    Rd r= (*(Rd*) vertices[0]);
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

struct DataSeg1  {
  static const int NbOfVertices =2;
  static const int NbOfVerticesCut =1;
  static const int NbOfFaces =0;
  static const int NbOfEdges =1;
  static const int NbOfTet =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  static const int ParaviewNumCell = 3;
  static const int NbOfRef = 2;
  static const int NvOnFace = 1;
  static const int NbSignPattern = 9;
  static const int NbNtCut = 1;
  static const int NbNtPatch = 1;

  typedef Vertex1 V;
  typedef  V::Rd Rd ;
  typedef Node1 Face;
  static R mesure(  V *  pv[NbOfVertices]) {
    return pv[1]->x-pv[0]->x;
  }
  typedef R0 RdHatBord;
  typedef R1 RdHat;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  { return R1(*nvb) ;}
  // static RdHat PBord(const int * nvb)  { return R1(*nvb) ;}

};
class Seg1: public GenericElement<DataSeg1>  {
public:

  Seg1() {}; // constructor empty for array


  void Gradlambda(R1 * GradL) const
  {
    GradL[1]= 1./mesure();
    GradL[0]=-GradL[1];
  }

  R1 toKref(const R1& P) const {
    const R &A =*vertices[0];
    const R &B =*vertices[1];

    R mes = fabs(B-A);
    return R1((P.x - A)/mes);
  }
  R1 toReferenceElement(const R1& P) const {
    const R &A =*vertices[0];
    const R &B =*vertices[1];

    R mes = fabs(B-A);
    return R1((P.x - A)/mes);
  }
  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);
    for (int i=1;i<nv;++i)
      r+=  Phat[i-1]*(*(Rd*) vertices[i]);
    return r;
  }

};
class BoundaryPoint1: public GenericElement<DataPoint1>  {
public:
  BoundaryPoint1() {}; // constructor empty for array

  Rd operator()(const RdHat & Phat) const {
    Rd r= (*(Rd*) vertices[0]);
    return r;
  }
};




#endif
