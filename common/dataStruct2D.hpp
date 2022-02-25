#ifndef _DATA_STRUCT_2D_HPP
#define _DATA_STRUCT_2D_HPP

#include "R1.hpp"
#include "R2.hpp"
#include "GenericVertex.hpp"
#include "GenericElement.hpp"



typedef double R;
typedef GenericVertex<R2> Vertex2;

struct DataTriangle2  {
  static const int NbOfVertices =3;
  static const int NbOfFaces =1;
  static const int NbOfEdges =3;
  static const int NbOfTet =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  static const int NbOfVerticesCut =2;
  static const int NbOfRef = 4;
  static const int ParaviewNumCell = 5;
  static const int NvOnFace = 3;
  static const int NbSignPattern = 27;
  static const int NbNtCut = 1;
  static const int NbNtPatch = 1;

  typedef Vertex2 V;
  typedef  V::Rd Rd ;

  static R mesure(  V *  pv[NbOfVertices]) {
    return det(*pv[0],*pv[1],*pv[2])*0.5;
  }
  typedef R2 RdHat;
  typedef R1 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  {
    return RdHat::KHat[nvb[0]]*(1-P.x)+R2::KHat[nvb[1]]*(P.x) ;}

};

struct DataSeg2  {
  static const int NbOfVertices =2;
  static const int NbOfEdges =1;
  static const int NbOfFaces =0;
  static const int NbOfTet =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  static const int NbOfRef = 2;
  static const int NbOfVerticesCut = 1;
  static const int nva = 1;
  static const int NvOnFace = 1;
  static const int NbSignPattern = 9;
  static const int NbNtCut = 1;
  static const int NbNtPatch = 1;

  typedef Vertex2 V;
  typedef  V::Rd Rd;
  static R mesure(  V *  pv[NbOfVertices]) {
    return R2(*pv[0],*pv[1]).norme();
  }
  typedef R1 RdHat;
  typedef R0 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord &P)  { return RdHat(*nvb) ;}
  // static RdHat PBord(const int * nvb)  { return RdHat(*nvb) ;}

  //static const int (* const nvface)[3];// = nvfaceSeg ;
  //static const int (* const nvedge)[2];//  = nvedgeSeg;

};

struct DataQuad2  {
  static const int NbOfVertices =4;
  static const int NbOfFaces =1;
  static const int NbOfEdges =4;
  static const int NbOfTet =0;
  static const int NbOfAdjElem = 4;
  static const int NbOfVertexOnHyperFace = 2;
  static const int NbOfVerticesCut =2;
  static const int NvOnFace = 4;
  static const int ParaviewNumCell = 9;
  static const int NbSignPattern = 81;
  static const int NbNtCut = 2;
  static const int NbNtPatch = 2;

  typedef Vertex2 V;
  typedef  V::Rd Rd ;

  static R mesure(  V *  pv[NbOfVertices]) {
    return det(*pv[0],*pv[1],*pv[2]);
  }
  typedef R2 RdHat;
  typedef R1 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  {
    return RdHat::KHat[nvb[0]]*(1-P.x)+R2::KHat[nvb[1]]*(P.x) ;}

};

class Triangle2: public GenericElement<DataTriangle2>
{
public:
  Triangle2() {}; // constructor empty for array
  Triangle2(Vertex * v0,int * iv,int r=0, double mss=UnSetMesure) {
    this->set(v0, iv, r , mss);
  }; // constructor empty for array


  R2 H(int i) const { ASSERTION(i>=0 && i <3);
    R2 E=Edge(i);return E.perp()/(2.*this->mesure());} // heigth

  void Gradlambda(R2 * GradL) const  {
    GradL[1]= H(1);
    GradL[2]= H(2);
    GradL[0]=-GradL[1]-GradL[2];
  }

  R2 toKref(const R2& P) const {
    R l[3];
    const R2 &A =*vertices[0];
    const R2 &B =*vertices[1];
    const R2 &C =*vertices[2];

    R2 PA(P,A), PB(P,B), PC(P,C);
    l[0] = 0.5/mes*( (PB ^ PC));
    l[1] = 0.5/mes*( (PC ^ PA));
    l[2] = 1 - l[0] - l[1];
    return l[0]*R2::KHat[0] + l[1]*R2::KHat[1] + l[2]*R2::KHat[2];
  }

  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);
    for (int i=1;i<nv;++i)
      r+=  Phat[i-1]*(*(Rd*) vertices[i]);
    return r;
  }

  R2 centroid() const {
    return 1./3*((*vertices[0])+(*vertices[1])+(*vertices[2]));
  }

  R2 toKref(const R1& P, int i) const;
  R mesureBord(int i) const ;


};
class BoundaryEdge2: public GenericElement<DataSeg2>
{
public:
  BoundaryEdge2() {}; // constructor empty for array
  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);
    for (int i=1;i<nv;++i)
      r+=  Phat[i-1]*(*(Rd*) vertices[i]);
    return r;
  }
};

class Quad2: public GenericElement<DataQuad2> {
public:
  Quad2() {}; // constructor empty for array
  Quad2(Vertex * v0,int * iv,int r=0, double mss=UnSetMesure) {
    this->set(v0, iv, r , mss);
  }; // constructor empty for array

  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0])
    + Phat[0]*(*(Rd*) vertices[1])
    + Phat[1]*(*(Rd*) vertices[3]);
    return r;
  }

  // R2 H(int i) const { ASSERTION(i>=0 && i <3);
  //   R2 E=Edge(i);return E.perp()/(2.*this->mesure());
  // } // heigth

  // void Gradlambda(R2 * GradL) const
  // {
  //   GradL[1]= H(1);
  //   GradL[2]= H(2);
  //   GradL[0]=-GradL[1]-GradL[2];
  // }

  R2 toKref(const R2& P) const {
    const R2 &A =*vertices[0];
    const R2 &B =*vertices[1];
    const R2 &C =*vertices[3];

    R2 AB(A,B), AC(A,C);
    R2 AP(A,P);
    double c1 = (AP,AB)/AB.norme2();
    double c2 = (AP,AC)/AC.norme2();
    return R2(c1,c2);

  }

  // R2 centroid() const {
  //   return 1./3*((*vertices[0])+(*vertices[1])+(*vertices[2]));
  // }

  // R2 toKref(const R1& P, int i) const;
  // R mesureBord(int i) const ;

};






#endif
