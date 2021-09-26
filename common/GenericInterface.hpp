#ifndef GENERICINTERFACE_HPP_
#define GENERICINTERFACE_HPP_
#include <iostream>
#include <cassert>
#include <bitset>
#include "HashTable.hpp"
#include "RNM.hpp"
#include "Label.hpp"
#include "../util/util.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"




#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

struct CutData2{
  typedef typename Mesh2::Element Element;

  bool cut = false;
  R2 vertices[2] = {R2(), R2()};
  byte sign[3];
  int edges[3] = {-1, -1, -1};

  R2 pointFromEdge(int i) const { assert(edges[i] != -1); return vertices[edges[i]];}
  CutData2(){assert(0);}
  CutData2(const Element& T, double* s) {
    for(int i=0;i<3;++i) sign[i] = fsign(s[i]);
    if(s[0] != s[1] || s[1] != s[2]) {
      cut = true;
      int idv = 0;
      for(int i=0; i<3;++i) {
        const Ubyte v0 = Element::nvedge[i][0], v1= Element::nvedge[i][1];
        if(sign[v0] != sign[v1]) {
          assert(idv < 2);
          edges[i] = idv;
          const R t = -s[v0]/(s[v1]-s[v0]);
          vertices[idv] = (1-t) * ((R2) T[v0]) + t * ((R2) T[v1]) ;
          idv++;
        }
      }
    }
  }
  CutData2(byte s[3]) {
    for(int i=0;i<3;++i) sign[i] = s[i];
  }
  CutData2(byte s[3], R2 point[2], int arr_edge[3]) : CutData2(s){
    if(s[0] != s[1] || s[1] != s[2]) {
      cut = true;
      vertices[0] = point[0];
      vertices[1] = point[1];
      edges[0] = arr_edge[0];
      edges[1] = arr_edge[1];
      edges[2] = arr_edge[2];
    }
  }
  bool edge_is_cut(int e) const {return (edges[e] != -1);}
  bool edge_isnot_cut(int e) const {return (edges[e] == -1);}


};

struct CutData3{
  typedef typename Mesh3::Element Element;

  bool cut = false;
  R3 vertices[4] = {R3(), R3(), R3(), R3()};
  byte sign[4];
  int edges[6] = {-1, -1, -1, -1, -1, -1};

  R3 pointFromEdge(int i) const { assert(edges[i] != -1); return vertices[edges[i]];}

  // CutData(const Element& T, double* s) {
  //   for(int i=0;i<3;++i) sign[i] = fsign(s[i]);
  //   if(s[0] != s[1] || s[1] != s[2]) {
  //     cut = true;
  //     int idv = 0;
  //     for(int i=0; i<3;++i) {
  //       const Ubyte v0 = Element::nvedge[i][0], v1= Element::nvedge[i][1];
  //       if(sign[v0] != sign[v1]) {
  //         assert(idv < 2);
  //         edges[i] = idv;
  //         const R t = -s[v0]/(s[v1]-s[v0]);
  //         vertices[idv] = (1-t) * ((R2) T[v0]) + t * ((R2) T[v1]) ;
  //         idv++;
  //       }
  //     }
  //   }
  // }
  CutData3(byte s[4]) {
    for(int i=0;i<4;++i) sign[i] = s[i];
  }
  CutData3(byte s[4], R3 point[3], int arr_edge[6]) : CutData3(s){
    if(s[0] != s[1] || s[1] != s[2] || s[0] != s[3]) {
      cut = true;
      vertices[0] = point[0];
      vertices[1] = point[1];
      vertices[2] = point[2];
      vertices[3] = point[3];

      edges[0] = arr_edge[0];
      edges[1] = arr_edge[1];
      edges[2] = arr_edge[2];
      edges[3] = arr_edge[3];
      edges[4] = arr_edge[4];
      edges[5] = arr_edge[5];
    }
  }

  bool edge_is_cut(int e) const { assert(0);return (edges[e] != -1);}
  bool edge_isnot_cut(int e) const {assert(0);return (edges[e] == -1);}

};

template<int d> struct TypeCutData {typedef CutData2 CutData;};
template<> struct TypeCutData<2> {typedef CutData2 CutData;};
template<> struct TypeCutData<3> {typedef CutData3 CutData;};


template<int N>
struct FaceInterface {

};

template<>
struct FaceInterface<2> : public  SortArray<Uint, 2>,
		    public Label {
  typedef SortArray<Uint, 2> FaceIdx;

  FaceInterface(const Uint& a0,const Uint &a1, int l=0)
  : FaceIdx(a0,a1), Label(l) {}
  FaceInterface(Uint *a, int l =0)   : FaceIdx(a), Label(l) {}
  FaceInterface() : FaceIdx(), Label(0) {}

};
template<>
struct FaceInterface<3> : public  SortArray<Uint, 3>,
		    public Label {
  typedef SortArray<Uint, 3> FaceIdx;


  FaceInterface(const Uint& a0,const Uint &a1,const Uint &a2, int l=0)
  : FaceIdx(a0,a1,a2), Label(l) {}
  FaceInterface(Uint *a, int l =0)   : FaceIdx(a), Label(l) {}
  FaceInterface() : FaceIdx(), Label(0) {}

};


template<typename M>
class GenericInterface  {

public :
  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  static const int nve = Mesh::nva;     // number of nodes on hyperface
  typedef typename Mesh::RefPatch RefPatch;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  // typedef SortArray<Uint, nve> FaceIdx;
  // typedef SortArray<Uint, nve> Face;
  typedef FaceInterface<nve> FaceIdx;
  typedef FaceInterface<nve> Face;

  typedef const FaceIdx* const_face_iterator;
  typedef       FaceIdx*       face_iterator;
  typedef const Rd* const_vertex_iterator;
  typedef       Rd*       vertex_iterator;

  typedef std::pair<Uint, Uint> RemumberVertexPairT; ///< Helper type to handle zero-vertexes


  //protected :
  KN<byte> ls_sign;
  const Mesh* backMesh;

  std::vector<FaceIdx>  faces_;            // All triangles of the interface.
  std::vector<Rd> vertices_;
  std::vector<Uint> element_of_face_;
  std::vector<Rd> outward_normal_;
  std::vector<bool> is_boundary_face_;     // True, if triangle face of one of the tetras
  std::map<int, int> face_of_element_;
private :
  const FaceIdx  make_face (const typename RefPatch::FaceIdx& ref_tri,
			    const typename Mesh::Element& K,
			    const double lset[Element::nv],
			    std::vector<RemumberVertexPairT>& renumber_zero_verts);
  const FaceIdx  make_face (int k,
			    const KN<R>& nodex,
			    const KN<R>& nodey);
  // const FaceIdx  make_face (int k,
  // 			    const Marker& marker);

  Rd make_normal(const typename Mesh::Element& K, const double lset[Element::nv]);
  // Rd make_normal_noGrad(const typename Mesh::Element& K, const double lset[Element::nv]);

public :
  GenericInterface() : backMesh(nullptr) { }
  GenericInterface(const Mesh & MM) : backMesh(&MM) { }
  GenericInterface(const Mesh & MM, const KN<double>& ls) : backMesh(&MM){
		make_patch(ls);
	}

  void make_patch (const Mesh & MM, const KN<double>& ls);
  void make_patch (const KN<double>& ls);
  // void make_patch (const KN<double>& ls, const KN<int>& cut_element, const KN<R>& nodex, const KN<R>& nodey);
  // void make_patch (const Marker& marker);


  Rd operator()(const int k, const int i) const {return vertices_[faces_[k][i]];}
  const Rd& operator()(const int i) const {return vertices_[CheckV(i)];}
  const FaceIdx& operator[](const int k) const {return faces_[CheckT(k)];}
  const Face& getFace(const int k) const {return faces_[CheckT(k)];}

  int CheckV(int i) const { ASSERTION(i>=0 && i < vertices_.size()); return i;}
  int CheckT(int i) const { ASSERTION(i>=0 && i < faces_.size()); return i;}


  Uint idxElementOfFace(const int k) const { return element_of_face_[k];}
  Uint nbElement  () const { return faces_.size(); }
  Rd normal(const int k) const { return outward_normal_[k];}
  bool isCut(const int k) const {
    return (face_of_element_.find(k) != face_of_element_.end());
  }
  // Uint nbVertices () const { return vertices_.size(); }

  virtual Rd mapToFace(const FaceIdx& f, const typename Element::RdHatBord x ) const = 0;
  virtual Rd computeDx(const FaceIdx& f) const = 0;
  virtual CutData getCutData(const int k) const = 0;
  const_face_iterator face_begin () const { return (faces_.begin()).base(); }
  const_face_iterator face_end   () const { return (faces_.end()).base(); }

  const_vertex_iterator vertex_begin () const { return (vertices_.begin()).base(); }
  const_vertex_iterator vertex_end   () const { return (vertices_.end()).base(); }

  #ifdef USE_MPI
  // int first_element() const { return MPIcf::my_rank();}
  // int next_element() const {return MPIcf::size();}
  // int last_element() const { return faces_.size();}
  virtual int first_element() const { return MPIcf::first_element(faces_.size());}
  virtual int next_element() const {  return MPIcf::next_element(faces_.size());}
  virtual int last_element() const {  return MPIcf::last_element(faces_.size());}

  #else
  int first_element() const { return 0;}
  int next_element() const {return 1;}
  int last_element() const { return faces_.size();}
  #endif

	double distance(Rd P) const;

	~GenericInterface() {}
private:
  GenericInterface(const GenericInterface &); // pas de construction par copie
  void operator=(const GenericInterface &);// pas affectation par copy


};


template<typename M>
void GenericInterface<M>::make_patch(const Mesh & MM, const KN<double>& ls) {
  backMesh = &MM;
  make_patch(ls);
}

template<typename M>
void GenericInterface<M>::make_patch(const KN<double>& ls) {

  assert(backMesh);
  faces_.resize( 0);                          // reinitialize arrays
  is_boundary_face_.resize( 0);
  vertices_.resize(0);
  element_of_face_.resize(0);
  outward_normal_.resize(0);
  face_of_element_.clear();

  const Mesh& cpyMesh = *backMesh ;
  copy_levelset_sign( ls, ls_sign);

  std::vector<RemumberVertexPairT> zero_vertex_uses; // used to renumber the zero_vertexes
  const Uint nb_vertex_K = Element::nv;
  double loc_ls[nb_vertex_K];
  byte   loc_ls_sign[nb_vertex_K];

   for (int k=0; k<backMesh->nbElmts(); k++) {                      // loop over elements

    const typename Mesh::Element & K(cpyMesh[k]);

    for (Uint i= 0; i < K.nv; ++i) {
      loc_ls_sign[i] = ls_sign[cpyMesh(K[i])];
      loc_ls     [i] = ls     [cpyMesh(K[i])];
    }

    const RefPatch& cut =  RefPatch::instance( loc_ls_sign);

    if (cut.empty()) continue;

    for (typename RefPatch::const_face_iterator it= cut.face_begin(), end= cut.face_end();
    	 it != end; ++it) {
         face_of_element_[k] = element_of_face_.size();
         faces_.push_back( make_face(*it, K, loc_ls, zero_vertex_uses));
         element_of_face_.push_back(k);
         outward_normal_.push_back(make_normal(K, loc_ls));
    }


   }

}


template<typename M>
typename GenericInterface<M>::Rd
GenericInterface<M>::make_normal (const typename Mesh::Element& K, const double lset[Element::nv]) {

  Rd grad[Element::nv];
  K.Gradlambda(grad);
  Rd normal_ls;
  for(int i=0; i<Element::nv;++i) {

    normal_ls += grad[i]*lset[i];
  }
  normal_ls /= normal_ls.norm();
  return normal_ls;
}

template<typename M>
const typename GenericInterface<M>::FaceIdx
GenericInterface<M>::make_face (const typename RefPatch::FaceIdx& ref_tri,
				  const typename Mesh::Element& K,
				  const double lset[Element::nv],
				  std::vector<RemumberVertexPairT>& zero_vertex_uses)
{

  Uint loc_vert_num;
  Uint triIdx[nve];

  for (Uint j= 0; j < nve; ++j) {
    loc_vert_num= ref_tri[j];
    if (loc_vert_num < K.nv) {                            // zero vertex

      const Uint idx = (*backMesh)(K[loc_vert_num]);
      Rd Q = (*backMesh)(idx);
      vertices_.push_back(Q);
      triIdx[j] = vertices_.size() - 1;

    }
    else { // genuine edge vertex

      const Ubyte i0 = Mesh::Element::nvedge[loc_vert_num - K.nv][0],
      	i1 = Mesh::Element::nvedge[loc_vert_num - K.nv][1];
      const double t = lset[i0]/(lset[i0] - lset[i1]);
      Rd Q = (1.0 - t) * ((Rd) K[i0]) + t * ((Rd) K[i1]); // linear interpolation
      vertices_.push_back(Q);
      triIdx[j] = vertices_.size() - 1;

    }
  }
  return FaceIdx(triIdx);
}




// // To build interface using imported mesh/data from matlab
// template<typename M>
// void GenericInterface<M>::make_patch(const KN<double>& lls_sign,
// 				     const KN<int>& cut_element,
// 				     const KN<R>& nodex,
// 				     const KN<R>& nodey) {
//
//   assert(backMesh);
//   faces_.resize( 0);                          // reinitialize arrays
//   is_boundary_face_.resize( 0);
//   vertices_.resize(0);
//   element_of_face_.resize(0);
//   outward_normal_.resize(0);
//   face_of_element_.clear();
//
//
//
//   const Mesh& cpyMesh = *backMesh ;
//   copy_levelset_sign( lls_sign, ls_sign);
//
//
//   std::vector<RemumberVertexPairT> zero_vertex_uses; // used to renumber the zero_vertexes
//   const Uint nb_vertex_K = Element::nv;
//   double loc_ls_sign[nb_vertex_K];
//
//
//   for (int k=0; k<cut_element.size(); k++) {
//
//     int idxk = cut_element(k);
//     const typename Mesh::Element & K(cpyMesh[idxk]);
//
//     for (Uint i= 0; i < K.nv; ++i) {
//       loc_ls_sign[i] = lls_sign[cpyMesh(K[i])];
//     }
//
//     face_of_element_[idxk] = element_of_face_.size();
//     faces_.push_back( make_face(k, nodex, nodey));
//     element_of_face_.push_back(idxk);
//     outward_normal_.push_back(make_normal_noGrad(K, loc_ls_sign));
//   }
//
// }

//
// template<typename M>
// const typename GenericInterface<M>::FaceIdx
// GenericInterface<M>::make_face (int k,
// 				const KN<R>& nodex,
// 				const KN<R>& nodey)
// {
//   Uint triIdx[nve];
//    for (Uint j= 0; j < nve; ++j) {
//      Rd Q(nodex(2*k+j), nodey(2*k+j));
//      vertices_.push_back(Q);
//      triIdx[j] = vertices_.size() - 1;
//    }
//   return FaceIdx(triIdx);
//  }





template<typename M>
double GenericInterface<M>::distance(Rd P) const {
	double l = 1e300;
	for(auto it = vertex_begin(); it!=vertex_end();++it) {
		Rd D(P,*it);
		l = min(l, D.norme2());
	}
	return sqrt(l);
};




#endif
