#ifndef BASE_INTERFACE_HPP_
#define BASE_INTERFACE_HPP_


#include <iostream>
#include <cassert>
#include <bitset>
#include "HashTable.hpp"
#include "RNM.hpp"
#include "Label.hpp"
#include "../util/util.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "cut_method.hpp"




#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif
// #include "GenericInterface.hpp"
#include "../FESpace/QuadratureFormular.hpp"

class FunFEMVirtual {
public :

double * data = nullptr;
KN_<double> v;


FunFEMVirtual () :v(data, 0) {}
FunFEMVirtual (int df) : data(new double[df]), v(data,df) {v = 0.;}
FunFEMVirtual (KN_<double>&u) : v(u) {}

virtual double eval(const int k, const R* x, int cu=0, int op=0) const  = 0;
virtual double eval(const int k, const R* x, const R t, int cu, int op, int opt) const = 0;
virtual double evalOnBackMesh(const int k, int dom, const R* x, int cu=0, int op=0) const  = 0;
virtual double evalOnBackMesh(const int k, int dom, const R* x, const R t, int cu, int op, int opt) const = 0;
virtual int idxElementFromBackMesh(int, int=0) const = 0;
};



template<int N>
struct FaceInterface {

};
template<>
struct FaceInterface<2> : public  SortArray<Uint, 2>, public Label {
  typedef SortArray<Uint, 2> FaceIdx;

  FaceInterface(const Uint& a0,const Uint &a1, int l=0)
  : FaceIdx(a0,a1), Label(l) {}
  FaceInterface(Uint *a, int l =0)   : FaceIdx(a), Label(l) {}
  FaceInterface() : FaceIdx(), Label(0) {}


};
template<>
struct FaceInterface<3> : public  SortArray<Uint, 3>, public Label {
  typedef SortArray<Uint, 3> FaceIdx;

  FaceInterface(const Uint& a0,const Uint &a1,const Uint &a2, int l=0)
  : FaceIdx(a0,a1,a2), Label(l) {}
  FaceInterface(Uint *a, int l =0)   : FaceIdx(a), Label(l) {}
  FaceInterface() : FaceIdx(), Label(0) {}

};


template<typename M>
class Interface  {

public :
  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;

  static const int nve = Rd::d;
  typedef FaceInterface<nve> Face;
  typedef SortArray<Ubyte, Element::Rd::d+1> ElementIdx;


  typedef const Face* const_face_iterator;
  typedef       Face*       face_iterator;
  typedef const Rd* const_vertex_iterator;
  typedef       Rd*       vertex_iterator;


  const Mesh* backMesh;

  std::vector<Face>  faces_;
  std::vector<Rd> vertices_;
  std::vector<Rd> outward_normal_;
  std::vector<Uint> element_of_face_;
  std::map<int, int> face_of_element_;
  std::vector<Uint> edge_of_node_;

public :
  Interface(const Mesh & MM) : backMesh(&MM) { }


  Rd operator()(const int k, const int i) const {return vertices_[faces_[k][i]];}
  const Rd& operator()(const int i) const {return vertices_[CheckV(i)];}
  const Face& operator[](const int k) const {return faces_[CheckT(k)];}
  const Face& getFace(const int k) const {return faces_[CheckT(k)];}

  int CheckV(int i) const { ASSERTION(i>=0 && i < vertices_.size()); return i;}
  int CheckT(int i) const { ASSERTION(i>=0 && i < faces_.size()); return i;}


  Uint idxElementOfFace(const int k) const { return element_of_face_[k];}
  Uint idxFaceOfElement(const int k) const {
    const auto it = face_of_element_.find(k);
    assert(it != face_of_element_.end());
    return it->second;}

  Uint nbElement  () const { return faces_.size(); }
  Rd normal(const int k) const { return outward_normal_[k];}
  bool isCut(const int k) const {
    return (face_of_element_.find(k) != face_of_element_.end());
  }
  const Element& get_element(int k) const{return (*backMesh)[k];}
  const Mesh& get_mesh() const {assert(backMesh);return *backMesh;}
  const_face_iterator face_begin () const { return (faces_.begin()).base(); }
  const_face_iterator face_end   () const { return (faces_.end()).base(); }

  const_vertex_iterator vertex_begin () const { return (vertices_.begin()).base(); }
  const_vertex_iterator vertex_end   () const { return (vertices_.end()).base(); }

  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(faces_.size());}
  virtual int next_element() const {  return MPIcf::next_element(faces_.size());}
  virtual int last_element() const {  return MPIcf::last_element(faces_.size());}

  #else
  int first_element() const { return 0;}
  int next_element() const {return 1;}
  int last_element() const { return faces_.size();}
  #endif

  virtual SignElement<Element> get_SignElement(int k) const =0;
  virtual Partition<Element> get_partition(int k) const = 0;
  virtual Partition<typename Element::Face> get_partition_face(const typename Element::Face& face, int k, int ifac) const = 0;

  virtual void cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part) const = 0;
  virtual R measure(const Face& f) const = 0;
  virtual R measure(int i) const {return measure(faces_[i]);};

  virtual bool isCutFace(int k, int ifac) const = 0;

  virtual Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x ) const = 0;
  // virtual Rd computeDx(const Face& f) const = 0;
  // virtual CutData getCutData(const int k) const = 0;

	~Interface() {}
private:
  Interface(const Interface &); // pas de construction par copie
  void operator=(const Interface &);// pas affectation par copy

};

template<typename Mesh>
class TimeInterface {
public:
	// typedef FunFEM<Mesh> Fun_h;
private:
	KN<Interface<Mesh>*> interface_;
	int n_;
  const QuadratureFormular1d& time_quadrature_;

public:

  // TimeInterface(int nt) : interface(nt), n(nt) {
  //   for(int i=0;i<n;++i){ interface[i] = nullptr;}
  // }

  TimeInterface(const QuadratureFormular1d& qTime) : interface_(qTime.n), n_(qTime.n), time_quadrature_(qTime) {
    for(int i=0;i<n_;++i){ interface_[i] = nullptr;}
  }
  TimeInterface(int nt) :  interface_(nt), n_(nt ), time_quadrature_(*Lobatto(exactLobatto_nPt(nt))) {
    for(int i=0;i<n_;++i){ interface_[i] = nullptr;}
  }


  void init(int i, const Mesh & Th, const FunFEMVirtual& ls);
  void init(const Mesh & Th, const KN<FunFEMVirtual>& ls);

  Interface<Mesh>* operator[](int i) const{
    assert(0 <= i && i < n_);
    return interface_[i];
  }
  Interface<Mesh>* operator()(int i) const{
    assert(0 <= i && i < n_);
    return interface_[i];
  }

  int size() const { return n_;}
  const QuadratureFormular1d& get_quadrature_time() const {return time_quadrature_;}

  ~TimeInterface(){
    for(int i=0;i<n_;++i){
      if(interface_[i]) delete interface_[i];
    }
  }

private:
  TimeInterface(const TimeInterface&);
  void operator=(const TimeInterface &);
};


#endif
