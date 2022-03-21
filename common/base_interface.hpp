#ifndef BASE_INTERFACE_HPP_
#define BASE_INTERFACE_HPP_

#include "GenericInterface.hpp"

class FunFEMVirtual {
public :

double * data = nullptr;
KN_<double> v;


FunFEMVirtual () :v(data, 0) {}
FunFEMVirtual (int df) : data(new double[df]), v(data,df) {v = 0.;}
FunFEMVirtual (KN_<double>&u) : v(u) {}

virtual double eval(const int k, const R* x, int cu=0, int op=0) const  = 0;
virtual double eval(const int k, const R* x, const R t, int cu, int op, int opt) const = 0;
virtual int idxElementFromBackMesh(int, int=0) const = 0;
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

  virtual void cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part) const = 0;
  virtual R measure(const Face& f) const = 0;
  virtual bool isCutFace(int k, int ifac) const = 0;

  // virtual Rd mapToFace(const Face& f, const typename Element::RdHatBord x ) const = 0;
  // virtual Rd computeDx(const Face& f) const = 0;
  // virtual CutData getCutData(const int k) const = 0;

	~Interface() {}
private:
  Interface(const Interface &); // pas de construction par copie
  void operator=(const Interface &);// pas affectation par copy

};

template<typename Mesh>
class Time_Interface {
public:
	// typedef FunFEM<Mesh> Fun_h;
private:
	KN<Interface<Mesh>*> interface;
	int n;

public:

  Time_Interface(int nt) : interface(nt), n(nt) {
    for(int i=0;i<n;++i){ interface[i] = nullptr;}
  }

  void init(int i, const Mesh & Th, const FunFEMVirtual& ls);
  void init(const Mesh & Th, const KN<FunFEMVirtual>& ls);

  Interface<Mesh>* operator[](int i) const{
    assert(0 <= i && i < n);
    return interface[i];
  }
  Interface<Mesh>* operator()(int i) const{
    assert(0 <= i && i < n);
    return interface[i];
  }

  int size() const { return n;}

  ~Time_Interface(){
    for(int i=0;i<n;++i){
      if(interface[i]) delete interface[i];
    }
  }

private:
  Time_Interface(const Time_Interface&);
  void operator=(const Time_Interface &);
};


#endif
