#include "GenericInterface.hpp"


template<typename M>
class Interface  {

public :
  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;

  static const int nve = Rd::d;
  typedef FaceInterface<nve> Face;

  typedef const Face* const_face_iterator;
  typedef       Face*       face_iterator;
  typedef const Rd* const_vertex_iterator;
  typedef       Rd*       vertex_iterator;


  const Mesh* backMesh;

  std::vector<Face>  faces_;            // All triangles of the interface.
  std::vector<Rd> vertices_;
  std::vector<Rd> outward_normal_;
  std::vector<Uint> element_of_face_;
  std::map<int, int> face_of_element_;
  std::vector<Uint> edge_of_node_;

private :
  // const Face  make_face (const typename RefPatch::FaceIdx& ref_tri,
	// 		    const typename Mesh::Element& K,
	// 		    const double lset[Element::nv],
	// 		    std::vector<RemumberVertexPairT>& renumber_zero_verts, int label);
  // const FaceIdx  make_face (int k,
	// 		    const KN<R>& nodex,
	// 		    const KN<R>& nodey, int label);
  // const FaceIdx  make_face (int k,
  // 			    const Marker& marker);

  // Rd make_normal(const typename Mesh::Element& K, const double lset[Element::nv]);
  // Rd make_normal_noGrad(const typename Mesh::Element& K, const double lset[Element::nv]);

public :
  // GenericInterface() : backMesh(nullptr) { }
  Interface(const Mesh & MM) : backMesh(&MM) { }
  // GenericInterface(const Mesh & MM, const KN<double>& ls, int label) : backMesh(&MM){
	// 	make_patch(ls, label);
	// }

  // void make_patch (const Mesh & MM, const KN<double>& ls, int label);
  // void make_patch (const KN<double>& ls, int label);
  // void add (const KN<double>& ls, int label);


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


  // virtual Rd mapToFace(const Face& f, const typename Element::RdHatBord x ) const = 0;
  // virtual Rd computeDx(const Face& f) const = 0;
  // virtual CutData getCutData(const int k) const = 0;

	~Interface() {}
private:
  Interface(const Interface &); // pas de construction par copie
  void operator=(const Interface &);// pas affectation par copy

};
