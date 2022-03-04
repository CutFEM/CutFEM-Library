#ifndef CUT_METHOD_HPP_
#define CUT_METHOD_HPP_

#include "../util/util.hpp"
#include "HashTable.hpp"
#include "dataStruct2D.hpp"
#include "dataStruct3D.hpp"
#include "geometry.hpp"
enum ElementSignEnum { AllElement, NegElement, PosElement, NoElement};



template<typename T>
class SignElement{

  typedef T Element;
  static const int nve = Element::nv;

  int sum_;
public:
  SignElement() : sum_(0) {}
  SignElement(const double ls[nve]) : sum_(0) {
    for (Ubyte i= 0; i < nve; ++i)
      sum_ += util::fsign( ls[i]);
  }
  SignElement(const byte ls[nve]) : sum_(0) {
    for (Ubyte i= 0; i < nve; ++i)
      sum_ += ls[i];
  }


  bool cut() const  {return abs(sum_) != nve;}
  byte sign() const { return (sum_ == nve) ? 1 : ((sum_ == -nve)? -1 : 0);}
  int sum() const {return sum_;}
  int nb_node_positif() const {
    assert(0);
  }
};


// This class that save information obout the vertix and edge cut
// from a pattern
template<typename T>
class SignPattern
{
  typedef T Element;
  static const int nve = Element::nv;
  static const int dim = Element::Rd::d;

protected:
public:
  Ubyte num_root_vert_;  ///< number of vertices, where the level set function is zero.
  Ubyte num_root_;       ///< number of roots of the level set function; num_root_vert <= num_root
  byte sign_[Element::nv];  //Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$
  Ubyte cut_simplex_[Element::nvc]; ///< local number with respect to the reference tetra of the object on the cut: [0,num_root_vert): vertex numbers in (0..2); [num_root_vert, num_root): edge numbers in (0..3). Both parts are sorted in increasing order.
  Ubyte cut_simplex_rep_[Element::nvc]; ///< local number of the object on the cut: (0..5)

  void compute_cuts ();
  void get_ordered_list_node(std::vector<int>& list_cut) const ;
  int find_other_cutEdge_on_face(const int e0, int& f0) const ;
  void ordered_list_node();

public:
  SignElement<T> sign_element;

  SignPattern () : num_root_vert_( 3), num_root_( 2) {} ///< Uninitialized default state
  SignPattern (const byte   ls[nve]) :sign_element(ls) { assign( ls); } ///< Assign the sign pattern on the vertices; throws if ls is identically 0.
  SignPattern (const double ls[nve]) :sign_element(ls) { assign( ls); } ///< Assign a sign pattern on the vertices; throws if ls is identically 0.
  void assign (const byte   ls[nve]); ///< Assign a sign pattern on the vertices; throws if ls is identically 0.
  void assign (const double ls[nve]); ///< Assign a sign pattern on the vertices; throws if ls is identically 0.

  byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

  bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.

  bool no_zero_vertex () const { return num_root_vert_ == 0; } ///< True, iff there is no vertex, in which ls vanishes.

  Ubyte num_cut_simplexes () const { return num_root_; } ///< Number of edges and vertices with a root of ls.
  Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the tetra that are roots of ls.

  bool is_nd () const { return num_root_vert_ == dim+1; }
  bool is_hyperPlan () const { return num_root_ > dim-1; }

  /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..5, and () returns an extended vertex number in 4..9.
  ///@{
  Ubyte operator[] (int i) const { return cut_simplex_[i]; }
  Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
  ///@}


};

template<typename T>
void SignPattern<T>::assign (const byte ls[Element::nv]) {
  num_root_vert_= num_root_= 0;

  byte sum= 0;
  for (Ubyte i= 0; i < Element::nv; ++i) sum += (sign_[i] = ls[i]);
  if (sum == Element::nv || sum == -Element::nv) return;
  compute_cuts ();
}

template<typename T>
void SignPattern<T>::assign (const double ls[Element::nv]) {
  num_root_vert_= num_root_= 0;

  byte sum= 0;
  for (Ubyte i= 0; i < Element::nv; ++i){
    byte ss = util::sign(ls[i]); //(ls[i] >0)? 1:-1;//
    sign_[i] = ss;
    sum += ss;//(sign_[i] = sign(ls[i]));
  }

  if (sum == Element::nv || sum == -Element::nv) // optimize the case of uncut tetras
    return;
  compute_cuts();
}

template<typename T>
void SignPattern<T>::compute_cuts () {
  for (Ubyte i= 0; i < Element::nv; ++i){
      if (sign( i) == 0) {
      cut_simplex_[num_root_vert_++]= i;
    }
    num_root_= num_root_vert_;
  }
  for (Ubyte i= 0; i < Element::ne; ++i){
    if (sign( Element::nvedge[i][0])*sign( Element::nvedge[i][1]) == -1) {
      cut_simplex_[num_root_++]= i;
    }
  }

  // std::cout << " num root\t" << static_cast<int>(num_root_) << std::endl;
  // for(int i=0;i<num_root_;++i) {
  //   std::cout << static_cast<int>(cut_simplex_[i]) << "\t" ;
  // }
  // std::cout << std::endl;

  if(Element::nv == 8) ordered_list_node(); // for hexa
  std::memcpy( cut_simplex_rep_, cut_simplex_, Element::nvc*sizeof(byte));
  for (int i= num_root_vert_; i < num_root_; ++i) cut_simplex_rep_[i]+= Element::nv;


  // std::cout << " num root\t" << static_cast<int>(num_root_) << std::endl;
  // for(int i=0;i<num_root_;++i) {
  //   std::cout << static_cast<int>(cut_simplex_[i]) << "\t" ;
  // }
  // std::cout << std::endl;

}

template<typename T>
void SignPattern<T>::ordered_list_node(){

  Ubyte list_cut[Element::nvc];
  // INITIAL EDGE && FACE
  const int e0 = this->cut_simplex_[0];
  const int f0 = T::faceOfEdge[e0][0];

  // FIND THE OTHER CUT ON THIS FACE
  int e_previous = e0;
  int e_next = -1;
  int f_next = f0;
  int i=0;
  while(e_next != e0) {
    list_cut[i] = static_cast<Ubyte>(e_previous);
    e_next = find_other_cutEdge_on_face(e_previous, f_next);
    e_previous = e_next;
    ++i;
  }
  assert(i == num_root_);
  for(int j=0;j<num_root_;++j) cut_simplex_[j] = list_cut[j];

}

template<typename T>
void SignPattern<T>::get_ordered_list_node(std::vector<int>& list_cut) const {

  // INITIAL EDGE && FACE
  const int e0 = this->cut_simplex_[0];
  const int f0 = T::faceOfEdge[e0][0];
  list_cut.resize(0);

  // FIND THE OTHER CUT ON THIS FACE
  int e_previous = e0;
  int e_next = -1;
  int f_next = f0;
  while(e_next != e0) {
    list_cut.push_back(e_previous);
    e_next = find_other_cutEdge_on_face(e_previous, f_next);
    e_previous = e_next;
  }
}

template<typename T>
int SignPattern<T>::find_other_cutEdge_on_face(const int e0, int& f0) const {
  int e_next = -1;
  int f_next = -1;

  for(int i=0;i<num_root_;++i) {
    // assert(this->cut_simplex_rep_[i] >= T::nv);

    int e = this->cut_simplex_[i];
    if(e == e0) continue;
    // check if edge e is on the face f0
    if(T::faceOfEdge[e][0] == f0 ){
      e_next = this->cut_simplex_[i];
      f_next = T::faceOfEdge[e][1];
      break;
    }
    else if(T::faceOfEdge[e][1] == f0) {
      e_next = this->cut_simplex_[i];
      f_next = T::faceOfEdge[e][0];
      break;
    }
  }

  assert(e_next != -1 && f_next != -1);
  f0 = f_next;
  return e_next;
}



template<typename E>
int instance_idx(const byte ls[E::nv]) {
  int s =0;
  for(int i=0,j=E::nv-1;i<E::nv;++i, --j){
    s += pow(3,i)*ls[j];
  }
  return  s;
}



// Class that save all the possible pattern of sign on an elements
// and create the cut surface
template<typename Element>
class RefPatch {
public:
  static const int nv = Element::nv;
  static const int dim = Element::Rd::d;

  typedef SortArray<Ubyte, dim> FaceIdx; ///< the vertices of a triangle of the cut: the tetra's vertices are denoted by 0..3, the edge-cuts by edge-num + 4, which is in 4..9.
  typedef const FaceIdx* const_face_iterator;
  typedef       FaceIdx*       face_iterator;

  /// \brief Initializes RefTetraPatchCL::instance_array_ below in the constructor by calling RefTetraPatchCL::instance() for all legal sign patterns.
  class InitializerCL
  {
  private:
    static int init_count_;

  public:
    InitializerCL ();
  };
  public:
    RefPatch () : size_( static_cast<Ubyte>( -1)), is_boundary_face_( 0) {
      InitializerCL ();
    } ///< Uninitialized default state
    RefPatch (const SignPattern<Element>& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignPattern<Element>& cut); ///< Assign a sign pattern on the vertices; returns the value of empty()


private:
  // Save all the possible pattern.
  // Each node has 3 possible signs (-1, 0, 1)
  static RefPatch<Element> instance_array_[Element::nb_sign_pattern];

  FaceIdx face_[4];    ///< at most two triangles
  Ubyte size_;             ///< number of triangles
  Ubyte is_boundary_face_; ///< true if the triangle is one of the tetra's faces.

  Ubyte num_elements (const SignPattern<Element>& cut) const {
    return cut.is_hyperPlan() ? cut.num_cut_simplexes() - (dim -1) : 0; }



public:
  bool  is_initialized () const { return size_ == 1; } ///< True after assign(...)
//
  ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
  static inline const RefPatch<Element>& instance (const byte   ls[Element::nv]) {
    int m0 = (Element::nb_sign_pattern -1)/2;
    RefPatch<Element>& instance = instance_array_[instance_idx<Element> ( ls) + m0];
    if ( !instance.is_initialized()) {
      instance.assign( SignPattern<Element>( ls));
    }
    return instance;
  }
  static inline const RefPatch<Element>& instance (const double ls[3])
  {
    byte ls_byte[nv];
    std::transform( ls + 0, ls + nv, ls_byte + 0, util::sign);
    return instance( ls_byte);
  }
  ///@}

  // bool is_boundary_triangle () const { return is_boundary_triangle_ == 1; } ///< true, iff the triangle is one of the tetra's faces.

  bool  empty () const { return size_ == 0; } ///< true, iff the area of the intersection is 0.
  size_t size () const { return size_; }      ///< Number of triangles, 0, 1, or 2

  ///@{ Random-access to the triangles
  const_face_iterator face_begin () const { return face_; }
  const_face_iterator face_end   () const { return face_ + size_; }
  ///@}
};

template<typename E>
RefPatch<E> RefPatch<E>::instance_array_[E::nb_sign_pattern];

template<typename E> int RefPatch<E>::InitializerCL::init_count_= 0;



// Class that take a sign pattern and build and array with the elements
// in the cut.
// The index given give vertices or edges
template<typename E>
class RefPartition {
public:
  static const int nv = E::nv;                         // cut creating triangles
  static const int dim = E::d;
  static const int nt_in_notcut_K = E::nb_ntcut;
  static const int max_nb_element = 20;
  static const int start_array = 10;
  typedef SortArray<Ubyte, dim+1> ElementIdx;           // the vertices of a triangle of the cut:

  typedef const ElementIdx* const_element_iterator;
  typedef       ElementIdx*       element_iterator;

  class InitializerCL {
  private:
    static int init_count_;

  public:
    InitializerCL ();
  };

  ElementIdx elements_[max_nb_element];                          // upper bound => Hexa cut made of 10 tetra
  element_iterator begin_;
  element_iterator end_;
  bool is_cut_ = false;

  //  all possible sign-patterns on the vertices
  // each node has 3 possible states (-1, 0, 1)
  static RefPartition<E> instance_array_[E::nb_sign_pattern];

  void AddElement (Ubyte v[E::d+1], int sign) {
    (sign == -1 ? *--begin_ : *end_++) = ElementIdx( v );
  }
  // e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
  void AddQuadrilateral (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, int sign) {
    Ubyte list_v1[] = {e0, e1, f1};
    AddElement( list_v1, sign);
    Ubyte list_v2[] = {e0, f1, f0};
    AddElement( list_v2, sign);
  }
  // e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
  void AddPentagone (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g1, int sign) {
    Ubyte list_v0[] = {e0, e1, g1};
    AddElement( list_v0, sign);

    Ubyte list_v1[] = {f0, e0, g1};
    AddElement(list_v1, sign);
    Ubyte list_v2[] = {f0, f1, g1};
    AddElement( list_v2, sign);
  }
  // e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
  void AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g0, Ubyte g1, int sign) {
    Ubyte list_v0[] = {e0, e1, f1, g1};
    AddElement( list_v0, sign);
    Ubyte list_v1[] = { e0, f0, f1, g1};
    AddElement( list_v1, sign);
    Ubyte list_v2[] = {e0, f0, g0, g1};
    AddElement( list_v2, sign);
  }

  void AddHexa (Ubyte e0, Ubyte e1, Ubyte e2, Ubyte e3, Ubyte e4, Ubyte e5, Ubyte e6, Ubyte e7,int sign) {
    Ubyte list_v0[] = {e0, e5, e6, e4};
    AddElement( list_v0, sign);
    Ubyte list_v1[] = { e0, e3, e5, e1};
    AddElement( list_v1, sign);
    Ubyte list_v2[] = {e3, e5, e6, e7};
    AddElement( list_v2, sign);
    Ubyte list_v3[] = {e0, e3, e6, e2};
    AddElement( list_v3, sign);
    Ubyte list_v4[] = {e0, e3, e5, e6};
    AddElement( list_v4, sign);
  }



  // If the intersection is quadrilateral, this returns the first uncut edge.
  // It is always one of the three first edges
  Ubyte first_uncut_edge (const SignPattern<E>& cut) const {
    return cut[0] == 1 ? 0 : (cut[1] == 2 ? 1 : 2); }

//
public :
  RefPartition () : begin_( elements_ + max_nb_element), end_( elements_) {

    InitializerCL ();
  }

  // Initialize with sign pattern on the vertices
  RefPartition (const SignPattern<E>& cut)  {assign( cut); }

  // Assign a sign pattern on the vertices; returns the value of is_uncut()
  bool assign (const SignPattern<E>& cut);

  // True after assign(...)
  bool is_initialized () const {return begin_ < end_;}

  ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
  static inline const RefPartition& instance (const byte   ls[E::nv]) {
    int m0 = (E::nb_sign_pattern -1)/2;
    RefPartition& instance = instance_array_[instance_idx<E> (ls) + m0];
    if ( !instance.is_initialized()) {
      instance.assign( SignPattern<E>( ls));
    }
    return instance;
  }
  static inline const RefPartition& instance (const double ls[E::nv])
  {
    byte ls_byte[E::nv];
    std::transform( ls + 0, ls + E::nv, ls_byte + 0, util::sign);
    return instance( ls_byte);
  }
  ///@}

  // True, iff the partition has exactly one tetra
  // bool is_uncut () const { return end_ == begin_ + nt_in_notcut_K; }
  bool is_uncut () const { return !is_cut_; }

  bool is_cut () const {
    return !is_uncut(); }
  // Sign of the tetra, to which t points
  int whatSign (const_element_iterator t) const { return t < elements_ + start_array ? -1 : 1; }
  int whatSign () const { return begin_ < elements_ + start_array ? -1 : 1; }


    ///@{ Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum
  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  { if(s == NoElement) return end_;
    else return s == PosElement ? elements_ + start_array : begin_; }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return s == NegElement ? elements_ + start_array : end_; }
    ///@}


};

template<typename E> RefPartition<E> RefPartition<E>::instance_array_[E::nb_sign_pattern];

template<typename E> int RefPartition<E>::InitializerCL::init_count_= 0;




class CutElement2 {
  typedef typename Mesh2::Rd Rd;
  typedef SortArray<int, 3> TriaIdx;

  public :
  static const int nv = 3;
  Rd  vertices[nv];
  int lab;

  CutElement2() {}
  CutElement2(const std::vector<Rd>& v0, const TriaIdx& iv,int r=0) {
    for(int i=0;i<nv;++i) vertices[i]=v0[iv[i]];
    lab=r;
  }
  const R2& operator[](int i) const {return vertices[i];}
  // R2& operator[](int i) {return vertices[i];}

  void print() const  {
    for(int i=0;i<nv;++i) {
      std::cout << vertices[i] << std::endl;
    }
  }
};

template<typename E>
class Virtual_Partition {
public:

  typedef typename E::Rd Rd;
  virtual void get_list_node (vector<Rd>& node, int s) const = 0;
  virtual CutElement2 get_element(int k) const = 0;
  virtual int nb_element() const = 0;
};

// Class that does computation on cut elements from refPartition
template<typename E>
class Partition : public Virtual_Partition<E>{

public :
  typedef SortArray<Ubyte, E::Rd::d+1> ElementIdx;           // the vertices of a triangle of the cut:
  typedef const ElementIdx* const_element_iterator;
  typedef       ElementIdx*       element_iterator;
  typedef typename E::Rd Rd;
  // typedef typename TypeCutData<Rd::d>::CutData CutData;
  typedef E Element;

  const RefPartition<E>& patch;
  const Element& T;
  double ls[E::nv];
  // const CutData* cutData = nullptr;

  Partition(const Element& t) : patch(RefPartition<E>::instance_array_[0]), T(t) {}

  Partition(const Element& t, const double lls[E::nv])
    : patch(RefPartition<E>::instance(lls)), T(t) {
      for(int i=0;i<E::nv;++i) ls[i] = lls[i];
  }
  // Partition(const Element& t, const CutData& data)
  //   : patch(RefPartition<E>::instance(data.sign)), T(t), cutData(&data) {
  //     for(int i=0;i<E::nv;++i) ls[i] = static_cast<double>(data.sign[i]);
  // }


  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  { return patch.element_begin(s); }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return patch.element_end(s); }
//
  int whatSign (const_element_iterator t) const { return patch.whatSign(t); }
  bool is_cut () const { return patch.is_cut(); }
  bool isnot_cut () const { return patch.is_uncut(); }



  ElementSignEnum what_part(const int the_domain) const {
    if(the_domain == -1) return AllElement;
    return (is_cut()) ? ( (the_domain == 0)? PosElement : NegElement  )
      : (((the_domain==0 && ls[0]>0) || (the_domain==1 && ls[0]<0))? AllElement : NoElement);
  }

//   double getEdgeLength() const {
//     return 1./3*(T.lenEdge(0)+ T.lenEdge(1)+T.lenEdge(2));
//   }


  R mesure(const_element_iterator it) const {
    if( patch.is_uncut() && E::nb_ntcut==1) return T.mesure();
    Rd N[Rd::d+1];

    for(int i=0; i<Rd::d+1;++i) {  // triangle or tet
      Uint idx = (*it)[i];
      if(idx < E::nv) N[i] = T[idx];
      else{
        const Ubyte v0= Element::nvedge[idx - E::nv][0], v1= Element::nvedge[idx - E::nv][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
      }
    }
    return geometry::mesure_simplex(N);
  }

//
//     return std::fabs(det(N[0],N[1],N[2]))*0.5;
//   }
//
//   R mesure(const int domain) const {
//
//     if( patch.is_uncut()) return T.mesure();
//
//     R mes = 0;
//     ElementSignEnum the_part = what_part(domain);
//
//     for(const_element_iterator it = element_begin(the_part);
//     it != element_end(the_part); ++it) {
//
//       Rd N[3];
//       for(int i=0; i<3;++i) {
//         Uint idx = (*it)[i];
//         if(idx < 3) {N[i] = T[idx];}
//         else {
//           if(cutData) {
//             N[i] = cutData->pointFromEdge(idx-3);
//           }
//           else {
//             const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
//             const R t = -ls[v0]/(ls[v1]-ls[v0]);
//             N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
//           }
//         }
//       }
//     mes += std::fabs(det(N[0],N[1],N[2]))*0.5;
//   }
//   return mes;
// }
//
//   R mesureEdge(int e, int dom) const {
//     assert(cutData);
//     if( patch.is_uncut()) return T.lenEdge(e);
//     int s = (dom == 0)? 1 : -1;
//
//     const Ubyte v0= Element::nvedge[e][0], v1= Element::nvedge[e][1];
//     // const R t = -ls[v0]/(ls[v1]-ls[v0]);
//     // R2 A = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
//     R2 A = cutData->pointFromEdge(e);
//     R2 B = (s * ls[v0] > 0) ? T[v0] : T[v1];
//     R2 E(A,B);
//     return E.norm();
//   }
//
//   Rd toEdge(int e, R1 ip, int dom) const {
//     assert(cutData);
//     if( patch.is_uncut()) assert(0);
//     int s = (dom == 0)? 1 : -1;
//
//     const Ubyte v0= Element::nvedge[e][0], v1= Element::nvedge[e][1];
//     R2 A = cutData->pointFromEdge(e);
//     R2 B = (s * ls[v0] > 0) ? T[v0] : T[v1];
//     R2 E = (1-ip.x) * A + ip.x * B ;
//
//     return E;
//   }
//
//   Rd toK(const_element_iterator it, const Rd Phat) const {
//
//     if( patch.is_uncut()) return T(Phat);
//     Rd N[3];
//     for(int i=0; i<3;++i) {
//       Uint idx = (*it)[i];
//       if(idx < 3) N[i] = T[idx];
//       else {
//         if(cutData) {
//           N[i] = cutData->pointFromEdge(idx-3);
//         }
//         else {
//           const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
//           const R t = -ls[v0]/(ls[v1]-ls[v0]);
//           N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
//         }
//       }
//     }
//     return  (1-Phat.x - Phat.y) * N[0] + Phat.x * N[1]
//     + Phat.y * N[2];
//
//   }

  Rd get_Vertex(const_element_iterator it, const int i) const {

    Uint idx = (*it)[i];
    if(idx < E::nv) return (Rd) T[idx];
    else {
      // if(cutData) {
      //   // return cutData->pointFromEdge(idx-E::nv);
      // }
      // else {
        const Ubyte v0= Element::nvedge[idx - E::nv][0], v1= Element::nvedge[idx - E::nv][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        return (1-t) * ((Rd)T[v0]) + t * ((Rd) T[v1]) ;
      // }
    }
  }

private:
  Rd get_vertex(const int idx) const {

    if(idx < E::nv) return (Rd) T[idx];
    else {
      // if(cutData) {
      //   // return cutData->pointFromEdge(idx-E::nv);
      // }
      // else {
        const Ubyte v0= Element::nvedge[idx - E::nv][0], v1= Element::nvedge[idx - E::nv][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        return (1-t) * ((Rd)T[v0]) + t * ((Rd) T[v1]) ;
      // }
    }
  }
  //   R get_LSvalue(const_element_iterator it, const int i) const {
//     Uint idx = (*it)[i];
//     if(idx < 3) return ls[idx];
//     else {
//       const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
//       const R t = -ls[v0]/(ls[v1]-ls[v0]);
//       return (1-t) * ls[v0] + t * ls[v1] ;
//     }
//   }
public:
void get_list_node (vector<Rd>& node, int s) const ;
CutElement2 get_element(int k) const { assert(0); return CutElement2();};
int nb_element() const {return patch.end_-patch.begin_;}

};


template<typename E>
class Physical_Partition : public Virtual_Partition<E>{

  typedef SortArray<int, E::Rd::d+1> ElementIdx;           // the vertices of a triangle of the cut:
  typedef const ElementIdx* const_element_iterator;
  typedef       ElementIdx*       element_iterator;
  typedef typename E::Rd Rd;
  typedef E Element;

public:
  const Element& T;



private:
  std::vector<Rd>       vertices_;
  std::vector<ElementIdx>  elements_idx;

public:


  Physical_Partition(const Element& TT) :T(TT) {}
  // void reset() {
  //   vertices_.resize(0);
  //   elements_idx.resize(0);
  // }
  void add_node(const Rd P) { vertices_.push_back(P);}
  void add_element(const ElementIdx& iv) {
    elements_idx.push_back(iv);
  }
  void erase_element(int k) {
    elements_idx.erase(elements_idx.begin()+k);
  }
  int nb_element() const {return elements_idx.size();}
  int nb_node() const {return vertices_.size();}

  CutElement2 get_element(int k) const { return CutElement2(vertices_,elements_idx[k]);}


  void get_list_node (vector<Rd>& node, int s) const {
    assert(0);
    node = vertices_;
  }
};








#endif
