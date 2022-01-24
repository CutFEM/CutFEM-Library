#ifndef INTERFACE3DN_HPP_
#define INTERFACE3DN_HPP_

#include "GenericInterface.hpp"
#include "Mesh3dn.hpp"
#include "marker.hpp"
//
// /*
//  *   Represents the reference tetra, which is cut by a linear level set function ls.
//  *   The values of the latter are prescribed on the vertices.
//  */
//
// class SignPatternTrait3 : public GSignPatternTrait<Tet>
// {
//
// public :
//   bool is_3d () const { return num_root_vert_ == 4; }
//   bool is_2d () const { return num_root_ > 2; }  ///< True, iff the intersection has positive area.
//   friend std::ostream& operator<< (std::ostream&, const SignPatternTrait3&); ///< Debug-output to a stream (dumps all members)
//
//   SignPatternTrait3() : GSignPatternTrait<Tet>() {}
//   SignPatternTrait3(const byte    ls[4]) : GSignPatternTrait<Tet>(ls) {}
//   SignPatternTrait3(const double  ls[4]) : GSignPatternTrait<Tet>(ls) {}
//
// };
//
//
// /*
//  *  Return a signed array-index for the possible 3^4 sign-patterns on the vertices of a tetra.
//  *  The index ranges from [-40..40].
//  */
// inline static byte instance_idx3 (const byte ls[4])
// {
//   return  27*ls[0] + 9*ls[1] + 3*ls[2] + ls[3];
// }
//
// inline static Ubyte instance_idx3 (const double ls[4])
// {
//   return  27*sign( ls[0]) + 9*sign( ls[1]) + 3*sign( ls[2]) + sign(ls[3]);
// }
//
//
//
// /*
//  * The triangle of the intersection of the reference-tet with a linear levelset-function.
//  * The class memoizes used sign-patterns
//  */
// class RefPatch3
// {
// public:
//   static const int nv = 3;
//   typedef SortArray<Ubyte, nv> FaceIdx;           // the vertices of a triangle of the cut:
//                                                      //the tetra's vertices are denoted by 0..3,
//                                                      //the edge-cuts by edge-num + 4
//   typedef const FaceIdx* const_face_iterator;
//   typedef       FaceIdx*       face_iterator;
//
//   class InitializerCL
//   {
//   private:
//     static int init_count_;
//
//   public:
//     InitializerCL ();
//   };
//
// private:
//   FaceIdx face_[2];      ///< at most two triangles
//   Ubyte size_;                 ///< number of triangles
//   Ubyte is_boundary_face_; ///< true if the triangle is one of the tetra's faces.
//
//   Ubyte num_elements (const SignPatternTrait3& cut) const {
//     return cut.is_2d() ? cut.num_cut_simplexes() - 2 : 0; }
//   FaceIdx MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const {
//     return FaceIdx( v0, v1,v2); }
//
//   static RefPatch3 instance_array_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices
//
//
// public:
//   RefPatch3 () : size_( static_cast<Ubyte>( -1)), is_boundary_face_( 0) {}
//
//   ///< Initialize with sign pattern on the vertices
//   RefPatch3 (const SignPatternTrait3& cut) { assign( cut); }
//
//   ///< Assign a sign pattern on the vertices; returns the value of empty()
//   bool assign (const SignPatternTrait3& cut);
//
//   ///< True after assign(...)
//   bool  is_initialized () const { return size_ <=2; }
//
//
//   ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
//   static inline const RefPatch3& instance (const byte   ls[4]) {
//     RefPatch3& instance = instance_array_[instance_idx3 ( ls) + 40];
//
//     if ( !instance.is_initialized()) {
//       instance.assign( SignPatternTrait3( ls));
//     }
//     return instance;
//   }
//   static inline const RefPatch3& instance (const double ls[4])
//   {
//     byte ls_byte[4];
//     std::transform( ls + 0, ls + 4, ls_byte + 0, sign);
//     return instance( ls_byte);
//   }
//   ///@}
//
//   ///< true, iff the area of the intersection is 0.
//   bool  empty () const { return size_ == 0; }
//
//   ///< Number of triangles, 0, 1, or 2
//   size_t size () const { return size_; }
//
//   const_face_iterator face_begin () const { return face_; }
//   const_face_iterator face_end   () const { return face_ + size_; }
//
// };
//



class Interface3 : public GenericInterface<Mesh3>
{
typedef CutData3 CutData;
public:
  Interface3() : GenericInterface<Mesh3>() {}
  Interface3(const Mesh & MM) : GenericInterface<Mesh3>(MM) {}
  Interface3(const Mesh & MM, const KN<double>& ls, int label=0)
  : GenericInterface<Mesh3>(MM, ls, label) {}

  R3 mapToFace(const FaceIdx& f, const R2 PHat ) const {
    return (1-PHat.x- PHat.y)*vertices_[f[0]] + PHat.x *vertices_[f[1]]
      +PHat.y*vertices_[f[2]] ;
  }
  R3 computeDx(const FaceIdx& f) const {
    R3 u(vertices_[f[0]],vertices_[f[1]]);
    R3 v(vertices_[f[0]],vertices_[f[2]]);
    return 0.5 * (u ^ v);
  }


  virtual CutData getCutData(const int k) const {
    assert(backMesh);
    const Element& K((*backMesh)[k]);
    byte loc_sign[4];
    R3 point[4];
    int array_edge[6] = {-1,-1,-1,-1,-1,-1};
    for(int i=0;i<Element::nv;++i) {
      loc_sign[i] = ls_sign((*backMesh)(k, i));
    }

    const typename Mesh3::RefPatch& cut =  RefPatch::instance( loc_sign);

    if(isCut(k)){
      int size = cut.size();
       int i0 = (size == 2)? -1 : 0;
       int idPoint = 0;
       for (typename RefPatch::const_face_iterator it= cut.face_begin(), end= cut.face_end();
       it != end; ++it, ++i0) {

         int idxEdge[Mesh::nva];
         for (Uint j= 0; j < Mesh::nva; ++j) {
           idxEdge[j] = (*it)[j] ;
           if(idxEdge[j]>3)  idxEdge[j] -= 4;
         }
         int idxFace =  face_of_element_.at(k)+ i0 ;
         for (Uint j= 0; j < Mesh::nva; ++j) {
           if(array_edge[idxEdge[j]]  == -1) {
             array_edge[idxEdge[j]] = idPoint;
             point[idPoint] = (*this)(idxFace,j );
             idPoint += 1;
           }
         }
       }
       return CutData(loc_sign, point, array_edge);
     }
     else
     return CutData(loc_sign);
   }

  // R3 computeDx(const Element&K, const int i) const {
  //   R3 u(K[Element::nvface[i][0]], K[Element::nvface[i][1]]);
  //   R3 v(K[Element::nvface[i][0]], K[Element::nvface[i][2]]);
  //   return 0.5 * (u ^ v);
  // }


  friend void write_paraview_vtu (std::ostream&, const Interface3&);

private:
  Interface3(const Interface3 &); // pas de construction par copie
  void operator=(const Interface3 &);// pas affectation par copy

};







class RefPartition3 {

public:
  static const int nv = 4;                         // tetra
  typedef SortArray<Ubyte, nv> TetraIdx;           // the vertices of a triangle of the cut:

  typedef const TetraIdx* const_element_iterator;
  typedef       TetraIdx*       element_iterator;

  typedef Mesh3::Element Element;



  class InitializerCL
  {
  private:
    static int init_count_;

  public:
    InitializerCL ();
  };


private :


  TetraIdx tetras_[6];                               // at most 6 tetra
  element_iterator begin_;
  element_iterator end_;

  TetraIdx MakeTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3) const {
    return  TetraIdx( v0, v1, v2, v3); }

  // The sequences grow away from tetras_+3
  void AddTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3, int sign)
        { (sign == -1 ? *--begin_ : *end_++) = TetraIdx( v0, v1, v2, v3); }
  // e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
  void AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g0, Ubyte g1, int sign) {
    AddTetra( e0, e1, f1, g1, sign);
    AddTetra( e0, f0, f1, g1, sign);
    AddTetra( e0, f0, g0, g1, sign);
  }

  // If the intersection is quadrilateral, this returns the first uncut edge.
  // It is always one of the three first edges
  Ubyte first_uncut_edge (const SignPatternTrait3& cut) const {
    return cut[0] == 1 ? 0 : (cut[1] == 2 ? 1 : 2); }


  // 81 = 3^4 = all possible sign-patterns on the vertices
  static RefPartition3 instance_array_[81];

public :
  RefPartition3 () : begin_( tetras_ + 3), end_( tetras_) {}

  // Initialize with sign pattern on the vertices
  RefPartition3 (const SignPatternTrait3& cut) { assign( cut); }

  // Assign a sign pattern on the vertices; returns the value of is_uncut()
  bool assign (const SignPatternTrait3& cut);

  // True after assign(...)
  bool is_initialized () const { return begin_ < end_; }

  ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
  static inline const RefPartition3& instance (const byte   ls[4]) {
    RefPartition3& instance = instance_array_[instance_idx3 ( ls) + 40];
    if ( !instance.is_initialized()) {
      instance.assign( SignPatternTrait3( ls));
    }
    return instance;
  }
  static inline const RefPartition3& instance (const double ls[4])
  {
    byte ls_byte[4];
    std::transform( ls + 0, ls + 4, ls_byte + 0, sign);
    return instance( ls_byte);
  }
  ///@}

  // True, iff the partition has exactly one tetra
  bool is_uncut () const { return end_ == begin_ + 1; }
  bool is_cut () const { return !is_uncut(); }
  // Sign of the tetra, to which t points
  int whatSign (const_element_iterator t) const { return t < tetras_ + 3 ? -1 : 1; }


    ///@{ Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum
  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  {if(s == NoElement) return end_;
    else return s == PosElement ? tetras_ + 3 : begin_; }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return s == NegElement ? tetras_ + 3 : end_; }
    ///@}
  // Element make_subElement(int k);


};


class Partition3 {

public :
  typedef SortArray<Ubyte, 4> TetraIdx;           // the vertices of a triangle of the cut:
  typedef const TetraIdx* const_element_iterator;
  typedef       TetraIdx*       element_iterator;
  typedef Mesh3::Element Element;
  typedef Mesh3::Rd Rd;
  typedef CutData3 CutData;


  const RefPartition3& patch;
  const Element& T;
  double ls[4];
  const CutData* cutData = nullptr;


  Partition3(const Element& t, const double lls[4])
    : patch(RefPartition3::instance(lls)), T(t) {
    for(int i=0;i<4;++i) ls[i] = lls[i];
  }

  Partition3(const Element& t, const CutData& data)
    : patch(RefPartition3::instance(data.sign)), T(t), cutData(&data) {
    for(int i=0;i<4;++i) ls[i] = static_cast<double>(data.sign[i]);
  }
  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  { return patch.element_begin(s); }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return patch.element_end(s); }

  int whatSign (const_element_iterator t) const { return patch.whatSign(t); }
  bool is_cut () const { return patch.is_cut(); }
  bool isnot_cut () const { return patch.is_uncut(); }
  // ElementSignEnum what_part(const int the_domain) const {
  //   return (is_cut()) ? ( (the_domain == 0)? PosElement : NegElement  )
  //     : AllElement;
  // }

  ElementSignEnum what_part(const int the_domain) const {
    return (is_cut()) ? ( (the_domain == 0)? PosElement : NegElement  )
      : (((the_domain==0 && ls[0]>0) || (the_domain==1 && ls[0]<0))? AllElement : NoElement);
  }

  double getEdgeLength() const {
    return 1./6*(T.lenEdge(0) + T.lenEdge(1) + T.lenEdge(2)
                +T.lenEdge(3) + T.lenEdge(4) + T.lenEdge(5));
  }

  R mesure(const_element_iterator it) const {

    if( patch.is_uncut()) return T.mesure();
    Rd N[4];

    for(int i=0; i<4;++i) {
      Uint idx = (*it)[i];
      if(idx < 4) {N[i] = T[idx];
      }
      else {
        if(cutData) {
          N[i] = cutData->pointFromEdge(idx-4);
        }
        else{
          const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
          const R t = -ls[v0]/(ls[v1]-ls[v0]);
          N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
        }
      }
    }
    R3 AB(N[0],N[1]);
    R3 AC(N[0],N[2]);
    R3 AD(N[0],N[3]);

    return std::fabs(det(AB,AC,AD)/6.);
  }



  R mesure(const int domain) const {

    if( patch.is_uncut()) return T.mesure();

    R mes = 0;
    ElementSignEnum the_part = what_part(domain);

    for(const_element_iterator it = element_begin(the_part);
    it != element_end(the_part); ++it) {

      Rd N[4];

      for(int i=0; i<4;++i) {
        Uint idx = (*it)[i];
        if(idx < 4) {N[i] = T[idx];
        }
        else {
          const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
          const R t = -ls[v0]/(ls[v1]-ls[v0]);
          N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;

        }
      }
      R3 AB(N[0],N[1]);
      R3 AC(N[0],N[2]);
      R3 AD(N[0],N[3]);

      mes += std::fabs(det(AB,AC,AD)/6.);
    }
    return mes;
  }


  R mesureEdge(int e, int dom) const { assert(0);}
  Rd toEdge(int e, R2 ip, int dom) const { assert(0);}

  void plot(const_element_iterator it, std::ostream& out) const {
    Rd N[4];

    for(int i=0; i<4;++i) {
      Uint idx = (*it)[i];
      if(idx < 4) N[i] = R3::KHat[idx];
      else {
        const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        N[i] = (1-t) * ((Rd) R3::KHat[v0]) + t * ((Rd) R3::KHat[v1]) ;
      }
    }

    out << N[0] << "\n " << N[1] << "\n \n"
    << N[0] << "\n " << N[2] << "\n \n"
    << N[0] << "\n " << N[3] << "\n \n"
    << N[1] << "\n " << N[2] << "\n \n"
    << N[1] << "\n " << N[3] << "\n \n"
    << N[2] << "\n " << N[3] << "\n \n";
  }

  Rd toKref(const_element_iterator it, const Rd Phat) const {

    if( patch.is_uncut()) return Phat;
    Rd N[4];
    for(int i=0; i<4;++i) {
      Uint idx = (*it)[i];
      if(idx < 4) N[i] = R3::KHat[idx];
      else {
        const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        N[i] = (1-t) * ((Rd) R3::KHat[v0]) + t * ((Rd) R3::KHat[v1]) ;
      }
    }
    return  (1-Phat.x - Phat.y - Phat.z) * N[0] + Phat.x * N[1]
    + Phat.y * N[2] + Phat.z * N[3];

  }

  Rd toK(const_element_iterator it, const Rd Phat) const {

    if( patch.is_uncut()) return T(Phat);
    Rd N[4];
    for(int i=0; i<4;++i) {
      Uint idx = (*it)[i];
      if(idx < 4) N[i] = T[idx];
      else {
        if(cutData) {
          N[i] = cutData->pointFromEdge(idx-4);
        }
        else {
          const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
          const R t = -ls[v0]/(ls[v1]-ls[v0]);
          N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
        }
      }
    }
    return  (1-Phat.x - Phat.y - Phat.z) * N[0] + Phat.x * N[1]
    + Phat.y * N[2] + Phat.z * N[3];

  }


  Rd get_Vertex(const_element_iterator it, const int i) const {

    Uint idx = (*it)[i];
    if(idx < 4) return (Rd) T[idx];
    else {
      if(cutData) {
        return cutData->pointFromEdge(idx-4);
      }
      else {
        const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        return (1-t) * ((Rd)T[v0]) + t * ((Rd) T[v1]) ;
      }
    }
  }

  R get_LSvalue(const_element_iterator it, const int i) const {
    Uint idx = (*it)[i];
    if(idx < 4) return ls[idx];
    else {
      const Ubyte v0= Element::nvedge[idx - 4][0], v1= Element::nvedge[idx - 4][1];
      const R t = -ls[v0]/(ls[v1]-ls[v0]);
      return (1-t) * ls[v0] + t * ls[v1] ;
    }
  }
};


template<int d> struct TypeInterface {typedef Interface2 Interface;};
template<> struct TypeInterface<2> {typedef Interface2 Interface;};
template<> struct TypeInterface<3> {typedef Interface3 Interface;};





#endif
