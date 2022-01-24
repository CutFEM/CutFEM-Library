#ifndef INTERFACE2DN_HPP_
#define INTERFACE2DN_HPP_

#include "GenericInterface.hpp"

#include "Mesh2dn.hpp"
// #include "marker.hpp"

class Interface2 : public GenericInterface<Mesh2>
{
  typedef CutData2 CutData;
public:
  Interface2() : GenericInterface<Mesh2>() {}
  Interface2(const Mesh & MM) : GenericInterface<Mesh2>(MM) {}
  Interface2(const Mesh & MM, const KN<double>& ls, int label = 0) : GenericInterface<Mesh2>(MM, ls, label) {}

  R2 mapToFace(const FaceIdx& f, const R1 PHat ) const  {
    return (1-PHat.x)*vertices_[f[0]] + PHat.x*vertices_[f[1]];
  }
  R2 computeDx(const FaceIdx& f) const {
    return Rd(vertices_[f[0]],vertices_[f[1]]);
  }
  virtual CutData getCutData(const int k) const {
    assert(backMesh);
    const Element& K((*backMesh)[k]);

    byte loc_sign[3];
    R2 point[2];
    int array_edge[3] = {-1,-1,-1};
    for(int i=0;i<Element::nv;++i) {

      loc_sign[i] = ls_sign((*backMesh)(k, i));
    }

    const typename Mesh2::RefPatch& cut =  RefPatch::instance( loc_sign);

    if(isCut(k)){
      int size = cut.size();
      int i0 = 0;
      int idPoint = 0;
      for (typename RefPatch::const_face_iterator it= cut.face_begin(), end= cut.face_end();
      it != end; ++it, ++i0) {

        int idxEdge[Mesh::nva];
        for (Uint j= 0; j < Mesh::nva; ++j) {
          idxEdge[j] = (*it)[j] ;
          if(idxEdge[j]>2)  idxEdge[j] -= 3;
          else assert(0); // no point on node
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

private:
  Interface2(const Interface2 &); // pas de construction par copie
  void operator=(const Interface2 &);// pas affectation par copy
};




class RefPartition2 {

public:
  static const int nv = 3;                         // tetra
  typedef SortArray<Ubyte, nv> TriaIdx;           // the vertices of a triangle of the cut:

  typedef const TriaIdx* const_element_iterator;
  typedef       TriaIdx*       element_iterator;

  typedef Mesh2::Element Element;



  class InitializerCL
  {
  private:
    static int init_count_;

  public:
    InitializerCL ();
  };


private :


  TriaIdx triangles_[4];                               // at most 3 triangles
  element_iterator begin_;
  element_iterator end_;

  TriaIdx MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const {
    return  TriaIdx( v0, v1, v2); }

  // The sequences grow away from tetras_+3
  void AddTriangle (Ubyte v0, Ubyte v1, Ubyte v2, int sign)
        { (sign == -1 ? *--begin_ : *end_++) = TriaIdx( v0, v1, v2); }
  // e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
  void AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, int sign) {
    AddTriangle( e0, e1, f0, sign);
    AddTriangle( e1, f0, f1, sign);
  }

  // If the intersection is quadrilateral, this returns the first uncut edge.
  // It is always one of the three first edges
  Ubyte first_uncut_edge (const SignPatternTrait2& cut) const {
    return cut[0] == 1 ? 0 : (cut[1] == 2 ? 1 : 2); }


  // 81 = 3^4 = all possible sign-patterns on the vertices
  static RefPartition2 instance_array_[27];

public :
  RefPartition2 () : begin_( triangles_ + 2), end_( triangles_) {}

  // Initialize with sign pattern on the vertices
  RefPartition2 (const SignPatternTrait2& cut) { assign( cut); }

  // Assign a sign pattern on the vertices; returns the value of is_uncut()
  bool assign (const SignPatternTrait2& cut);

  // True after assign(...)
  bool is_initialized () const { return begin_ < end_; }

  ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
  static inline const RefPartition2& instance (const byte   ls[3]) {
    RefPartition2& instance = instance_array_[instance_idx2 (ls) + 13];
    if ( !instance.is_initialized()) {
      instance.assign( SignPatternTrait2( ls));
    }
    return instance;
  }
  static inline const RefPartition2& instance (const double ls[3])
  {
    byte ls_byte[3];
    std::transform( ls + 0, ls + 3, ls_byte + 0, sign);
    return instance( ls_byte);
  }
  ///@}

  // True, iff the partition has exactly one tetra
  bool is_uncut () const { return end_ == begin_ + 1; }
  bool is_cut () const { return !is_uncut(); }
  // Sign of the tetra, to which t points
  int whatSign (const_element_iterator t) const { return t < triangles_ + 2 ? -1 : 1; }


    ///@{ Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum
  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  { if(s == NoElement) return end_;
    else return s == PosElement ? triangles_ + 2 : begin_; }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return s == NegElement ? triangles_ + 2 : end_; }
    ///@}
  // Element make_subElement(int k);


};



class Partition2 {

public :
  typedef SortArray<Ubyte, 3> TriaIdx;           // the vertices of a triangle of the cut:
  typedef const TriaIdx* const_element_iterator;
  typedef       TriaIdx*       element_iterator;
  typedef Mesh2::Element Element;
  typedef Mesh2::Rd Rd;
  typedef CutData2 CutData;


  const RefPartition2& patch;
  const Element& T;
  double ls[3];
  const CutData* cutData = nullptr;


  Partition2(const Element& t, const double lls[3])
    : patch(RefPartition2::instance(lls)), T(t) {
    for(int i=0;i<3;++i) ls[i] = lls[i];
  }
  Partition2(const Element& t, const CutData& data)
    : patch(RefPartition2::instance(data.sign)), T(t), cutData(&data) {
    for(int i=0;i<3;++i) ls[i] = static_cast<double>(data.sign[i]);
  }


  const_element_iterator element_begin (ElementSignEnum s = AllElement) const
  { return patch.element_begin(s); }
  const_element_iterator element_end (ElementSignEnum s = AllElement) const
  { return patch.element_end(s); }

  int whatSign (const_element_iterator t) const { return patch.whatSign(t); }
  bool is_cut () const { return patch.is_cut(); }
  bool isnot_cut () const { return patch.is_uncut(); }
  // ElementSignEnum what_part(const int the_domain) const {
  //   return (is_cut()) ? ( (the_domain == 0)? PosElement : NegElement  ) :AllElement;
  // }
  ElementSignEnum what_part(const int the_domain) const {
    if(the_domain == -1) return AllElement;
    return (is_cut()) ? ( (the_domain == 0)? PosElement : NegElement  )
      : (((the_domain==0 && ls[0]>0) || (the_domain==1 && ls[0]<0))? AllElement : NoElement);
  }

  double getEdgeLength() const {
    return 1./3*(T.lenEdge(0)+ T.lenEdge(1)+T.lenEdge(2));
  }

  R mesure(const_element_iterator it) const {

    if( patch.is_uncut()) return T.mesure();
    Rd N[3];

    for(int i=0; i<3;++i) {
      Uint idx = (*it)[i];
      if(idx < 3) {N[i] = T[idx];
      }
      else {
        if(cutData) {
          N[i] = cutData->pointFromEdge(idx-3);
        }
        else {
          const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
          const R t = -ls[v0]/(ls[v1]-ls[v0]);
          N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
        }
      }
    }

    return std::fabs(det(N[0],N[1],N[2]))*0.5;
  }

  R mesure(const int domain) const {

    if( patch.is_uncut()) return T.mesure();

    R mes = 0;
    ElementSignEnum the_part = what_part(domain);

    for(const_element_iterator it = element_begin(the_part);
    it != element_end(the_part); ++it) {

      Rd N[3];
      for(int i=0; i<3;++i) {
        Uint idx = (*it)[i];
        if(idx < 3) {N[i] = T[idx];}
        else {
          if(cutData) {
            N[i] = cutData->pointFromEdge(idx-3);
          }
          else {
            const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
            const R t = -ls[v0]/(ls[v1]-ls[v0]);
            N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
          }
        }
      }
    mes += std::fabs(det(N[0],N[1],N[2]))*0.5;
  }
  return mes;
}

  R mesureEdge(int e, int dom) const {
    assert(cutData);
    if( patch.is_uncut()) return T.lenEdge(e);
    int s = (dom == 0)? 1 : -1;

    const Ubyte v0= Element::nvedge[e][0], v1= Element::nvedge[e][1];
    // const R t = -ls[v0]/(ls[v1]-ls[v0]);
    // R2 A = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
    R2 A = cutData->pointFromEdge(e);
    R2 B = (s * ls[v0] > 0) ? T[v0] : T[v1];
    R2 E(A,B);
    return E.norm();
  }

  Rd toEdge(int e, R1 ip, int dom) const {
    assert(cutData);
    if( patch.is_uncut()) assert(0);
    int s = (dom == 0)? 1 : -1;

    const Ubyte v0= Element::nvedge[e][0], v1= Element::nvedge[e][1];
    R2 A = cutData->pointFromEdge(e);
    R2 B = (s * ls[v0] > 0) ? T[v0] : T[v1];
    R2 E = (1-ip.x) * A + ip.x * B ;

    return E;
  }

  // Rd toKref(const_element_iterator it, const Rd Phat) const {
  //
  //   if( patch.is_uncut()) return Phat;
  //   Rd N[3];
  //   for(int i=0; i<3;++i) {
  //     Uint idx = (*it)[i];
  //     if(idx < 3) N[i] = R2::KHat[idx];
  //     else {
  //       const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
  //       const R t = -ls[v0]/(ls[v1]-ls[v0]);
  //       N[i] = (1-t) * ((Rd) R2::KHat[v0]) + t * ((Rd) R2::KHat[v1]) ;
  //     }
  //   }
  //   return  (1-Phat.x - Phat.y) * N[0] + Phat.x * N[1] + Phat.y * N[2];
  // }


  Rd toK(const_element_iterator it, const Rd Phat) const {

    if( patch.is_uncut()) return T(Phat);
    Rd N[3];
    for(int i=0; i<3;++i) {
      Uint idx = (*it)[i];
      if(idx < 3) N[i] = T[idx];
      else {
        if(cutData) {
          N[i] = cutData->pointFromEdge(idx-3);
        }
        else {
          const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
          const R t = -ls[v0]/(ls[v1]-ls[v0]);
          N[i] = (1-t) * ((Rd) T[v0]) + t * ((Rd) T[v1]) ;
        }
      }
    }
    return  (1-Phat.x - Phat.y) * N[0] + Phat.x * N[1]
    + Phat.y * N[2];

  }

  Rd get_Vertex(const_element_iterator it, const int i) const {

    Uint idx = (*it)[i];
    if(idx < 3) return (Rd) T[idx];
    else {
      if(cutData) {
        return cutData->pointFromEdge(idx-3);
      }
      else {
        const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
        const R t = -ls[v0]/(ls[v1]-ls[v0]);
        return (1-t) * ((Rd)T[v0]) + t * ((Rd) T[v1]) ;
      }
    }
  }

  R get_LSvalue(const_element_iterator it, const int i) const {
    Uint idx = (*it)[i];
    if(idx < 3) return ls[idx];
    else {
      const Ubyte v0= Element::nvedge[idx - 3][0], v1= Element::nvedge[idx - 3][1];
      const R t = -ls[v0]/(ls[v1]-ls[v0]);
      return (1-t) * ls[v0] + t * ls[v1] ;
    }
  }


};


#endif
