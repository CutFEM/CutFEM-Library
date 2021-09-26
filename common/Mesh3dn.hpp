#ifndef MESH3DN_HPP_
#define MESH3DN_HPP_

#include <cstdlib>
#include "dataStruct3D.hpp"
#include "GenericMesh.hpp"
using namespace std;

template<typename M> class GenericInterface;
class Interface3;
class TimeInterface3;
class SignPatternTrait3;
class RefPatch3;
class RefPartition3;
class Partition3;

class Mesh3 : public GenericMesh<Tet,Triangle3,Vertex3> {

public:
  typedef SignPatternTrait3 SignPattern;
  typedef RefPatch3 RefPatch;
  typedef RefPartition3 RefPartition;
  typedef Partition3 Partition;
  typedef GenericInterface<Mesh3> Interface;

  // const Mesh3 * backMesh = nullptr;

  Mesh3(){}
  Mesh3(const string);
  Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb);
  Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb);  // surface mesh
  Mesh3(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz);
  Mesh3(const Interface& gamma);
  Mesh3(TimeInterface3& gamma);
  Mesh3(const Mesh3&, std::string whatToDo);





  // virtual Uint idxGlobalVertexBackMesh(const Uint i) const {
  //   Uint j =  idxVertexInBackMesh(i);
  //   return (backMesh) ? backMesh->idxGlobalVertex(j) : j;
  // }



  friend void write_paraview_vtu (std::ostream&, const Mesh3&);

private:
  int load(const string & filename);
  void readmsh(ifstream & f);


  Mesh3(const Mesh3 &); // pas de construction par copie
  void operator=(const Mesh3 &);// pas affectation par copy
};



/*
 *   Represents the reference tetra, which is cut by a linear level set function ls.
 *   The values of the latter are prescribed on the vertices.
 */

class SignPatternTrait3 : public GSignPatternTrait<Tet>
{

public :
  bool is_3d () const { return num_root_vert_ == 4; }
  bool is_2d () const { return num_root_ > 2; }  ///< True, iff the intersection has positive area.
  friend std::ostream& operator<< (std::ostream&, const SignPatternTrait3&); ///< Debug-output to a stream (dumps all members)

  SignPatternTrait3() : GSignPatternTrait<Tet>() {}
  SignPatternTrait3(const byte    ls[4]) : GSignPatternTrait<Tet>(ls) {}
  SignPatternTrait3(const double  ls[4]) : GSignPatternTrait<Tet>(ls) {}

};


/*
 *  Return a signed array-index for the possible 3^4 sign-patterns on the vertices of a tetra.
 *  The index ranges from [-40..40].
 */
inline static byte instance_idx3 (const byte ls[4])
{
  return  27*ls[0] + 9*ls[1] + 3*ls[2] + ls[3];
}

inline static Ubyte instance_idx3 (const double ls[4])
{
  return  27*sign( ls[0]) + 9*sign( ls[1]) + 3*sign( ls[2]) + sign(ls[3]);
}



/*
 * The triangle of the intersection of the reference-tet with a linear levelset-function.
 * The class memoizes used sign-patterns
 */
class RefPatch3
{
public:
  static const int nv = 3;
  typedef SortArray<Ubyte, nv> FaceIdx;           // the vertices of a triangle of the cut:
                                                     //the tetra's vertices are denoted by 0..3,
                                                     //the edge-cuts by edge-num + 4
  typedef const FaceIdx* const_face_iterator;
  typedef       FaceIdx*       face_iterator;

  class InitializerCL
  {
  private:
    static int init_count_;

  public:
    InitializerCL ();
  };

private:
  FaceIdx face_[2];      ///< at most two triangles
  Ubyte size_;                 ///< number of triangles
  Ubyte is_boundary_face_; ///< true if the triangle is one of the tetra's faces.

  Ubyte num_elements (const SignPatternTrait3& cut) const {
    return cut.is_2d() ? cut.num_cut_simplexes() - 2 : 0; }
  FaceIdx MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const {
    return FaceIdx( v0, v1,v2); }

  static RefPatch3 instance_array_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices


public:
  RefPatch3 () : size_( static_cast<Ubyte>( -1)), is_boundary_face_( 0) {}

  ///< Initialize with sign pattern on the vertices
  RefPatch3 (const SignPatternTrait3& cut) { assign( cut); }

  ///< Assign a sign pattern on the vertices; returns the value of empty()
  bool assign (const SignPatternTrait3& cut);

  ///< True after assign(...)
  bool  is_initialized () const { return size_ <=2; }


  ///@{ Recommended access to the triangles for a given sign-pattern; memoizes the result.
  static inline const RefPatch3& instance (const byte   ls[4]) {
    RefPatch3& instance = instance_array_[instance_idx3 ( ls) + 40];

    if ( !instance.is_initialized()) {
      instance.assign( SignPatternTrait3( ls));
    }
    return instance;
  }
  static inline const RefPatch3& instance (const double ls[4])
  {
    byte ls_byte[4];
    std::transform( ls + 0, ls + 4, ls_byte + 0, sign);
    return instance( ls_byte);
  }
  ///@}

  ///< true, iff the area of the intersection is 0.
  bool  empty () const { return size_ == 0; }

  ///< Number of triangles, 0, 1, or 2
  size_t size () const { return size_; }

  const_face_iterator face_begin () const { return face_; }
  const_face_iterator face_end   () const { return face_ + size_; }

};



#endif
