
#ifndef MESH1DN_HPP_
#define MESH1DN_HPP_

#include "GenericMesh.hpp"
#include <cstdlib>

typedef GenericVertex<R1> Vertex1;

// class SignPatternTrait1;
class RefPatch2;
// class RefPartition1;
// class Partition1;

class Mesh1 : public GenericMesh<Seg1, BoundaryPoint1, Vertex1> {
 public:
   // typedef SignPatternTrait1 SignPattern;
   typedef RefPatch2 RefPatch;
   // typedef RefPartition1 RefPartition;
   // typedef Partition1 Partition;

   Mesh1(int nx, R orx, R lx); // build structured mesh
   Mesh1(const char *);        //
   const Element *Find(R1 P, R1 &Phat, bool &outside,
                       const Element *tstart) const;

 private:
   int load(const std::string &filename);
   Mesh1(const Mesh1 &);          // pas de construction par copie
   void operator=(const Mesh1 &); // pas affectation par copy
};

#endif
