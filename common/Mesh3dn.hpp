#ifndef MESH3DN_HPP_
#define MESH3DN_HPP_

#include <cstdlib>
#include "dataStruct3D.hpp"
#include "GenericMesh.hpp"
using namespace std;


class Mesh3 : public GenericMesh<Tet,Triangle3,Vertex3> {

public:

  Mesh3(){}
  Mesh3(const string);
  Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb);
  Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb);  // surface mesh
  Mesh3(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz);

private:
  void readmsh(ifstream & f);


  Mesh3(const Mesh3 &); // pas de construction par copie
  void operator=(const Mesh3 &);// pas affectation par copy
};


class MeshHexa : public GenericMesh<Hexa,Quad3,Vertex3>
{
public:

  MeshHexa(int nx, int ny, int nz, R orx, R ory,R orz, R lx, R ly,R lz);  // build structured mesh
private:
  MeshHexa(const MeshHexa &);                             // no copy constructor
  void operator=(const MeshHexa &);                    // no copy allowed
};



typedef Mesh3     MeshT3;
typedef MeshHexa  MeshQ3;

#endif
