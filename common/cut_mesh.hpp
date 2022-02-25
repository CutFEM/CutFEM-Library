#ifndef _CUT_MESH_HPP
#define _CUT_MESH_HPP



#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"


// 1  2 3 4    Omega
// 1     3     Omega0
//    2    4   Omega1
template<typename T,typename B,typename V>
class Cut_Mesh : public GenericMesh<T,B,V> {

  typedef GenericMesh<T,B,V> Mesh;
  typedef T Element;

  const Mesh& Th;

  std::vector<std::map<int, int>> idx_in_background_mesh;   // [domain](idxK_cutMesh) -> idxK_backMesh
  std::vector<std::map<int,int>>  idx_from_background_mesh; // [domain](idxK_backMesh) -> idxK_cutMesh

  std::vector<int> nt_paraview;   // nb of element for mesh in paraview

public:
  Cut_Mesh(const Mesh& th) : Mesh() , Th(th) { }

  // Give the background mesh and a sign Function defined on the mesh nodes
  // Will create 2 subdomains
  Cut_Mesh(const Mesh& th, const Rn& level_set) : Mesh() , Th(th) ,
  idx_in_background_mesh(2), idx_from_background_mesh(2), nt_paraview(2) {
    this->add(level_set,  1, 0);
    this->add(level_set, -1, 1);
  }


private:
  void add(const Rn& level_set, int sign_domain, int dom_id);

  bool check_exist(int k, int dom) const {
    const auto it = idx_from_background_mesh[dom].find(k);
    if(it == idx_from_background_mesh[dom].end()) return false;
    else return true;
  }




  // void print
  // virtual DataFENodeDF GenericMesh<T,B,V>::BuildDFNumbering(


};


typedef Cut_Mesh<Triangle2,BoundaryEdge2,Vertex2> Cut_Mesh2;
typedef Cut_Mesh<Tet,Triangle3,Vertex3> Cut_Mesh3;












#endif
