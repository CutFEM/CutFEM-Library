#include "paraview.hpp"

template<> void Paraview<Mesh2>::ParaviewMesh::build(const ActiveMesh<Mesh>&  Vh) {
  return buildNoCut(Vh);
}
template<> void Paraview<MeshQuad2>::ParaviewMesh::build(const ActiveMesh<Mesh>&  Vh) {
  return buildNoCut(Vh);
}
template<> void Paraview<MeshHexa>::ParaviewMesh::build(const ActiveMesh<Mesh>&  Vh) {
  return buildCut(Vh);
}
template<> void Paraview<Mesh3>::ParaviewMesh::build(const ActiveMesh<Mesh>& Vh) {
  return buildCut(Vh);
}

// template<> void Paraview<Mesh2>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
//   return buildNoCut(Vh, levelSet);
// }
// template<> void Paraview<MeshQuad2>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
//   return buildNoCut(Vh, levelSet);
// }
// template<> void Paraview<MeshHexa>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
//   return buildCut(Vh, levelSet);
// }
// template<> void Paraview<Mesh3>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
//   return buildCut(Vh, levelSet);
// }
