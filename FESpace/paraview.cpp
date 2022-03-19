#include "paraview.hpp"

template<> void ParaviewCut<Mesh2>::ParaviewMesh::build(const Cut_Mesh<Mesh>&  Vh) {
  return buildNoCut(Vh);
}
template<> void ParaviewCut<MeshQuad2>::ParaviewMesh::build(const Cut_Mesh<Mesh>&  Vh) {
  return buildNoCut(Vh);
}
template<> void ParaviewCut<MeshHexa>::ParaviewMesh::build(const Cut_Mesh<Mesh>&  Vh) {
  return buildCut(Vh);
}
template<> void ParaviewCut<Mesh3>::ParaviewMesh::build(const Cut_Mesh<Mesh>& Vh) {
  return buildCut(Vh);
}

template<> void Paraview<Mesh2>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
  return buildNoCut(Vh, levelSet);
}
template<> void Paraview<MeshQuad2>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
  return buildNoCut(Vh, levelSet);
}
template<> void Paraview<MeshHexa>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
  return buildCut(Vh, levelSet);
}
template<> void Paraview<Mesh3>::ParaviewMesh::build(const FESpace & Vh, Fun_h* levelSet) {
  return buildCut(Vh, levelSet);
}
