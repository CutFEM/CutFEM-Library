#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "baseProblem.hpp"
#include "paraview.hpp"

typedef FunFEM<Mesh2> Fun_h;
typedef Mesh2 Mesh;
typedef ActiveMeshT2 CutMesh;
typedef FESpace2   Space;
typedef CutFESpaceT2 CutSpace;

R shift = 0.;
R fun_levelSet(const R2 P, const int i) {
  return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - 0.5;
}

int main(int argc, char** argv )
{


  // Mesh2 Kh("../mesh/FittedMesh.msh");
  Mesh Kh("../mesh/square2.msh");
  Mesh Th(100, 100, -1., -1., 2., 2.);
  Space Lh(Kh, DataFE<Mesh2>::P1);
  Space Lh2(Th, DataFE<Mesh2>::P1);
  Fun_h levelSet(Lh, fun_levelSet);
  Fun_h levelSet2(Lh2, fun_levelSet);

  Paraview<Mesh> writer(Kh, "staticDropMesh.vtk");
  writer.add(levelSet, "levelSet" , 0, 1);
  Paraview<Mesh> writer2(Th, "levelSet.vtk");
  writer2.add(levelSet2, "levelSet" , 0, 1);



return 0;
};
