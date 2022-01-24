#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "problem.hpp"
int main(int argc, char** argv )
{

// std::string f = argv[1];

// Mesh2 Th(f.c_str());
Mesh3 Th(5, 10, 5, 0.,0.,0., 1., 2., 1.);
FESpace3 Vh(Th);
// std::cout << " P1 Space " << std::endl;
// Vh.info();
// FESpace2 Vhh(Th, DataFE<Mesh2>::P2);
// std::cout << " P2 Space " << std::endl;
// Vhh.info();
Paraview3 writer(Vh, "showMesh3Output.vtk");

// system("open -a /Applications/ParaView-5.7.0.app showMeshOutput.vtk" );

return 0;
};
