#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "problem.hpp"
int main(int argc, char** argv )
{

std::string f = argv[1];

Mesh2 Th(f.c_str());
FESpace2 Vh(Th);
std::cout << " P1 Space " << std::endl;
Vh.info();
FESpace2 Vhh(Th, DataFE<Mesh2>::P2);
std::cout << " P2 Space " << std::endl;
Vhh.info();
Paraview2 writer(Vh, "showMeshOutput.vtk");

system("open -a /Applications/ParaView-5.7.0.app showMeshOutput.vtk" );

return 0;
};
