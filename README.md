# CutFEM-Library

to create theD Darcy library to reproduce results

1 - Build the library
"cd cpp; mkdir build; cd build; cmake ..; make -j4;"
2 - Go to the python folder
cd ../../python/darcy;


** set the variable in the cutFEMConfig.in
** need cmake
** set variable
  -  CMAKE_C_COMPILER
  -  CMAKE_CXX_COMPILER   (mpi if used)
** install lapacke
** install one library to solve linear system
  - MUMPS or UMFPACK
** If mumps, you will have to link to scotch or parmetis
** Set path in FindUMPFPACK or FindMUMPS and FindLAPACK


** create folder build
** run cmake .. in the licutfem/build folder
