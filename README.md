# cutfem

to install 

** set the variable in the cutFEMConfig.in
** need cmake
** set variable
  -  CMAKE_C_COMPILER
  -  CMAKE_CXX_COMPILER   (mpi if used)
** install lapacke
** install one library to solve linear system
  - MUMPS or UMPPACK

** Set path in FindUMPFPACK or FindMUMPS and FindLAPACK


** create folder build
** run cmake .. in the licutfem/build folder
