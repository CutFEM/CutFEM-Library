# CutFEM-Library

To create the Darcy problem to reproduce results using Python:

1 - Build the library
"cd cpp; mkdir build; cd build; cmake ..; make -j4;"
2 - Go to the python folder
"cd ../../python/darcy;"
3 - run darcy.py
"python3 darcy.py"


Note : It is important to turn off the options for finding libraries if they are not installed, otherwise the compilation will not succeed.

Note : When compiling the python library, the MPI option should be turned off.

Note : If you want to run the tests, you have to first download the submodule Catch2
"git submodule update --init --cpp/test/Catch2"

On another hand, to update all modules (test, solver etc) one can do
"git submodule update --init --recursive"


TO INSTALL UMFPACK
1) add the module SuiteSparse
"git submodule update --init --SuiteSparse"
2) Go to UMFPACK folder
"make local; make install;
3) Fixe the cmake/FindUMFPACK.cmake


Note : 

