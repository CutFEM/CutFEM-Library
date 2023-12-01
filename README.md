# Introduction 
This library aim to implement and test CutFEM algorithms on simple examples.


# CutFEM-Library
1. add a build folder in the root folder.
2. add the library that you want to link.
    If you run sequential, then just enable USE_UMFPACK
3. Set the option of the building procedure
4. "cd build; cmake ..; make -j4;"

# Running the example using Python
To create the Darcy problem to reproduce results using Python:
1 - All option can be OFF appart from CUTFEM_BUILD_PYTHON_WRAPPER
2. "cmake ..; make -j4;"
2 - Go to the python folder
"cd ../../python/darcy;"
3 - In the darcy_wrapper.py file one must link to the library. Verify that the good extension is used.
4 - run darcy.py


Note : It is important to turn off the options for finding libraries if they are not installed, otherwise the compilation will not succeed.

Note : When compiling the python library, the MPI option should be turned off.

On another hand, to update all modules (test, solver etc) one can do
"git submodule update --init --recursive"


# TO INSTALL A SUBMODULE
git submodule update --init --recusive submodule_name

# TO INSTALL UMFPACK
1. add the module SuiteSparse
"git submodule update --init --SuiteSparse"
2. Go to UMFPACK folder
"make local; make install;
3. Fixe the cmake/FindUMFPACK.cmake



