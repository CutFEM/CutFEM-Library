# About 
This library contains implementations of CutFEM algorithms on simple numerical examples. It is maintained and developed by Thomas Frachon, Sebastian Myrbäck and Erik Nilsson at the Department of Mathematics, KTH Royal Institute of Technology, Stockholm, Sweden.


# Installing the CutFEM Library

1. **Install dependencies**
    * For default configuration, you need to install Cmake and UMFPACK. You need to be able to run C++20, so make sure you have a compatible compiler.
    * To use e.g. MUMPS, you need to install mumps and mpi.

2. **Clone**
    ```
    git clone https://github.com/CutFEM/CutFEM-Library.git
    ```
3. **Configure**
    ```
    cd CutFEM-Library
    mkdir build
    cd build
    cmake ..
    ```
4. **Compile** \
    If you want to compile all files in the [cpp/example](cpp/example) folder, write ```make``` or ```make -jX``` for building with ```X``` processors. 

5. **Other settings** \
    If you want to configure the library with other settings, such as being able to use mumps for solving linear systems, you can change the settings in the [CMakeLists.txt](CMakeLists.txt) file.


# Reproducing data

### *A High-Order Conservative Cut Finite Element Method for Problems in Time-Dependent Domains*, 2024, S. Myrbäck, S. Zahedi. [Link here](https://www.sciencedirect.com/science/article/pii/S0045782524005012)

This paper presents a conservative space-time CutFEM for solving convection-diffusion equations in bulk, or bulk-surface domains.

The files are located in [cpp/example/convection_diffusion](cpp/example/convection_diffusion/) and consist of two files: [bulk.cpp](cpp/example/convection_diffusion/bulk.cpp), and [coupled.cpp](cpp/example/convection_diffusion/coupled.cpp).

To run them you first need to build them, by writing e.g. ```make -j4 bulk```, unless you've already built them when installing the library above.

To execute and run the files, there are several options/parameters to choose from:

#### Bulk

1. **Example**: ```circle/kite```
2. **Method**: ```conservative/non_conservative```
3. **Polynomial order in space**: ```1/2/3```
4. **Polynomial order in time**: ```0/1/2/3```
5. **Ghost-penalty stabilization technique**: ```fullstab/macro```

These options are set by defining preprocessor variables before the main function, e.g. ```#define circle```. Other user variables, such as stabilization constant, macroelement parameter, and number of quadrature nodes in time and space, are changed inside the main function. The code is then executed as:

```
./bin/bulk
```


#### Coupled

1. **Method**: ```conservative/non_conservative```
2. **Polynomial order in space**: ```1/2```
3. **Polynomial order in time**: ```0/1/2/3```
4. **Ghost-penalty stabilization technique**: ```fullstab/macro```

Stabilization constants, macroelement parameters, number of quadrature nodes in time and space, etc, are changed inside the main function. The code is then executed as:

```
./bin/coupled
```
Convergence history is printed both in the terminal and in the files "output_bulk.dat" and "output_coupled.dat". 


### *A divergence preserving cut finite element method for Darcy flow*, 2024, Thomas Frachon, Erik Nilsson, Sara Zahedi
To create the Darcy problem to reproduce results using Python:
1 - All options can be OFF apart from CUTFEM_BUILD_PYTHON_WRAPPER
2. "cmake ..; make -j4;"
2 - Go to the python folder
"cd ../../python/darcy;"
3 - In the darcy_wrapper.py file one must link to the library. Verify that the good extension is used.
4 - run darcy.py


Note: It is important to turn off the options for finding libraries if they are not installed, otherwise the compilation will not succeed.

Note: When compiling the python library, the MPI option should be turned off.

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



