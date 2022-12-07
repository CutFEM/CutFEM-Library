# CutFEM-Library

To create theD Darcy library to reproduce results

1 - Build the library
"cd cpp; mkdir build; cd build; cmake ..; make -j4;"
2 - Go to the python folder
"cd ../../python/darcy;"
3 - run darcy.py
"python3 darcy.py"


Note : It is important to turn off the options for finding libraries if they are not installed, otherwise the compilation will not succeed.

Note : If you want to run the tests, you have to first download the submodule Catch2
"git submodule update --init --recursive"


