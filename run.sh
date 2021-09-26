#!/bin/bash

echo -n "How many processors ?"
read  nproc
cd common/build/;
make -j 4;
cd ../../parallel/build/;
make -j 4;
cd ../../FESpace/build/;
make -j 4;
cd ../../solver/build/;
make -j 4;
cd ../../problem/build/;
make -j 4;
cd ../../build/;
make -j 4;
mpirun -np $nproc ./bin/main
