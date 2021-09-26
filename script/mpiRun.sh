#!/bin/bash

echo -n "How many processors ?"
read  nproc
cd common/build/;
make -j 4;
cd ../../FESpace/build/;
make -j 4;
cd ../../problem/build/;
make -j 4;
cd ../../build/;
make -j 4 main;
echo

mpirun -np $nproc ./bin/main
