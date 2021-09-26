#!/bin/bash

cd ../common/build/;
make -j 4;
cd ../../FESpace/build/;
make -j 4;
cd ../../problem/build/;
make -j 4;
cd ../../build/;
make -j 4 main;
echo

./bin/main
