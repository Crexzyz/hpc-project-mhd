#!/bin/bash

for i in {1..20}
do
    mpirun -n 4 ./build/mpi.out >> build/mpi.txt
done