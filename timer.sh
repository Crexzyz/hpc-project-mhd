#!/bin/bash

for i in {1..20}
do
    ./build/omp.out >> build/omp.txt
done