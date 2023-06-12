#!/bin/bash

for i in {1..20}
do
    ./build/serial.out >> build/serial.txt
done