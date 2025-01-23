#!/bin/bash
cd ../src/openMP_x64/
gcc -m64 -msse -mavx -O0 -no-pie pst_OpenMP_64c.c -o pst_OpenMP_64 -lm -fopenmp
