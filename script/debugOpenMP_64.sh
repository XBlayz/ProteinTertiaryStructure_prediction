#!/bin/bash
cd ../src/openMP_x64/
for f in $(ls *64.nasm); do
	nasm -f elf64 $f
done
gcc -Wall -pg -m64 -msse -mavx -O0 -no-pie *64.o pst_OpenMP_64c.c -o pst_OpenMP_64_debug -lm -fopenmp
