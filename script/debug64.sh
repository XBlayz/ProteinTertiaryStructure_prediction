#!/bin/bash
cd ../src/x64/
for f in $(ls *64.nasm); do
	nasm -f elf64 $f
done
gcc -Wall -pg -m64 -msse -mavx -O0 -no-pie *64.o pst64c.c -o pst64_debug -lm
