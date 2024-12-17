#!/bin/bash
cd ../src/x86/;
for f in $(ls *32.nasm); do
	nasm -f elf32 $f;
done;
gcc -m32 -msse -O0 -no-pie *32.o pst32c.c -o pst32 -lm
