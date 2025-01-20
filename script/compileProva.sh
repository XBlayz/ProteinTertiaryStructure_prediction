#!/bin/bash
cd ../test/
nasm -f elf64 prova.nasm
gcc -m64 -msse -mavx -O0 -no-pie sseutils64.o prova.o prova_assembly.c -o prova -lm
./prova