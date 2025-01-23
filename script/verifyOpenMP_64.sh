#!/bin/bash
cd ../test/
if [ ! -f vOpenMP_64 ]; then
    gcc verify_output_OpenMP_64.c -o vOpenMP_64
fi
./vOpenMP_64