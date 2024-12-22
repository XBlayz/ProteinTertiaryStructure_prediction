#!/bin/bash
cd ../test/
if [ ! -f v32 ]; then
    gcc verify_output32.c -o v32
fi
./v32