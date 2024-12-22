#!/bin/bash
cd ../test/
if [ ! -f v64 ]; then
    gcc verify_output64.c -o v64
fi
./v64