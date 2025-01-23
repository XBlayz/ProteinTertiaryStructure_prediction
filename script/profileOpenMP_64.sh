#!/bin/bash
cd ../src/openMP_x64/
./pst_OpenMP_64_debug -seq ../../script/input/seq_256.ds2 -to 20 -alpha 1 -k 1 -sd 3
gprof pst_OpenMP_64_debug gmon.out > profile.txt