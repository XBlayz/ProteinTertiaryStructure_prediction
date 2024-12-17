#!/bin/bash
cd ../src/x64/;
./pst64_debug -seq ../../script/input/seq_256.ds2 -to 20 -alpha 1 -k 1 -sd 3;
gprof pst64 gmon.out > profile.txt