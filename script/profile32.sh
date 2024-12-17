#!/bin/bash
cd ../src/x86/;
./pst32_debug -seq ../../script/input/seq_256.ds2 -to 20 -alpha 1 -k 1 -sd 3;
gprof pst32 gmon.out > profile.txt