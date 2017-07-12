#/bin/bash

g++ -c prepareKernel_2DSqu.cpp -o out.o
g++ out.o -L/opt/intel/Compiler/11.0/074/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -o a.out
