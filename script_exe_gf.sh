#!/bin/sh
#
         module load mkl intel gcc/8.2
         gfortran -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl qbithm.f
         sbatch job-qbithm
