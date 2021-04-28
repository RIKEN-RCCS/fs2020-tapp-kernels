#!/bin/bash
#PJM -N STREAMLIKE-PATTERN1
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:10:00"
#PJM --mpi "max-proc-per-node=1"
#PJM --rsc-list "node=1"
#PJM -j
#PJM -S

module list
set -x
date
hostname

TMPDIR=${HOME}/tmp/check_NICAM.st_pattern1
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *

SRC_DIR=${HOME}/fs2020_kernels/src/streamlike_pattern1
cp -rp ${SRC_DIR}/* ./

OPTIMIZE="-Kfast,openmp,ocl -falign-loops -Icommon/include "
FOPTIMIZE="${OPTIMIZE} -Cpp -fw -DSINGLE "


frt -c ${FOPTIMIZE} mod_precision.f90
frt -c ${FOPTIMIZE} mod_streamlike.f90
frt -c ${FOPTIMIZE} main.f90
#	fcc -c ${COPTIMIZE} common/src/report.c
frt ${FOPTIMIZE} *.o

export OMP_NUM_THREADS=12
export FLIB_TRACEBACK_MEM_SIZE=128
time ./a.out
