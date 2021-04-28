#!/bin/bash
#PJM -N NICAM-DIFFUSION
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

TMPDIR=${HOME}/tmp/check_NICAM.diffusion
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *

SRC_DIR=${HOME}/fs2020_kernels/src/NICAM.diffusion.20170220
cp -rp ${SRC_DIR}/* ./

OPTIMIZE="-Kfast,openmp,ocl -falign-loops -Icommon/include "
FOPTIMIZE="${OPTIMIZE} -Cpp -fw -DUSE_FAPP -DSINGLE "


frt ${FOPTIMIZE} postK_nicam_diffusion.f90
#	fcc -c ${COPTIMIZE} common/src/report.c
#	frt ${COPTIMIZE} adventure.o report.o

export OMP_NUM_THREADS=12
export FLIB_TRACEBACK_MEM_SIZE=128
time ./a.out
