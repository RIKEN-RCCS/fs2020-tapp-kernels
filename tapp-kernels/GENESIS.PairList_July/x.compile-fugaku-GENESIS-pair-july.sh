#!/bin/bash
#PJM -N GENE-PAIRLIST-JULY
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

TMPDIR=${HOME}/tmp/check_GENESIS-pairlist-july
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *

SRC_DIR=${HOME}/fs2020_kernels/src/GENESIS.PairList_July
cp -rp ${SRC_DIR}/* ./

gunzip data_pairlist_july.gz

OPTIMIZE="-Kfast,openmp,ocl,simd=2,swp_strong -falign-loops -Icommon/include "
FOPTIMIZE="${OPTIMIZE} -Knoalias -Cpp -fw "
COPTIMIZE="${OPTIMIZE} -DDISABLE_VALIDATION "


frt -c ${FOPTIMIZE}  initialize_data.f90
frt -c ${FOPTIMIZE}  main.f90
frt -c ${FOPTIMIZE}  kernel_pairlist_july.f90
frt -c ${FOPTIMIZE}  sub_data_io.f90
#	fcc -c ${COPTIMIZE} common/src/report.c

frt ${FOPTIMIZE} initialize_data.o main.o kernel_pairlist_july.o sub_data_io.o

export OMP_NUM_THREADS=12
export FLIB_TRACEBACK_MEM_SIZE=128
time ./a.out
