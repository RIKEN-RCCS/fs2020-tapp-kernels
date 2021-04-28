#!/bin/bash
#PJM -N FFB-SPMMV
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

TMPDIR=${HOME}/tmp/check_FFB-SPMMV
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *

SRC_DIR=${HOME}/fs2020_kernels/src/FFB.spmmv_vec8.20170509
cp -rp ${SRC_DIR}/* ./
gunzip data_ffb_spmmv.gz

OPTIMIZE="-Kfast,openmp,ocl -Icommon/include "
FOPTIMIZE="${OPTIMIZE} -Cpp -fw "
COPTIMIZE="${OPTIMIZE} -DDISABLE_VALIDATION -DCLEAN_CACHE "

frt -c ${FOPTIMIZE}  main.f90
frt -c ${FOPTIMIZE}  ax4.f
frt -c ${FOPTIMIZE}  initialize_data.f90
frt -c ${FOPTIMIZE}  sub_data_io.f90
fcc -c ${COPTIMIZE} common/src/report.c

frt ${FOPTIMIZE} main.o ax4.o initialize_data.o	sub_data_io.o report.o

export OMP_NUM_THREADS=12
time ./a.out

ls -l
exit
