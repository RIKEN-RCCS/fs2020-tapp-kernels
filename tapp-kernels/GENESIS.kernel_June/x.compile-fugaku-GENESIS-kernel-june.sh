#!/bin/bash
#PJM -N GENE-KERNEL-JUNE
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

TMPDIR=${HOME}/tmp/check_GENESIS-kernel-june
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *


SRC_DIR=${HOME}/fs2020_kernels/src/GENESIS.kernel_June

cp -rp ${SRC_DIR}/* ./
gunzip data_kernel_June.gz

OPTIMIZE="-Kfast,openmp,ocl,simd=2,nounroll,swp_strong -falign-loops -Icommon/include "
FOPTIMIZE="${OPTIMIZE} -Knoalias -Cpp -fw "
#	FOPTIMIZE="${OPTIMIZE} -Knoalias -Cpp -fw -Nquickdbg=subchk -g "
COPTIMIZE="${OPTIMIZE} -DDISABLE_VALIDATION "


frt -c ${FOPTIMIZE}  initialize_common.f90
frt -c ${FOPTIMIZE}  main.f90
frt -c ${FOPTIMIZE}  kernel_June_1.F90
frt -c ${FOPTIMIZE}  sub_data_io.f90
fcc -c ${COPTIMIZE} common/src/report.c

frt ${FOPTIMIZE} initialize_common.o main.o kernel_June_1.o sub_data_io.o report.o

#	export OMP_NUM_THREADS=12
export OMP_NUM_THREADS=3
time ./a.out
pwd
#	ls -l
exit

export OMP_NUM_THREADS=4
time ./a.out
pwd
#	ls -l
#	exit
