#!/bin/bash
#PJM -N GAMERA-KERNEL
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "node=1"
#PJM -j
#PJM -S

module list
date
hostname

TMPDIR=${HOME}/tmp/check_GAMERA
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *.o

SRC_DIR=${HOME}/fs2020_kernels/src/GAMERA.TIMER_COMP_MATVEC_IF
cp -rp ${SRC_DIR}/* ./
gunzip data_file_77.gz

OPTIMIZE="-Kfast,openmp,ocl -Kloop_nofission -Icommon/include -DREAL_4=real -DMAXCOLOR=100 "
FOPTIMIZE="${OPTIMIZE} -Cpp "
COPTIMIZE="${OPTIMIZE} -DDISABLE_VALIDATION "

set -x
frt -c ${FOPTIMIZE}  main.F
frt -c ${FOPTIMIZE}  -fs -o cal.o  cal_amat_tet10_s_inline_simd_color_loopdiv_r3_proposed1.F
frt -c ${FOPTIMIZE}  sub_data_io.F

frt ${FOPTIMIZE} main.o cal.o sub_data_io.o
export OMP_NUM_THREADS=12
time ./a.out

exit
