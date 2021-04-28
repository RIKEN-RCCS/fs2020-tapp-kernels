#!/bin/bash
#PJM -N NICAM-HADV-FLUX
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

TMPDIR=${HOME}/tmp/check_NICAM.hadv_flux
mkdir -p ${TMPDIR}
cd ${TMPDIR}/
if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm *

SRC_DIR=${HOME}/fs2020_kernels/src/NICAM.Horizontal_Adv_flux.20170418
cp -rp ${SRC_DIR}/* ./

OPTIMIZE="-Kfast,openmp,ocl,loop_nofission -falign-loops -Icommon/include  -K__control=0x4 "
FOPTIMIZE="${OPTIMIZE} -Cpp -fw -DUSE_FAPP -DSINGLE "


frt -c ${FOPTIMIZE} postK_nicam.Horizon_adv_fluxmod_src_tracer.f90
frt -c ${FOPTIMIZE} main.f90
frt ${FOPTIMIZE} *.o

ls -lt

export OMP_NUM_THREADS=12
export FLIB_TRACEBACK_MEM_SIZE=128
time ./a.out

