This file is best viewed in markdown viewer

# the kernel codes from Priority Issue Target Applications

                                                          First  draft 2021/3/21 K.Mikami
                                                          Second draft 2021/4/27 K.Mikami


# [1] About this software

This software is the collection of the kernel codes from Priority Issue Target Applications,
i.e.  TAPP-kernels.
A brief description of the targer applications can be found in the Fugaku Web page at
https://www.r-ccs.riken.jp/jp/fugaku
as well as the detail descriptions linked from the same page.
The TAPP-kernels represent the typical computational work load of a subset of the target applications.

The kernels are setup to complete the timer sections in less than 0.01 seconds on Fugaku's 12-core CMG
per iteration.
Such timing granularity was chosen in order to run the kernels on simulator platforms
within reasonable wall time. The simulator platforms typically require orders of magnitude
longer wall time compared to real platforms.
The kernel code users may change the source code to adjust the wall time when measuring on the real platforms
by simply increasing the number of iterations of the outermost loop structure.
The timer sections are sandwitched between PROF_START and PROF_STOP.

The kernels are expected to run in single process multi thread mode such as 12 threads on Fugaku CMG.
Approproate compile options should be provided for each of the kernels.
For reference purpose, the compiler options tested during the early codesign stage are included
in the file "tested_early_options.sh".
Note that the options shown in "tested_early_options.sh" are also taken from the early design stage
compilers around 2017 and certain options are obsolete by today.
If the kernel must be run in single thread, certain changes must be made in the source code.
Such changes are listed in section [4] below.

The authors of the original target applications have generously agreed to share them as the
open software based on the license condition described bottom in the hope to provide the example
material for the systems development.

Remark that TAPP-kernels were taken out of the early codesign-in-progress applications for evaluation purpose,
and the source codes do not match the final optimzed version that achieved the official performance
numbers published in the above Web page.


# [2] Material

The kernel codes are stored in the directory structure shown below.

<pre>
tapp-kernels/				: source code top directory which includes below:
Adventure.region0.20170821/		: Adventure DomainFEM (region0)    kernel directory
Adventure.region1.20160808/		: Adventure CoarseMatVec (region1) kernel directory
Adventure.region2.20170727/		: Adventure CoarseMatVec (region2) kernel directory
FFB.callap_kernel2.20160805/		: FFB callap_kernel2 kernel directory
FFB.spmmv_vec8.20170509/		: FFB spmmv_vec8 kernel directory
GAMERA.TIMER_COMP_MATVEC_IF/		: GAMERA Element-by-Element kernel directory
GENESIS.PairList_July/			: GENESIS PairList Dec-July kernel directory
GENESIS.PairList_June/			: GENESIS PairList June kernel directory
GENESIS.kernel_July/			: GENESIS Nonb15F Dec-July	 kernel directory
GENESIS.kernel_June/			: GENESIS Nonb15F June1 kernel directory
NICAM.Horizontal_Adv_flux.20170418/  	: NICAM Horizontal_Adv_flux kernel directory
NICAM.Horizontal_Adv_limiter.20170407/	: NICAM Horizontal_Adv_limiter  kernel directory
NICAM.Vertical_Adv_limiter.20160902/	: NICAM Vertical_Adv_limiter kernel directory
NICAM.diffusion.20170220/		: NICAM diffusion	 kernel directory
NICAM.divdamp.20170420/			: NICAM divdamp kernel directory
NICAM.vi_rhow_solver.20160902/		: NICAM vi_rhow_solver kernel directory
QCD.ddd_in_s_.20170621/			: LQCD ddd_in_s kernel directory
QCD.jinv_ddd_in_s_.20170720/		: LQCD jinv_ddd_in_s kernel directory
streamlike_pattern1/			: NICAM stream like computing kernel 1 directory
streamlike_pattern2/			: NICAM stream like computing kernel 2 directory
streamlike_pattern3/			: NICAM stream like computing kernel 3 directory

Readme.md				: This Readme file
LICENSE.GPLv3				: license information.GPLv3.  See [5] of this Readme.md
LICENSE.LGPLv3				: license information.LGPLv3.  See [5] of this Readme.md
</pre>
<!-- -->

# [3] Instructions to setup multi thread source.

The kernel source files listed in previous [2] are all setup for OpenMP multi thread run by default.
There should not be any changes required for compilation with compilers supporting OpenMP.
For testers' convenience, each directory contains the compilation and execution shell script
for the supercomputer Fugaku as examples.
Thees shell script files are named as "x.compile-fugaku-${tapp-kernel-name}.sh"
While these example have been tested on Fugaku, some compile options are imported
from the compilers at early development stages and do not appear in the official user manual.
Note that such options are not mandatory, and are included only as the supplementary information.


# [4] Instructions to setup single thread source.

In order to run the kernels in on single thread mode,
some of the kernels need source code changes as explained below.

## 3. GENESIS Nonb15F June1
###   1) line 60 of initialize_data.f90
        gparam%nthread = omp_get_num_threads() -> gparam%nthread=1
###   2) line 75 of kernel_June_1.f90
        id = omp_get_thread_num() -> id=0

## 4. GENESIS Nonb15F Dec-July
###   1) line 70 of initialize_data.f90
        gparam%nthread = omp_get_num_threads() -> gparam%nthread=1
###   2) line 85 of kernel_july.f90
        id = omp_get_thread_num() -> id=0

## 5. GENESIS PairList June
###   1) line 56 of initialize_data.f90
        gparam%nthread = omp_get_num_threads() -> gparam%nthread=1
###   2) line 64 of kernel_pairlist_June_1.f90
        id = omp_get_thread_num() -> id=0

## 6. GENESIS PairList Dec-July
###   1) line 59 of initialize_data.f90
        gparam%nthread = omp_get_num_threads() -> gparam%nthread=1
###   2) line 64 of kernel_pairlist_july.f90
        id = omp_get_thread_num() -> id=0

## 12. NICAM Vertical_Adv_limiter
###   1) line 261, set num_threads as 1
        num_threads = omp_get_num_threads()→  num_threads = 1

## 13. NICAM vi_rhow_solve
###   1) line 309, set num_threads as 1
        num_threads = omp_get_num_threads()→  num_threads = 1

## 18. Adventure CoarseMatVec(region1)
###   1) line 65, set thread_num as 0
        thread_num = omp_get_thread_num() → thread_num = 0

## 22. LQCD ddd_in_s
###   1) Add -D_OPENMP in compiler option
###   2) in qws.cc line 65, set iam as 0
        int iam = omp_get_thread_num(); → int iam = 0;
###   3) in qws.cc line 68, set thmax as 1
        thmax = omp_get_num_threads(); →  thmax = 1;

## 23. LQCD jinv_ddd_in_s
###   1) Add -D_OPENMP in compiler option
###   2) in qws.cc line 102, set iam as 0
        iam = omp_get_thread_num(); → iam = 0;
###   3) in qws.cc line 106, set thamx as 1
        thmax = omp_get_num_threads(); →  thmax = 1;

## end of Instructions 


# [5] Copyright notice and LICENSE information

This collection of the kernel codes from Priority Issue Target Applications,
named TAPP-kernels hereafter,
is distributed under the GNU Lesser General Public License version 3.

  Copyright (C) 2021 Flagship 2020 Project, RIKEN Center for Computational Science

  TAPP-kernels is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  TAPP-kernels is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TAPP-kernels -- see the files
  LICENSE.LGPLv3 and LICENSE.GPLv3.
  If not, see https://www.gnu.org/licenses/.

