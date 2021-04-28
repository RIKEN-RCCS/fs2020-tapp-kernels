/*
!===============================================================================
!
! Copyright (C) 2014,2015 Yoshifumi Nakamura
!
! _CODENAME_ = ???
!
! This file is part of _CODENAME_
!
! _CODENAME_ is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! _CODENAME_ is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with _CODENAME_.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------------
*/
#include "qws.h"
#include "qwsintrin.h"
#include "report.h"
#include "profiler.h"

#include "init.h"
#include "tools.h"

extern int vold, vols, rank, nx, ny, nz, nt, nxh, nxd, nxs;
extern double kappa2, kappa, mkappa;
extern int prof_flag;
float* result_arr;
int result_arr_len;
unsigned char *malloc_p;
size_t         malloc_len = 0;
extern "C" void qws_init_(int* , int* , int* , int* , int* , int* , int* , int* , int* );
extern "C" void bicgstab_precdd_s_(scs_t* , scs_t* , double* , int* , int* , int* , int*);

int main()
{
  PROF_INIT;
  PROF_START_ALL;
  
  int lx=NX;
  int ly=NY;
  int lz=NZ;
  int lt=NT;
  int npe_f[4];
  npe_f[0]=1;
  npe_f[1]=1;
  npe_f[2]=1;
  npe_f[3]=1;

  double tol_s= -1.0;

  int fbc_f[4];
  fbc_f[0]=1;
  fbc_f[1]=1;
  fbc_f[2]=1;
  fbc_f[3]=-1;
  int pce_f=0;
  int pco_f=1;
  
  int block_size[4];
  block_size[0]=1;
  block_size[1]=2;
  block_size[2]=2;
  block_size[3]=2;

  __attribute__((aligned(256))) static unsigned char malloc_area[MALLOCAREA];
  malloc_p = malloc_area;
  
  qws_init_(&lx, &ly, &lz, &lt, npe_f, fbc_f, &pce_f, &pco_f, block_size);

  kappa = (double)0.05;
  kappa2 = kappa * kappa;
  mkappa = - kappa;

  int iter = 1;
  int maxiter = 1;
  int nsap = 4;
  int nm = 2;

  __attribute__((aligned(64))) static scs_t b_s[VOLS*2];
  __attribute__((aligned(64))) static scs_t t_s[VOLS*2];

  fill_rand((float*)b_s, (float*)(b_s)+sizeof(b_s)/sizeof(float));
  fill_rand((float*)t_s, (float*)(t_s)+sizeof(t_s)/sizeof(float));

  prof_flag = 0;
  bicgstab_precdd_s_(t_s, b_s, &tol_s, &iter, &maxiter, &nsap, &nm);

  prof_flag = 1;
  bicgstab_precdd_s_(t_s, b_s, &tol_s, &iter, &maxiter, &nsap, &nm);

  result_arr = (float*)t_s;
  result_arr_len = sizeof(t_s)/sizeof(float);
//result_arr_len = 55296;

#ifdef DS_TO_DOUBLE
  report_validation(get_ss_r8(result_arr, result_arr_len), -4.180054191179794e+307, 0.01);
#else
  report_validation(get_ss_r4(result_arr, result_arr_len),  4.593723632812500e+03, 0.01);
#endif

  PROF_STOP_ALL;
  PROF_FINALIZE;

  return 0;
}
