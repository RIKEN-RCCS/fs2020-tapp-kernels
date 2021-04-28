#include "qws.h"
#define _PFI1 1
#define _PFI2 1
#define _PFI3 2
#define _PFI4 1
#define _PFI5 1
#include "qwsintrin.h"
#ifndef EML_LIB
#include <complex>
#else
#include "eml_lib.h"
#endif
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "time.h"
#include "tools.h"

#ifdef __cplusplus
extern "C"{
#endif

  using namespace std;
  extern int vols;
  extern int volse;
#pragma omp threadprivate(volse)
  extern void prec_ddd_s_(scs_t* out, scs_t* in, int* nsap, int* nm);

  //  extern void check_timing_ (const char *);

  void bicgstab_precdd_s_(scs_t* x, scs_t* b, double* tol, int* conviter, int* maxiter, int* nsap, int* nm){
    __attribute__((aligned(64))) static scs_t *q, *r, *p , *t, *r0;
    int ret=0;
    if( q==0) ret = posix_memalign((void **)&q,256,sizeof(scs_t)*vols*2);
    if( r==0) ret = posix_memalign((void **)&r,256,sizeof(scs_t)*vols*2);
    if( p==0) ret = posix_memalign((void **)&p,256,sizeof(scs_t)*vols*2);
    if( t==0) ret = posix_memalign((void **)&t,256,sizeof(scs_t)*vols*2);
    if( r0==0) ret = posix_memalign((void **)&r0,256,sizeof(scs_t)*vols*2);
    float bnorm, rnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< float > rho0, rho, beta, omega, alpha, ctmp;
    //rvecs_t rvd0, rvd1, rvd2, rvd3, rvd4, rvd5;
    //rvecs_t xr, xi, rr, ri, pr, pi;
    float ar, ai, br, bi ,cr, ci;
    int i, j, iter;

    _BCG_PRECDDS_TIC_;

#if 0
    //------------------------------------------------------------------ flexible start
    //    _PREC_DDD_S_TIC_;
    prec_ddd_s_(q, x, nsap, nm);
    //    _PREC_DDD_S_TOC_;
    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          r[i].ccs[j].v[v]  = b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
          p[i].ccs[j].v[v]  = r[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = r[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v]*b[i].ccs[j].v[v];
          rtmp1 += r[i].ccs[j].v[v]*r[i].ccs[j].v[v];
        }
      }
    }
    redu[0] = rtmp0;
    redu[1] = rtmp1;
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<float>(redu[1],0);
#else
    //------------------------------------------------------------------ x=0 start
    rtmp0 = (float)0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          x[i].ccs[j].v[v]  = 0;
          r[i].ccs[j].v[v]  = b[i].ccs[j].v[v];
          p[i].ccs[j].v[v]  = b[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = b[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v]*b[i].ccs[j].v[v];
        }
      }
    }
    redu[0] = rtmp0;
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<float>(redu[0],0);
#endif

    for (iter=0; iter<(*maxiter);iter++){
      _BCG_PRECDDS_ITER_TIC_;
      // q = Ap
      //    _PREC_DDD_S_TIC_;
      prec_ddd_s_(q, p, nsap, nm);
      //    _PREC_DDD_S_TOC_;

      // alpha = rho0 / <r0,q>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vols*2; i++){
        if (i+_PFI1 < volse*2) {
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((q+i+_PFI1 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((q+i+_PFI1 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((q+i+_PFI1 )->c_prefetch)[128]),0,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * q[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * q[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * q[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * q[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      _BCG_PRECDDS_ITER_REDUC2_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC2_TOC_;
#endif
      ctmp = complex<float>(redu[0], redu[1]);
      alpha = rho0 / ctmp;

      // x = x + alpha p
      // r = r - alpha q
      ar=alpha.real();
      ai=alpha.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
        if (i+_PFI2 < volse*2) {
          __builtin_prefetch(&(((p+i+_PFI2 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((p+i+_PFI2 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((p+i+_PFI2 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((q+i+_PFI2 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((q+i+_PFI2 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((q+i+_PFI2 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[0]),1,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[64]),1,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[128]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[0]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[64]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[128]),1,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  p[i].cs[j][0].v[v] - ai *  p[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  p[i].cs[j][1].v[v] + ai *  p[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
      for(i=0; i<vols*2; i++){
        if (i+_PFI3 < volse*2) {
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[128]),0,1);
        }
        for(j=0; j<24; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      _BCG_PRECDDS_ITER_REDUC1_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC1_TOC_;
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
      if (sqrt(rnorm/bnorm) < *tol){
	_BCG_PRECDDS_ITER_TOC_;
	break;
      }

      // t = Ar
      //      _PREC_DDD_S_TIC_;
      prec_ddd_s_(t, r, nsap, nm);
      //      _PREC_DDD_S_TOC_;

      // omega = <t,r> / |t|^2
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1, rtmp2)
      for(i=0; i<vols*2; i++){
        if (i+_PFI4 < volse*2) {
          __builtin_prefetch(&(((r+i+_PFI4 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r+i+_PFI4 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r+i+_PFI4 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((t+i+_PFI4 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((t+i+_PFI4 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((t+i+_PFI4 )->c_prefetch)[128]),0,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += (t[i].cs[j][0].v[v] * t[i].cs[j][0].v[v] + t[i].cs[j][1].v[v] *  t[i].cs[j][1].v[v]);
            rtmp1 += (t[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + t[i].cs[j][1].v[v] *  r[i].cs[j][1].v[v]);
            rtmp2 += (t[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - t[i].cs[j][1].v[v] *  r[i].cs[j][0].v[v]);
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      redu[2] = rtmp2;
#ifdef _MPI_
      _BCG_PRECDDS_ITER_REDUC3_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC3_TOC_;
#endif
      omega = complex<float>(redu[1]/redu[0], redu[2]/redu[0] );

      // x = x + omega r
      // r = r - omega t
      ar=omega.real();
      ai=omega.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
        if (i+_PFI2 < volse*2) {
          __builtin_prefetch(&(((t+i+_PFI2 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((t+i+_PFI2 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((t+i+_PFI2 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[0]),1,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[64]),1,1);
          __builtin_prefetch(&(((x+i+_PFI2 )->c_prefetch)[128]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[0]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[64]),1,1);
          __builtin_prefetch(&(((r+i+_PFI2 )->c_prefetch)[128]),1,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  r[i].cs[j][0].v[v] - ai *  r[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  r[i].cs[j][1].v[v] + ai *  r[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  t[i].cs[j][0].v[v] + ai *  t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  t[i].cs[j][1].v[v] - ai *  t[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
      for(i=0; i<vols*2; i++){
        if (i+_PFI3 < volse*2) {
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r+i+_PFI3 )->c_prefetch)[128]),0,1);
        }
        for(j=0; j<24; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      _BCG_PRECDDS_ITER_REDUC1_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC1_TOC_;
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
      if (sqrt(rnorm/bnorm) < *tol){
	_BCG_PRECDDS_ITER_TOC_;
	break;
      }


      // rho = <r0,r>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vols*2; i++){
        if (i+_PFI1 < volse*2) {
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r0+i+_PFI1 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((r+i+_PFI1 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r+i+_PFI1 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r+i+_PFI1 )->c_prefetch)[128]),0,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      _BCG_PRECDDS_ITER_REDUC2_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC2_TOC_;
#endif
      rho = complex<float>(redu[0],redu[1]);

      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;


      // p = p - omega q
      // p = r + beta p
      ar = omega.real();
      ai = omega.imag();
      br = beta.real();
      bi = beta.imag();
#pragma omp parallel for private(i, j, cr, ci)
      for(i=0; i<vols*2; i++){
        if (i+_PFI5 < volse*2) {
          __builtin_prefetch(&(((q+i+_PFI5 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((q+i+_PFI5 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((q+i+_PFI5 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((r+i+_PFI5 )->c_prefetch)[0]),0,1);
          __builtin_prefetch(&(((r+i+_PFI5 )->c_prefetch)[64]),0,1);
          __builtin_prefetch(&(((r+i+_PFI5 )->c_prefetch)[128]),0,1);
          __builtin_prefetch(&(((p+i+_PFI5 )->c_prefetch)[0]),1,1);
          __builtin_prefetch(&(((p+i+_PFI5 )->c_prefetch)[64]),1,1);
          __builtin_prefetch(&(((p+i+_PFI5 )->c_prefetch)[128]),1,1);
        }
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            cr = p[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            ci = p[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
            p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
            p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
          }
        }
      }
      _BCG_PRECDDS_ITER_TOC_;
    }//iter
    *conviter = iter;
    _BCG_PRECDDS_TOC_;
  }//bicgstab


#ifdef __cplusplus
}
#endif
