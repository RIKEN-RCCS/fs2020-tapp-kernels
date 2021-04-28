#define _PFI 1
#define _PFIX 2
#define _PFIT 1
#include "qws.h"
#include "wilson_s.h"
#include "clover_s.h"
#include "mult_all.h"
#include "ddd_out_s_0_inline.h"
#include <omp.h>

#include "time.h"
#include "tools.h"

#ifdef __cplusplus
extern "C"{
#endif

  extern pglus_t __restrict__ glus __attribute__((aligned(_ALIGN_SIZE)));
  extern pclvs_t __restrict__ clvs __attribute__((aligned(_ALIGN_SIZE)));
  extern projscs1_t * __restrict__ xfs_send;
  extern projscs1_t * __restrict__ xfs_recv;
  extern __restrict__ projscs1_t *xbs_send, *xbs_recv;
  extern __restrict__ projscs_t  *yfs_send, *yfs_recv;
  extern __restrict__ projscs_t  *ybs_send, *ybs_recv;
  extern __restrict__ projscs_t  *zfs_send, *zfs_recv;
  extern __restrict__ projscs_t  *zbs_send, *zbs_recv;
  extern __restrict__ projscs_t  *tfs_send, *tfs_recv;
  extern __restrict__ projscs_t  *tbs_send, *tbs_recv;

  extern double kappa, mkappa;
  extern int nxs, ny, nz, nt;
  extern int vols;
  extern int thmax;
  extern int iam;
#pragma omp threadprivate(iam)
  extern int pt;
  extern int npe[4];
  extern double fbc[4][2];

  void xbound(int req, int prec);
  void xbound_wait(int req, int prec);
  void xbound_send_waitall(int prec);
  void xbound_recv_waitall(int prec);

  //---------------------------------------------------------------------------------------- preprocess mult D for boundary
    //
    // pack data for send
    //
  void ddd_out_pre_s_noprl_(scs_t* __restrict__  in, int* __restrict__  idomain) {

    _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_;

    int xdmn, ydmn, zdmn, tdmn;

    xdmn = 1 - *idomain;
    if (npe[1] == 1) { ydmn = *idomain; } else { ydmn = 1 - *idomain; }
    if (npe[2] == 1) { zdmn = *idomain; } else { zdmn = 1 - *idomain; }
    if (npe[3] == 1) { tdmn = *idomain; } else { tdmn = 1 - *idomain; }

    const pglus_t __restrict__ gx = &glus[vols*0 + NDIM*vols*xdmn];
    const pglus_t __restrict__ gy = &glus[vols*1 + NDIM*vols*ydmn];
    const pglus_t __restrict__ gz = &glus[vols*2 + NDIM*vols*zdmn];
    const pglus_t __restrict__ gt = &glus[vols*3 + NDIM*vols*tdmn];


    int num_t = thmax/6;
    int num_t2;
    int nmodnum = thmax%6;
    int nmodx,nmody,nmodz,nmodt;
    num_t2 =num_t*2;
    nmodx= (ny * nz * nt)%num_t2;
    nmody= (nxs * nz * nt)%num_t;
    nmodz= nt%num_t;
    nmodt= (nxs * ny * nz)%(thmax-num_t2-num_t*2);
    if (nmodx == 0 && nmody == 0 && nmodz == 0 && nmodt == 0 && nmodnum == 0 && nt < thmax){
      if (iam < num_t2) {
        int xs = nt*nz*ny/(num_t2) * iam;
        int xe = nt*nz*ny/(num_t2) * (iam + 1);
   

        for (int i0 = xs; i0 < xe; ++i0) {
          if ( i0+_PFIX != xe) {
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(xfs_send[i0+_PFIX].c[0][0][0]),1,1);
            __builtin_prefetch(&(xbs_send[i0+_PFIX].c[0][0][0]),1,1);
#endif
          }
      // X-forward
          {
            int ixf = i0*nxs;
            const scs_t* __restrict__ in_ixf = in+ixf+vols*xdmn;
            __mult_x_forw_pre_3_(xfs_send[i0],(*(in_ixf)));
          }
      // X-backward
          { 
            int ixb = nxs-1 + i0*nxs;
            const __restrict__ pglus_t gx_ixb = gx+ixb;
            const scs_t* __restrict__ in_ixb = in+ixb+vols*xdmn;
            projscs1_t a __attribute__((aligned(_ALIGN_SIZE)));
            __mult_x_back_pre_3_(a,(*(in_ixb)));
            __mult_udag_y_3_(xbs_send[i0],a,(*(gx_ixb)));
          }
        }
      } //if num_t

      if (iam > num_t2-1 && iam < num_t2+num_t) {
        int ys = nt*nz/num_t * (iam - num_t2);
        int ye = nt*nz/num_t * (iam - num_t2 + 1);
        int i2=ys*nxs;
        for (int i0 = ys; i0 < ye; ++i0) {
          for (int x = 0; x < nxs; ++x,++i2) {
            //
            // Y-forward send
            //
            int iyf = x + nxs*ny*i0;
            int iyb = x +  nxs*ny*i0 + (ny-1)*nxs;
            if ( x !=nxs-1) {
#if VLENS == 16
#elif VLENS == 8
              __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(yfs_send[i2+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(yfs_send[i2+_PFI].cv[2][0][0]),1,1);
              __builtin_prefetch(&(ybs_send[i2+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(ybs_send[i2+_PFI].cv[2][0][0]),1,1);
#endif
            } else if(i0 != ye-1) {
#if VLENS == 16
#elif VLENS == 8
              __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(yfs_send[i2+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(yfs_send[i2+_PFI].cv[2][0][0]),1,1);
              __builtin_prefetch(&(ybs_send[i2+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(ybs_send[i2+_PFI].cv[2][0][0]),1,1);
#endif
            }
            const scs_t* __restrict__ in_iyf = in+iyf+vols*ydmn;
            __mult_y_forw_pre_(yfs_send[i2],(*(in_iyf)));
            //
            // Y-backward send
            //
            const __restrict__ pglus_t gy_iyb = gy+iyb;
            const scs_t* __restrict__ in_iyb = in+iyb+vols*ydmn;
            __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;
            __mult_y_back_pre_(a,(*(in_iyb)));
            __mult_udag_y_(ybs_send[i2],a,(*(gy_iyb)));
          }
        }
      } //if num_t

      if (iam > num_t2+num_t-1 && iam < num_t2+num_t*2) {
        int nsize=nt/num_t ;
        int zs = nsize * (iam-num_t2-num_t);
        int ze = zs + nsize ;
        for (int t = zs; t < ze; ++t) {
          int i4=nxs*ny*t;
          for (int i0 = 0; i0 < nxs*ny; ++i0,++i4) {
            //
            // Z-forward send
            //
            int izf = i0 + t*nxs*ny*nz;
            int izb = i0 + t*nxs*ny*nz + (nz-1)*nxs*ny;
            if ( i0 !=nxs*ny-1) {
#if VLENS == 16
#elif VLENS == 8
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(zfs_send[i4+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(zfs_send[i4+_PFI].cv[2][0][0]),1,1);
              __builtin_prefetch(&(zbs_send[i4+_PFI].cv[0][0][0]),1,1);
              __builtin_prefetch(&(zbs_send[i4+_PFI].cv[2][0][0]),1,1);
#endif
            }
            const scs_t* __restrict__ in_izf = in+izf+vols*zdmn;
            __mult_z_forw_pre_(zfs_send[i4],(*(in_izf)));
            //
            // Z-backward send
            //

            const scs_t* __restrict__ in_izb = in+izb+vols*zdmn;
            const __restrict__ pglus_t gz_izb = gz+izb;
            __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

            __mult_z_back_pre_(a,(*(in_izb)));
            __mult_udag_y_(zbs_send[i4],a,(*(gz_izb)));
          }
        }
      } //if num_t

      if (iam > num_t2+num_t*2-1 ) {
        int t_threads= thmax-num_t*2-num_t2;
        int ts = nz*ny*nxs/(t_threads) * (iam- thmax+t_threads);
        int te = nz*ny*nxs/(t_threads) * (iam- thmax+t_threads+1);
        for (int i0 = ts; i0 < te; ++i0) {
          int itb = i0 + (nt-1)*nxs*ny*nz;
            if ( i0 !=te-1) {
#if VLENS == 16
#elif VLENS == 8
              __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[143]),0,1);
              __builtin_prefetch(&(tfs_send[i0+_PFIT].cv[0][0][0]),1,1);
              __builtin_prefetch(&(tfs_send[i0+_PFIT].cv[2][0][0]),1,1);
              __builtin_prefetch(&(tfs_send[i0+_PFIT].c[2][1][1][VLENS-1]),1,1);
              __builtin_prefetch(&(tbs_send[i0+_PFIT].cv[0][0][0]),1,1);
              __builtin_prefetch(&(tbs_send[i0+_PFIT].cv[2][0][0]),1,1);
              __builtin_prefetch(&(tbs_send[i0+_PFIT].c[2][1][1][VLENS-1]),1,1);
#endif
            } else{
#if VLENS == 16
#elif VLENS == 8
#endif
            }
            //
            // T-forward send
            //
            const scs_t* __restrict__ in_itf = in+i0+vols*tdmn;
            __mult_t_forw_pre_(tfs_send[i0],(*(in_itf)));
            //
            // T-backward send
            //
            const scs_t* __restrict__ in_itb = in+itb+vols*tdmn;
            const __restrict__ pglus_t gt_itb = gt+itb;
            __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;
            __mult_t_back_pre_(a,(*(in_itb)));
            __mult_udag_y_(tbs_send[i0],a,(*(gt_itb)));
        }
      } //if num_t
    } else{
#pragma omp for collapse(3) nowait
      for (int t = 0; t < nt; ++t) {
        for (int z = 0; z < nz; ++z) {
          for (int y = 0; y < ny; ++y) {

            int it = nxs*ny*nz*t;
            int iz = nxs*ny*z;
            int iy = nxs*y;

            // X-forward
            {
              int ixf = 0 + iy + iz + it;
              int i0  = y + ny*z + ny*nz*t;

              __mult_x_forw_pre_3_(xfs_send[i0],in[vols*xdmn+ixf]);
            }

            // X-backward
            { 
              int ixb = nxs-1 + iy + iz + it;
              int i1  = y + ny*z + ny*nz*t;
              projscs1_t a __attribute__((aligned(_ALIGN_SIZE)));

              __mult_x_back_pre_3_(a,in[vols*xdmn+ixb]);
              __mult_udag_y_3_(xbs_send[i1],a,(*(gx + ixb)));
            }
          }
        }
      }

#pragma omp for collapse(3) nowait
      for (int t = 0; t < nt; ++t) {
        for (int z = 0; z < nz; ++z) {
          for (int x = 0; x < nxs; ++x) {

            int it = nxs*ny*nz*t;
            int iz = nxs*ny*z;

            //
            // Y-forward send
            //
            {
              int iyf = x + iz + it;
              int i2  = x + nxs*z + nxs*nz*t;

              __mult_y_forw_pre_(yfs_send[i2],in[vols*ydmn+iyf]);
            }

            //
            // Y-backward send
            //
            {
              int iyb = x + iz + it + (ny-1)*nxs;
              int i3  = x + nxs*z + nxs*nz*t;

              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

              __mult_y_back_pre_(a,in[vols*ydmn+iyb]);
              __mult_udag_y_(ybs_send[i3],a,(*(gy + iyb)));
            }
          }
        }
      }

#pragma omp for collapse(3) nowait
      for (int t = 0; t < nt; ++t) {
        for (int y = 0; y < ny; ++y) {
          for (int x = 0; x < nxs; ++x) {

            int it = nxs*ny*nz*t;
            int iy = nxs*y;

            //
            // Z-forward send
            //
            { 
              int izf = x + iy + it;
              int i4  = x + nxs*y + nxs*ny*t;

              __mult_z_forw_pre_(zfs_send[i4],in[vols*zdmn+izf]);
            }

            //
            // Z-backward send
            //
            {
              int izb = x + iy + it + (nz-1)*nxs*ny;
              int i5  = x + nxs*y + nxs*ny*t;

              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

              __mult_z_back_pre_(a,in[vols*zdmn+izb]);
              __mult_udag_y_(zbs_send[i5],a,(*(gz + izb)));
            }
          }
        }
      }

#pragma omp for collapse(3) nowait
      for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
          for (int x = 0; x < nxs; ++x) {

            int iz = nxs*ny*z;
            int iy = nxs*y;

            //
            // T-forward send
            //
            {
              int itf = x + iy + iz;
              int i6  = x + nxs*y + nxs*ny*z;

              __mult_t_forw_pre_(tfs_send[i6],in[vols*tdmn+itf]);
            }

            //
            // T-backward send
            //
            {
              int itb = x + iy + iz + (nt-1)*nxs*ny*nz;
              int i7  = x + nxs*y + nxs*ny*z;

              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

              __mult_t_back_pre_(a,in[vols*tdmn+itb]);
              __mult_udag_y_(tbs_send[i7],a,(*(gt + itb)));
            }
          }
        }
      }
    }

#pragma omp barrier

    _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_;
    _COMM_TIC_;
#ifdef RDC
#pragma omp master
#else
#pragma omp single nowait
#endif
    {
      if (*idomain == 0) {
        memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
      } else {
        xbound(0,4);
      }
      if (*idomain == 1) {
        memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
      } else {
        xbound(1,4);
      }
      xbound(2,4);
      xbound(3,4);
      xbound(4,4);
      xbound(5,4);
      xbound(6,4);
      xbound(7,4);
    }
    _COMM_TOC_;
#pragma omp barrier
  }
  void ddd_out_pre_s_(scs_t* __restrict__  in, int* __restrict__  idomain) {
#pragma omp parallel
{
    ddd_out_pre_s_noprl_(in, idomain);
}
}


    //
    // pack data for send
    //
#ifndef RDC
  void ddd_out_pre_s_noprl_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain) {


    int xdmn, ydmn, zdmn, tdmn;

    xdmn = 1 - *idomain;
    if (npe[1] == 1) { ydmn = *idomain; } else { ydmn = 1 - *idomain; }
    if (npe[2] == 1) { zdmn = *idomain; } else { zdmn = 1 - *idomain; }
    if (npe[3] == 1) { tdmn = *idomain; } else { tdmn = 1 - *idomain; }

    const pglus_t __restrict__ gx = &glus[vols*0 + NDIM*vols*xdmn];
    const pglus_t __restrict__ gy = &glus[vols*1 + NDIM*vols*ydmn];
    const pglus_t __restrict__ gz = &glus[vols*2 + NDIM*vols*zdmn];
    const pglus_t __restrict__ gt = &glus[vols*3 + NDIM*vols*tdmn];


    int num_t = thmax/6;
    int num_t2;
    int nmodnum = thmax%6;
    int nmodx,nmody,nmodz,nmodt;
//    if(nt ==2){
//       num_t2 =num_t*2;
//    }else{
//       num_t2 =num_t*3;
//    }
       num_t2 =num_t*2;
    nmodx= (ny * nz * nt)%num_t2;
    nmody= (nxs * nz * nt)%num_t;
    nmodz= nt%num_t;
    nmodt= (nxs * ny * nz)%(thmax-num_t2-num_t*2);
    if (nmodx == 0 && nmody == 0 && nmodz == 0 && nmodt == 0 && nmodnum == 0 && nt < thmax){
//  if (iam==0)printf("fj dbg2\n");

    if (iam < num_t2) {
    int xs = nt*nz*ny/(num_t2) * iam;
    int xe = nt*nz*ny/(num_t2) * (iam + 1);
   

//#pragma omp for collapse(3) nowait
//    for (int t = 0; t < nt; ++t) {
//    for (int z = 0; z < nz; ++z) {
//    for (int y = 0; y < ny; ++y) {
    for (int i0 = xs; i0 < xe; ++i0) {
//      int it = nxs*ny*nz*t;
//      int iz = nxs*ny*z;
//      int iy = nxs*y;
       if ( i0+_PFIX != xe) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*xdmn+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gx+(i0+_PFIX)*nxs+nxs-1)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(xfs_send[i0+_PFIX].c[0][0][0]),1,1);
      __builtin_prefetch(&(xbs_send[i0+_PFIX].c[0][0][0]),1,1);
#endif
      }
      // X-forward
      {
        //int x = nxs-1;
//        int ixf = 0 + iy + iz + it;
        int ixf = i0*nxs;
//        int i0  = y + ny*z + ny*nz*t;
         const scs_t* __restrict__ in_ixf = in+ixf+vols*xdmn;
        __mult_x_forw_pre_3_(xfs_send[i0],(*(in_ixf)));
      }
      // X-backward
      { 
        //int x = 0;
//        int ixb = nxs-1 + iy + iz + it;
        int ixb = nxs-1 + i0*nxs;
//        int i1  = y + ny*z + ny*nz*t;
        const __restrict__ pglus_t gx_ixb = gx+ixb;
        const scs_t* __restrict__ in_ixb = in+ixb+vols*xdmn;
        projscs1_t a __attribute__((aligned(_ALIGN_SIZE)));
        __mult_x_back_pre_3_(a,(*(in_ixb)));
//        __mult_udag_y_3_(xbs_send[i1],a,(*(gx + ixb)));
        __mult_udag_y_3_(xbs_send[i0],a,(*(gx_ixb)));
      }
    }
    } //if num_t
//    }
//    }
    if (iam > num_t2-1 && iam < num_t2+num_t) {
    int ys = nt*nz/num_t * (iam - num_t2);
    int ye = nt*nz/num_t * (iam - num_t2 + 1);
//    int ys = nt*nz/num_t * (iam - 2*num_t);
//    int ye = nt*nz/num_t * (iam - 2*num_t + 1);
    int i2=ys*nxs;
//    int i2=(ys*nt*nz)/num_t;
//#pragma omp for collapse(3) nowait
//    for (int t = 0; t < nt; ++t) {
//    for (int z = 0; z < nz; ++z) {
    for (int i0 = ys; i0 < ye; ++i0) {
      for (int x = 0; x < nxs; ++x,++i2) {
//        int it = nxs*ny*nz*t;
//        int iz = nxs*ny*z;
        //
        // Y-forward send
        //
          //int y = ny-1;
//          int iyf = x + iz + it;
          int iyf = x + nxs*ny*i0;
          int iyb = x +  nxs*ny*i0 + (ny-1)*nxs;
       if ( x !=nxs-1) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyf+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(yfs_send[i2+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(yfs_send[i2+_PFI].cv[2][0][0]),1,1);
      __builtin_prefetch(&(ybs_send[i2+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(ybs_send[i2+_PFI].cv[2][0][0]),1,1);
#endif
            } else if(i0 != ye-1) {
#if VLENS == 16
#elif VLENS == 8
//      __builtin_prefetch(&((( in+vols*ydmn+(i0+1)*nxs*ny+_PFI-1)->c_prefetch)[0]),0,1);
//      __builtin_prefetch(&((( in+vols*ydmn+(i0+1)*nxs*ny+_PFI-1)->c_prefetch)[64]),0,1);
//      __builtin_prefetch(&((( in+vols*ydmn+(i0+1)*nxs*ny+_PFI-1)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyf+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*ydmn+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gy+iyb+nxs*ny-nxs+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(yfs_send[i2+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(yfs_send[i2+_PFI].cv[2][0][0]),1,1);
      __builtin_prefetch(&(ybs_send[i2+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(ybs_send[i2+_PFI].cv[2][0][0]),1,1);
#endif
      }
         const scs_t* __restrict__ in_iyf = in+iyf+vols*ydmn;
          __mult_y_forw_pre_(yfs_send[i2],(*(in_iyf)));
        //
        // Y-backward send
        //
          //int y = 0;
//          int iyb = x + iz + it + (ny-1)*nxs;
         const __restrict__ pglus_t gy_iyb = gy+iyb;
         const scs_t* __restrict__ in_iyb = in+iyb+vols*ydmn;
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;
          __mult_y_back_pre_(a,(*(in_iyb)));
          __mult_udag_y_(ybs_send[i2],a,(*(gy_iyb)));
      }
    }
    } //if num_t
//    }
    if (iam > num_t2+num_t-1 && iam < num_t2+num_t*2) {
    int nsize=nt/num_t ;
    int zs = nsize * (iam-num_t2-num_t);
    int ze = zs + nsize ;
    for (int t = zs; t < ze; ++t) {
//    int t = nt/num_t * (iam - num_t*3);
//    int size_z =  ny*nxs/num_t;
//    int ze = zs+size_z;
    int i4=nxs*ny*t;
//#pragma omp for collapse(3) nowait
//    for (int t = 0; t < nt; ++t) {
//    for (int y = 0; y < ny; ++y) {
//    for (int x = 0; x < nxs; ++x) {
      for (int i0 = 0; i0 < nxs*ny; ++i0,++i4) {
//        int it = nxs*ny*nz*t;
//        int iy = nxs*y;
        //
        // Z-forward send
        //
       int izf = i0 + t*nxs*ny*nz;
//       int i4  = t*nxs*ny+i0;
       int izb = i0 + t*nxs*ny*nz + (nz-1)*nxs*ny;
       if ( i0 !=nxs*ny-1) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izf)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*zdmn+_PFI+izb)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gz+izb+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(zfs_send[i4+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(zfs_send[i4+_PFI].cv[2][0][0]),1,1);
      __builtin_prefetch(&(zbs_send[i4+_PFI].cv[0][0][0]),1,1);
      __builtin_prefetch(&(zbs_send[i4+_PFI].cv[2][0][0]),1,1);
#endif
      }
          //int z = nz-1;
//          int izf = x + iy + it;
//          int i4  = x + nxs*y + nxs*ny*t;
         const scs_t* __restrict__ in_izf = in+izf+vols*zdmn;
          __mult_z_forw_pre_(zfs_send[i4],(*(in_izf)));
        //
        // Z-backward send
        //
          //int z = 0;
//          int izb = x + iy + it + (nz-1)*nxs*ny;
//          int i5  = x + nxs*y + nxs*ny*t;

         const scs_t* __restrict__ in_izb = in+izb+vols*zdmn;
         const __restrict__ pglus_t gz_izb = gz+izb;
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

          __mult_z_back_pre_(a,(*(in_izb)));
          __mult_udag_y_(zbs_send[i4],a,(*(gz_izb)));
      }
      }
    } //if num_t
//    }
//    }

//#pragma omp for collapse(3) nowait
//    for (int z = 0; z < nz; ++z) {
//    for (int y = 0; y < ny; ++y) {
//      for (int x = 0; x < nxs; ++x) {
    if (iam > num_t2+num_t*2-1 ) {
    int t_threads= thmax-num_t*2-num_t2;
    int ts = nz*ny*nxs/(t_threads) * (iam- thmax+t_threads);
    int te = nz*ny*nxs/(t_threads) * (iam- thmax+t_threads+1);
//    int ts = nz*ny*nxs/(num_t*2) * (iam- num_t*4);
//    int te = nz*ny*nxs/(num_t*2) * (iam -num_t*4 + 1);
      for (int i0 = ts; i0 < te; ++i0) {
          int itb = i0 + (nt-1)*nxs*ny*nz;
//       int iz = nxs*ny*z;
//       int iy = nxs*y;
       if ( i0 !=te-1) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*tdmn+i0+_PFIT)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( in+vols*tdmn+itb+_PFIT)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&(tfs_send[i0+_PFIT].cv[0][0][0]),1,1);
      __builtin_prefetch(&(tfs_send[i0+_PFIT].cv[2][0][0]),1,1);
      __builtin_prefetch(&(tfs_send[i0+_PFIT].c[2][1][1][VLENS-1]),1,1);
      __builtin_prefetch(&(tbs_send[i0+_PFIT].cv[0][0][0]),1,1);
      __builtin_prefetch(&(tbs_send[i0+_PFIT].cv[2][0][0]),1,1);
      __builtin_prefetch(&(tbs_send[i0+_PFIT].c[2][1][1][VLENS-1]),1,1);
#endif
            } else{
#if VLENS == 16
#elif VLENS == 8
#endif
      }
        //
        // T-forward send
        //
          //int t = nt-1;
//          int itf = x + iy + iz;
//          int i6  = x + nxs*y + nxs*ny*z;
         const scs_t* __restrict__ in_itf = in+i0+vols*tdmn;
          __mult_t_forw_pre_(tfs_send[i0],(*(in_itf)));
        //
        // T-backward send
        //
//          int i7  = x + nxs*y + nxs*ny*z;
         const scs_t* __restrict__ in_itb = in+itb+vols*tdmn;
         const __restrict__ pglus_t gt_itb = gt+itb;
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;
          __mult_t_back_pre_(a,(*(in_itb)));
          __mult_udag_y_(tbs_send[i0],a,(*(gt_itb)));
      }
    } //if num_t
//    }
//    }
    } else{
#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {

      int it = nxs*ny*nz*t;
      int iz = nxs*ny*z;
      int iy = nxs*y;

      // X-forward
      {
        //int x = nxs-1;

        int ixf = 0 + iy + iz + it;
        int i0  = y + ny*z + ny*nz*t;

        __mult_x_forw_pre_3_(xfs_send[i0],in[vols*xdmn+ixf]);

      }

      // X-backward
      { 
        //int x = 0;

        int ixb = nxs-1 + iy + iz + it;
        int i1  = y + ny*z + ny*nz*t;
        projscs1_t a __attribute__((aligned(_ALIGN_SIZE)));

        __mult_x_back_pre_3_(a,in[vols*xdmn+ixb]);
        __mult_udag_y_3_(xbs_send[i1],a,(*(gx + ixb)));

      }
    }
    }
    }

#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < nxs; ++x) {

        int it = nxs*ny*nz*t;
        int iz = nxs*ny*z;

        //
        // Y-forward send
        //
        {
          //int y = ny-1;

          int iyf = x + iz + it;
          int i2  = x + nxs*z + nxs*nz*t;

          __mult_y_forw_pre_(yfs_send[i2],in[vols*ydmn+iyf]);

        }

        //
        // Y-backward send
        //
        {
          //int y = 0;

          int iyb = x + iz + it + (ny-1)*nxs;
          int i3  = x + nxs*z + nxs*nz*t;

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

          __mult_y_back_pre_(a,in[vols*ydmn+iyb]);
          __mult_udag_y_(ybs_send[i3],a,(*(gy + iyb)));

        }
      }
    }
    }

#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nxs; ++x) {

        int it = nxs*ny*nz*t;
        int iy = nxs*y;

        //
        // Z-forward send
        //
        { 
          //int z = nz-1;

          int izf = x + iy + it;
          int i4  = x + nxs*y + nxs*ny*t;

          __mult_z_forw_pre_(zfs_send[i4],in[vols*zdmn+izf]);


        }

        //
        // Z-backward send
        //
        {
          //int z = 0;

          int izb = x + iy + it + (nz-1)*nxs*ny;
          int i5  = x + nxs*y + nxs*ny*t;

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

          __mult_z_back_pre_(a,in[vols*zdmn+izb]);
          __mult_udag_y_(zbs_send[i5],a,(*(gz + izb)));


        }
      }
    }
    }

#pragma omp for collapse(3) nowait
    for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nxs; ++x) {

       int iz = nxs*ny*z;
       int iy = nxs*y;

        //
        // T-forward send
        //
        {
          //int t = nt-1;

          int itf = x + iy + iz;
          int i6  = x + nxs*y + nxs*ny*z;

          __mult_t_forw_pre_(tfs_send[i6],in[vols*tdmn+itf]);

        }

        //
        // T-backward send
        //
        {
          //int t = 0;

          int itb = x + iy + iz + (nt-1)*nxs*ny*nz;
          int i7  = x + nxs*y + nxs*ny*z;

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a;

          __mult_t_back_pre_(a,in[vols*tdmn+itb]);
          __mult_udag_y_(tbs_send[i7],a,(*(gt + itb)));

        }

      }
    }
    }
    }

#pragma omp barrier

#pragma omp single nowait
    {
      if (*idomain == 0) {
        memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
      } else {
        xbound(0,4);
      }
      if (*idomain == 1) {
        memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
      } else {
        xbound(1,4);
      }
      xbound(2,4);
      xbound(3,4);
      xbound(4,4);
      xbound(5,4);
      xbound(6,4);
      xbound(7,4);
    }
  }
#endif
//RDC

#ifndef RDC
  void ddd_out_pre_s_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain) {
#pragma omp parallel
{
    ddd_out_pre_s_noprl_no_timer_(in, idomain);
}
}
#endif
//RDC



  //---------------------------------------------------------------------------------------- postprocess mult D for boundary
  void ddd_out_pos_s_noprl_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {

    const float lfactor = factor;
    const int lnxs = nxs, lny = ny, lnz = nz, lnt = nt;
    const int lvols = vols;
#define factor lfactor
#define nxs lnxs
#define ny lny
#define nz lnz
#define nt lnt
#define vols lvols
    const pglus_t __restrict__ gx = &glus[vols*0 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gy = &glus[vols*1 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gz = &glus[vols*2 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gt = &glus[vols*3 + NDIM*vols*(*idomain)];
    static scs_t *tmp;
#ifdef RDC
#pragma omp master
#else
#pragma omp single
#endif
    {
      if( tmp==0)  posix_memalign((void **)&tmp,256,sizeof(scs_t)*vols);
    }
#ifdef RDC
#pragma omp barrier
#endif

    _COMM_TIC_;
#ifdef RDC
#pragma omp master
#else
#pragma omp single
#endif
    {
      xbound_wait(0,4);
      xbound_wait(1,4);
      xbound_wait(2,4);
      xbound_wait(3,4);
      xbound_wait(4,4);
      xbound_wait(5,4);
      xbound_wait(6,4);
      xbound_wait(7,4);
    }
#ifdef RDC
#pragma omp barrier
#endif
    _COMM_TOC_;
    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_;
    int nmodx,nmody,nmodz,nmodt;
    nmodx= (ny * nz * nt)%thmax;
    nmody= ( nz * nt)%thmax;
    nmodz= ( ny * nt)%thmax;
    nmodt= ( ny * nz)%thmax;
    if (nmodx == 0 && nmody == 0 && nmodz == 0 && nmodt == 0 &&  nt < thmax){

      int xs = nt*nz*ny/thmax * iam;
      int xe = nt*nz*ny/thmax * (iam + 1);
      int ys = nt*nz/thmax * iam ;
      int i2=ys*nxs;
      const scs_t *tmp_PF_iyf = tmp + ys * nxs*ny + (ny-1)*nxs;
      const scs_t *tmp_PF_iyb = tmp + ys * nxs*ny;
      const pglus_t gy_PF_i0 = gy+ys*nxs*ny+(ny-1)*nxs;
      for (int ix = xs; ix < xe; ++ix) {
        int ixb = ix*nxs;
        int ixf = ix*nxs+nxs-1;
        if ( ix+_PFI != xe) {
#if VLENS == 16
#elif VLENS == 8
          __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[0]),0,1);
          __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[64]),0,1);
          __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[128]),0,1);
          __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[143]),0,1);
          __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[0]),1,1);
          __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[64]),1,1);
          __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[128]),1,1);
          __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[0]),1,1);
          __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[64]),1,1);
          __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[128]),1,1);
          __builtin_prefetch(&(xfs_recv[ix+_PFI].c[0][0][0]),0,1);
          __builtin_prefetch(&(xbs_recv[ix+_PFI].c[0][0][0]),0,1);
#endif
        } else{
#if VLENS == 16
#elif VLENS == 8
          __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[0]),0,1);
          __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[64]),0,1);
          __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[128]),0,1);
          __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[143]),0,1);
          __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[0]),1,1);
          __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[64]),1,1);
          __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[128]),1,1);
          __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[0]),1,1);
          __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[64]),1,1);
          __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[128]),1,1);
          __builtin_prefetch(&(yfs_recv[i2].cv[0][0][0]),0,1);
          __builtin_prefetch(&(yfs_recv[i2].cv[2][0][0]),0,1);
          __builtin_prefetch(&(ybs_recv[i2].cv[0][0][0]),0,1);
          __builtin_prefetch(&(ybs_recv[i2].cv[2][0][0]),0,1);
#endif
        }
        for(int x = 0; x < nxs; ++x) {
          for(int cs = 0; cs < 24; cs++) {
            for(int j = 0; j < VLENS; j++) {
              tmp[ixb+x].ccs[cs].v[j] = 0.0f;
            }
          }
        }
        projscs1_t *xbs_recvi = xbs_recv + ix;
        __mult_x_back_pst_3_(*(tmp+ixb),(*xbs_recvi));
        __attribute__((aligned(_ALIGN_SIZE))) projscs1_t ua;
        __mult_u_y_3_(ua,(*(xfs_recv + ix)),(*(gx + ixf)));
        __mult_x_forw_pst_3_(*(tmp+ixf),ua);
      }

      int ye = nt*nz/thmax * (iam  + 1);
      int tyxs = nt*ny*nxs/thmax * iam;
      int t = tyxs/(nxs*ny);
      int zs = tyxs%(nxs*ny);
      int ze = zs+nt*ny*nxs/thmax;
      int i4=nxs*ny*t+zs;
      const pglus_t gz_PF_i0 = gz+zs+nxs*ny*nz*t+(nz-1)*nxs*ny;
      for (int i1 = ys; i1 < ye; ++i1) {
        for (int x = 0; x < nxs; ++x,++i2) {
          int iyf = i1 * nxs*ny + (ny-1)*nxs +x;
          int iyb = i1 * nxs*ny + x;
          if ( x !=nxs-1) {
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&(yfs_recv[i2+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(yfs_recv[i2+_PFI].cv[2][0][0]),0,1);
            __builtin_prefetch(&(ybs_recv[i2+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(ybs_recv[i2+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&(zfs_recv[i4].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zfs_recv[i4].cv[2][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4].cv[2][0][0]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(yfs_recv + i2)),(*(gy + iyf)));
          __mult_y_forw_pst_(*(tmp+iyf),ua);
          projscs_t *ybs_recvi = ybs_recv + i2;
          __mult_y_back_pst_(*(tmp+iyb),(*ybs_recvi));

        }
      }
      //Z
      int ts = nz*ny*nxs/thmax * iam;
      const pglus_t gt_PF_i0 = gt + ts + nz*ny*nxs;
      const pclvs_t cl_PF_itf   = clvs + ts + nz*ny*nxs + vols*(*idomain);
      const pclvs_t cl_PF_itb   = clvs + ts + vols*(*idomain);
      const scs_t *out_PF_itf = out + ts + nz*ny*nxs;
      const scs_t *out_PF_itb = out + ts;
#pragma omp barrier
      if((ny*nxs) % (ze-zs) == 0){
        for (int i1 = zs; i1 < ze; ++i1,++i4) {
          int izf = i1+nxs*ny*nz*t+(nz-1)*nxs*ny;
          int izb = i1+nxs*ny*nz*t;
          if ( i1+_PFI !=ze) {
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[2][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(tfs_recv[ts].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tfs_recv[ts].cv[2][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[ts].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[ts].cv[2][0][0]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[128]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + izf)));
          __mult_z_forw_pst_(*(tmp+izf),ua);
          projscs_t *zbs_recvi = zbs_recv + i4;
          __mult_z_back_pst_(*(tmp+izb),(*zbs_recvi));
        }
      }else{
        int iz=0;
        for (int i1 = zs; i1 < ze; ++i1,++i4) {
          int izf = i1+iz+nxs*ny*nz*t+(nz-1)*nxs*ny;
          int izb = i1+iz+nxs*ny*nz*t;
          int _PFIZ = _PFI;
          if(i1%(nxs*ny)==(nxs*ny-1) ){
             _PFIZ = nxs*ny*nz - nxs*ny+1;
          }
          if ( i1+_PFI !=ze) {
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[2][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(tfs_recv[ts].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tfs_recv[ts].cv[2][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[ts].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[ts].cv[2][0][0]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[128]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + izf)));
          __mult_z_forw_pst_(*(tmp+izf),ua);
          projscs_t *zbs_recvi = zbs_recv + i4;
          __mult_z_back_pst_(*(tmp+izb),(*zbs_recvi));
          if(i1%(nxs*ny) == (nxs*ny-1) ){
            iz = -nxs*ny;
            t++;
          }
        }
      }
      
      //T
      int te = nz*ny*nxs/thmax * (iam + 1);
      if(nt ==2){
#pragma omp barrier
        for (int i1 = ts; i1 < te; ++i1) {
          int itf = i1 + nz*ny*nxs;
          if ( i1+_PFI !=te) {
#if VLENS == 16
#elif VLENS == 8
            __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[143]),0,1);
            __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[0]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[64]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[128]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[192]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[256]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[320]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[383]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[448]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[512]),0,1);
            __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[575]),0,1);
            __builtin_prefetch(&(tfs_recv[i1+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tfs_recv[i1+_PFI].cv[2][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[i1+_PFI].cv[0][0][0]),0,1);
            __builtin_prefetch(&(tbs_recv[i1+_PFI].cv[2][0][0]),0,1);
            __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[128]),1,1);
            __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[0]),1,1);
            __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[64]),1,1);
            __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[128]),1,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          float tbc_fwd = ((float)fbc[3][0])*0.5f;
          __mult_u_y_(ua,(*(tfs_recv + i1)),(*(gt + itf)));
          __mult_t_forw_pst_bc_(*(tmp+itf),ua,tbc_fwd);
          __mult_clvs( tmp[itf].cv, clvs[itf + vols*(*idomain)].cv);
#pragma loop norecurrence
          for (int c = 0; c < 3; c++) {
            for (int s = 0; s < 4; s++) {
              for (int ri = 0; ri < 2; ri++) {
                for (int j = 0; j < VLENS; j++) {
                  ((float(*)[4][2][VLENS])((out + itf)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+itf)->c)[c][s][ri][j] * factor;
                }
              }
            }
          }
          projscs_t *tbs_recvi = tbs_recv + i1;
          float tbc_bwd = ((float)fbc[3][1])*0.5f;
          __mult_t_back_pst_bc_(*(tmp+i1),(*tbs_recvi),tbc_bwd);
          __mult_clvs( tmp[i1].cv, clvs[i1 + vols*(*idomain)].cv);
#pragma loop norecurrence
          for (int c = 0; c < 3; c++) {
            for (int s = 0; s < 4; s++) {
              for (int ri = 0; ri < 2; ri++) {
                for (int j = 0; j < VLENS; j++) {
                  ((float(*)[4][2][VLENS])((out + i1)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+i1)->c)[c][s][ri][j] * factor;
                }
              }
            }
          }
        }
      }else{
#pragma omp barrier
        for (int i1 = ts; i1 < te; ++i1) {
          int i0 = i1 + nz*ny*nxs*(nt-1);
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          float tbc_fwd = ((float)fbc[3][0])*0.5f;
          __mult_u_y_(ua,(*(tfs_recv + i1)),(*(gt + i0)));
          __mult_t_forw_pst_bc_(*(tmp+i0),ua,tbc_fwd);
          i0 = i1;
          projscs_t *tbs_recvi = tbs_recv + i1;
          float tbc_bwd = ((float)fbc[3][1])*0.5f;
          __mult_t_back_pst_bc_(*(tmp+i0),(*tbs_recvi),tbc_bwd);
        }
        int cls = nz*ny*nxs*nt/thmax * iam;
        int cle = nz*ny*nxs*nt/thmax * (iam + 1);
#pragma omp barrier
        for (int i0 = cls; i0 < cle; ++i0) {
          __mult_clvs( tmp[i0].cv, clvs[i0 + vols*(*idomain)].cv);
#pragma loop norecurrence
          for (int c = 0; c < 3; c++) {
            for (int s = 0; s < 4; s++) {
              for (int ri = 0; ri < 2; ri++) {
                for (int j = 0; j < VLENS; j++) {
                  ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+i0)->c)[c][s][ri][j] * factor;
                }
              }
            }
          }
        }
      }
#pragma omp barrier
    }else{
#pragma omp  for collapse(4)
      for (int t = 0; t < nt; t++) {
        for (int z = 0; z < nz; z++) {
          for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nxs; x++) {

              if (x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x==0 || y == 0 || z == 0 || t == 0 ) {

                int it = nxs*ny*nz*t;
                int iz = nxs*ny*z;
                int iy = nxs*y;
                int i0 = x + iy + iz + it;

                scs_t tmp __attribute__((aligned(_ALIGN_SIZE))) ;

                for(int c = 0; c < 3; c++) {
                  for(int s = 0; s < 4; s++) {
                    for(int ri = 0; ri < 2; ri++) {
                      for(int j = 0; j < VLENS; j++) {
                        tmp.c[c][s][ri][j] = 0.0f;
                      }
                    }
                  }
                }

                //
                // X-forward
                //
                if ( x == nxs-1 ) { _mult_forw_x_recv_; }

                //
                // X-backward
                //
                if ( x == 0 ) { _mult_back_x_recv_; }

                //
                // Y-forward 
                //
                if ( y == ny-1 ) { _mult_forw_y_recv_; }

                //
                // Y-backward
                //
                if ( y == 0 ) { _mult_back_y_recv_; }

                //
                // Z-forward
                //
                if ( z == nz-1 ) { _mult_forw_z_recv_; }

                //
                // Z-backward
                //
                if ( z == 0 ) { _mult_back_z_recv_; }

                //
                // T-forward
                //
                if ( t == nt-1 ) { _mult_forw_t_recv_; }

                //
                // T-backward
                //
                if ( t == 0 ) { _mult_back_t_recv_; }

                __mult_clvs( tmp.cv, clvs[i0 + vols*(*idomain)].cv);

#pragma loop norecurrence
                for (int c = 0; c < 3; c++) {
                  for (int s = 0; s < 4; s++) {
                    for (int ri = 0; ri < 2; ri++) {
                      for (int j = 0; j < VLENS; j++) {
                        ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp.c))[c][s][ri][j] * factor;
                      }
                    }
                  }
                }
              } // if edge sites
            }
          }
        }
      }
    }

    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_;
    _COMM_TIC_;
#ifdef RDC
#pragma omp master
#else
#pragma omp single
#endif
    {
      xbound_send_waitall(4);
    }
#ifdef RDC
#pragma omp barrier
#endif
    _COMM_TOC_;
  }
#undef factor 
#undef nxs
#undef ny 
#undef nz 
#undef nt 
#undef vols 
  void ddd_out_pos_s_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
#pragma omp parallel
  {
    ddd_out_pos_s_noprl_(out, in, idomain, factor);
  }
}

#ifndef RDC
  void ddd_out_pos_s_noprl_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {

    const float lfactor = factor;
    const int lnxs = nxs, lny = ny, lnz = nz, lnt = nt;
    const int lvols = vols;
#define factor lfactor
#define nxs lnxs
#define ny lny
#define nz lnz
#define nt lnt
#define vols lvols
    const pglus_t __restrict__ gx = &glus[vols*0 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gy = &glus[vols*1 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gz = &glus[vols*2 + NDIM*vols*(*idomain)];
    const pglus_t __restrict__ gt = &glus[vols*3 + NDIM*vols*(*idomain)];
    static scs_t *tmp;
#pragma omp single
    {
    if( tmp==0)  posix_memalign((void **)&tmp,256,sizeof(scs_t)*vols);
    }

#pragma omp single
    {
      xbound_wait(0,4);
      xbound_wait(1,4);
      xbound_wait(2,4);
      xbound_wait(3,4);
      xbound_wait(4,4);
      xbound_wait(5,4);
      xbound_wait(6,4);
      xbound_wait(7,4);
    }
    int nmodx,nmody,nmodz,nmodt;
    nmodx= (ny * nz * nt)%thmax;
    nmody= ( nz * nt)%thmax;
    nmodz= ( ny * nt)%thmax;
    nmodt= ( ny * nz)%thmax;
    if (nmodx == 0 && nmody == 0 && nmodz == 0 && nmodt == 0 &&  nt < thmax){
//  if (iam==0)printf("fj dbg\n");

    int xs = nt*nz*ny/thmax * iam;
    int xe = nt*nz*ny/thmax * (iam + 1);
    int ys = nt*nz/thmax * iam ;
    int i2=ys*nxs;
    const scs_t *tmp_PF_iyf = tmp + ys * nxs*ny + (ny-1)*nxs;
    const scs_t *tmp_PF_iyb = tmp + ys * nxs*ny;
    const pglus_t gy_PF_i0 = gy+ys*nxs*ny+(ny-1)*nxs;
    for (int ix = xs; ix < xe; ++ix) {
         int ixb = ix*nxs;
         int ixf = ix*nxs+nxs-1;
       if ( ix+_PFI != xe) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gx+ixf+_PFI*nxs)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+ixf+_PFI*nxs)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+ixb+_PFI*nxs)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(xfs_recv[ix+_PFI].c[0][0][0]),0,1);
      __builtin_prefetch(&(xbs_recv[ix+_PFI].c[0][0][0]),0,1);
#endif
            } else{
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gy_PF_i0)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp_PF_iyf)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp_PF_iyb)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(yfs_recv[i2].cv[0][0][0]),0,1);
      __builtin_prefetch(&(yfs_recv[i2].cv[2][0][0]),0,1);
      __builtin_prefetch(&(ybs_recv[i2].cv[0][0][0]),0,1);
      __builtin_prefetch(&(ybs_recv[i2].cv[2][0][0]),0,1);
#endif
      }
         for(int x = 0; x < nxs; ++x) {
           for(int cs = 0; cs < 24; cs++) {
             for(int j = 0; j < VLENS; j++) {
               tmp[ixb+x].ccs[cs].v[j] = 0.0f;
             }
           }
         }
         projscs1_t *xbs_recvi = xbs_recv + ix;
         __mult_x_back_pst_3_(*(tmp+ixb),(*xbs_recvi));
          __attribute__((aligned(_ALIGN_SIZE))) projscs1_t ua;
          __mult_u_y_3_(ua,(*(xfs_recv + ix)),(*(gx + ixf)));
          __mult_x_forw_pst_3_(*(tmp+ixf),ua);
      }

    int ye = nt*nz/thmax * (iam  + 1);
    int tyxs = nt*ny*nxs/thmax * iam;
    int t = tyxs/(nxs*ny);
    int zs = tyxs%(nxs*ny);
    int ze = zs+nt*ny*nxs/thmax;
//    int i4=nxs*ny*t+zs;
    int i4=nxs*ny*t+zs;
    const pglus_t gz_PF_i0 = gz+zs+nxs*ny*nz*t+(nz-1)*nxs*ny;
    for (int i1 = ys; i1 < ye; ++i1) {
      for (int x = 0; x < nxs; ++x,++i2) {
          int iyf = i1 * nxs*ny + (ny-1)*nxs +x;
          int iyb = i1 * nxs*ny + x;
          if ( x !=nxs-1) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gy+iyf+_PFI)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+iyf+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+iyb+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(yfs_recv[i2+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(yfs_recv[i2+_PFI].cv[2][0][0]),0,1);
      __builtin_prefetch(&(ybs_recv[i2+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(ybs_recv[i2+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gz_PF_i0)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&(zfs_recv[i4].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zfs_recv[i4].cv[2][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4].cv[2][0][0]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(yfs_recv + i2)),(*(gy + iyf)));
          __mult_y_forw_pst_(*(tmp+iyf),ua);
          projscs_t *ybs_recvi = ybs_recv + i2;
          __mult_y_back_pst_(*(tmp+iyb),(*ybs_recvi));

      }
      }
      //Z
//    int ze = zs+nt*ny*nxs/thmax;
    int ts = nz*ny*nxs/thmax * iam;
    const pglus_t gt_PF_i0 = gt + ts + nz*ny*nxs;
    const pclvs_t cl_PF_itf   = clvs + ts + nz*ny*nxs + vols*(*idomain);
    const pclvs_t cl_PF_itb   = clvs + ts + vols*(*idomain);
    const scs_t *out_PF_itf = out + ts + nz*ny*nxs;
    const scs_t *out_PF_itb = out + ts;
#pragma omp barrier
//  if (iam==0)printf("fj dbg2\n");
    if((ny*nxs) % (ze-zs) == 0){
//  if (iam==0)printf("fj dbg\n");
      for (int i1 = zs; i1 < ze; ++i1,++i4) {
          int izf = i1+nxs*ny*nz*t+(nz-1)*nxs*ny;
          int izb = i1+nxs*ny*nz*t;
          if ( i1+_PFI !=ze) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFI)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+izf+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[2][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(tfs_recv[ts].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tfs_recv[ts].cv[2][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[ts].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[ts].cv[2][0][0]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[128]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + izf)));
          __mult_z_forw_pst_(*(tmp+izf),ua);
          projscs_t *zbs_recvi = zbs_recv + i4;
          __mult_z_back_pst_(*(tmp+izb),(*zbs_recvi));
      }
      }else{
//  if (iam==0)printf("fj dbg3\n");
          int iz=0;
      for (int i1 = zs; i1 < ze; ++i1,++i4) {
          int izf = i1+iz+nxs*ny*nz*t+(nz-1)*nxs*ny;
          int izb = i1+iz+nxs*ny*nz*t;
          int _PFIZ = _PFI;
          if(i1%(nxs*ny)==(nxs*ny-1) ){
             _PFIZ = nxs*ny*nz - nxs*ny+1;
          }
          if ( i1+_PFI !=ze) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gz+izf+_PFIZ)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+izf+_PFIZ)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+izb+_PFIZ)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zfs_recv[i4+_PFI].cv[2][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(zbs_recv[i4+_PFI].cv[2][0][0]),0,1);
#endif
          } else{
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gt_PF_i0)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((cl_PF_itf )->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((cl_PF_itf)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((cl_PF_itb )->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((cl_PF_itb)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(tfs_recv[ts].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tfs_recv[ts].cv[2][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[ts].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[ts].cv[2][0][0]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( out_PF_itf)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( out_PF_itb)->c_prefetch)[128]),0,1);
#endif
          }
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
          __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + izf)));
          __mult_z_forw_pst_(*(tmp+izf),ua);
          projscs_t *zbs_recvi = zbs_recv + i4;
          __mult_z_back_pst_(*(tmp+izb),(*zbs_recvi));
          if(i1%(nxs*ny) == (nxs*ny-1) ){
            iz = -nxs*ny;
            t++;
          }
      }
      }
      
      //T
    int te = nz*ny*nxs/thmax * (iam + 1);
    if(nt ==2){
#pragma omp barrier
    for (int i1 = ts; i1 < te; ++i1) {
       int itf = i1 + nz*ny*nxs;
          if ( i1+_PFI !=te) {
#if VLENS == 16
#elif VLENS == 8
      __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&((( gt+itf+_PFI)->c_prefetch)[143]),0,1);
      __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+itf+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( tmp+i1+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(((clvs+itf+vols*(*idomain)+_PFI)->c_prefetch)[0]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[64]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[128]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[192]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[256]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[320]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[383]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[448]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[512]),0,1);
      __builtin_prefetch(&(((clvs+i1+vols*(*idomain)+_PFI)->c_prefetch)[575]),0,1);
      __builtin_prefetch(&(tfs_recv[i1+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tfs_recv[i1+_PFI].cv[2][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[i1+_PFI].cv[0][0][0]),0,1);
      __builtin_prefetch(&(tbs_recv[i1+_PFI].cv[2][0][0]),0,1);
      __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( out+itf+_PFI)->c_prefetch)[128]),1,1);
      __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[0]),1,1);
      __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[64]),1,1);
      __builtin_prefetch(&((( out+i1+_PFI)->c_prefetch)[128]),1,1);
#endif
          }
//       int ibf = i1;
       __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
       float tbc_fwd = ((float)fbc[3][0])*0.5f;
       __mult_u_y_(ua,(*(tfs_recv + i1)),(*(gt + itf)));
       __mult_t_forw_pst_bc_(*(tmp+itf),ua,tbc_fwd);
       __mult_clvs( tmp[itf].cv, clvs[itf + vols*(*idomain)].cv);
#pragma loop norecurrence
      for (int c = 0; c < 3; c++) {
      for (int s = 0; s < 4; s++) {
      for (int ri = 0; ri < 2; ri++) {
      for (int j = 0; j < VLENS; j++) {
 	      ((float(*)[4][2][VLENS])((out + itf)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+itf)->c)[c][s][ri][j] * factor;
       }
           }
         }
         }
       projscs_t *tbs_recvi = tbs_recv + i1;
       float tbc_bwd = ((float)fbc[3][1])*0.5f;
       __mult_t_back_pst_bc_(*(tmp+i1),(*tbs_recvi),tbc_bwd);
       __mult_clvs( tmp[i1].cv, clvs[i1 + vols*(*idomain)].cv);
#pragma loop norecurrence
      for (int c = 0; c < 3; c++) {
      for (int s = 0; s < 4; s++) {
      for (int ri = 0; ri < 2; ri++) {
      for (int j = 0; j < VLENS; j++) {
 	      ((float(*)[4][2][VLENS])((out + i1)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+i1)->c)[c][s][ri][j] * factor;
       }
           }
         }
         }
      }
    }else{
#pragma omp barrier
    for (int i1 = ts; i1 < te; ++i1) {
       int i0 = i1 + nz*ny*nxs*(nt-1);
       __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;
       float tbc_fwd = ((float)fbc[3][0])*0.5f;
       __mult_u_y_(ua,(*(tfs_recv + i1)),(*(gt + i0)));
       __mult_t_forw_pst_bc_(*(tmp+i0),ua,tbc_fwd);
//          if ( t == 0 ) { _mult_back_t_recv_; }
       i0 = i1;
       projscs_t *tbs_recvi = tbs_recv + i1;
       float tbc_bwd = ((float)fbc[3][1])*0.5f;
       __mult_t_back_pst_bc_(*(tmp+i0),(*tbs_recvi),tbc_bwd);
      }
    int cls = nz*ny*nxs*nt/thmax * iam;
    int cle = nz*ny*nxs*nt/thmax * (iam + 1);
#pragma omp barrier
    for (int i0 = cls; i0 < cle; ++i0) {
       __mult_clvs( tmp[i0].cv, clvs[i0 + vols*(*idomain)].cv);
#pragma loop norecurrence
      for (int c = 0; c < 3; c++) {
      for (int s = 0; s < 4; s++) {
      for (int ri = 0; ri < 2; ri++) {
      for (int j = 0; j < VLENS; j++) {
 	      ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp+i0)->c)[c][s][ri][j] * factor;
      }
      }
      }
      }
      }
      }
#pragma omp barrier
   }else{
#pragma omp  for collapse(4)
    for (int t = 0; t < nt; t++) {
    for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nxs; x++) {

        if (x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x==0 || y == 0 || z == 0 || t == 0 ) {

          int it = nxs*ny*nz*t;
          int iz = nxs*ny*z;
          int iy = nxs*y;
          int i0 = x + iy + iz + it;

          scs_t tmp __attribute__((aligned(_ALIGN_SIZE))) ;

          for(int c = 0; c < 3; c++) {
          for(int s = 0; s < 4; s++) {
            for(int ri = 0; ri < 2; ri++) {
              for(int j = 0; j < VLENS; j++) {
                tmp.c[c][s][ri][j] = 0.0f;
              }
            }
          }
          }

          //
          // X-forward
          //
          if ( x == nxs-1 ) { _mult_forw_x_recv_; }

          //
          // X-backward
          //
          if ( x == 0 ) { _mult_back_x_recv_; }

          //
          // Y-forward 
          //
          if ( y == ny-1 ) { _mult_forw_y_recv_; }

          //
          // Y-backward
          //
          if ( y == 0 ) { _mult_back_y_recv_; }

          //
          // Z-forward
          //
          if ( z == nz-1 ) { _mult_forw_z_recv_; }

          //
          // Z-backward
          //
          if ( z == 0 ) { _mult_back_z_recv_; }

          //
          // T-forward
          //
          if ( t == nt-1 ) { _mult_forw_t_recv_; }

          //
          // T-backward
          //
          if ( t == 0 ) { _mult_back_t_recv_; }

          __mult_clvs( tmp.cv, clvs[i0 + vols*(*idomain)].cv);

#pragma loop norecurrence
          for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int ri = 0; ri < 2; ri++) {
	      for (int j = 0; j < VLENS; j++) {
		((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp.c))[c][s][ri][j] * factor;
	      }
            }
          }
          }
        } // if edge sites
      }
      }
      }
      }
    }

#pragma omp single
    {
      xbound_send_waitall(4);
    }
  }
#endif
//RDC
#undef factor 
#undef nxs
#undef ny 
#undef nz 
#undef nt 
#undef vols 
  void ddd_out_pos_s_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
#ifndef RDC
#pragma omp parallel
{
    ddd_out_pos_s_noprl_no_timer_(out, in, idomain, factor);
}
#endif
}


#ifdef __cplusplus
}
#endif
