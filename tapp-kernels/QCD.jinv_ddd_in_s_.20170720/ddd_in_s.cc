//#define PREFETCH
#define _PFI 1
#include "qws.h"
#include "clover_s.h"
#include "mult_all.h"
#include "prefetch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(_ALIGN_SIZE))) pglus_t glus;
  extern __attribute__((aligned(_ALIGN_SIZE))) pclvs_t clvs;
  extern double kappa, mkappa;
  extern int nxs, ny, nz, nt;
  extern int vols;
  extern int thmax;
  extern int iam;
#pragma omp threadprivate(iam)

  //---------------------------------------------------------------------------------------- mult D in bulk
  void ddd_in_s_noprl(scs_t* __restrict__ out, const scs_t* __restrict__ in, const int* DEO) {

    const double lkappa = kappa, lmkappa = mkappa;
    const int lnxs = nxs, lny = ny, lnz = nz, lnt = nt;
    const int lvols = vols;

#define kappa lkappa
#define mkappa lmkappa
#define nxs lnxs
#define ny lny
#define nz lnz
#define nt lnt
#define vols lvols

    const __restrict__ pglus_t gx = &glus[vols*0 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gy = &glus[vols*1 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gz = &glus[vols*2 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gt = &glus[vols*3 + NDIM*vols*(*DEO)];
    const __restrict__ pclvs_t cl = &clvs[              vols*(*DEO)];

    //
    // X
    //
      int ntzs,ntze;
      float ff[3][4][2][VLENS*2] __attribute__((aligned(4)));
      float ff_in[3][4][2][VLENS*2] __attribute__((aligned(4)));
      float ff_g [3][3][2][VLENS*2] __attribute__((aligned(4)));
#ifdef _OPENMP
      int nsize,nmod;
      int thmod,PF_iyf,PF_t,PF_i0;
//    iam = omp_get_thread_num();
      nmod = (nt*nz)%thmax;
      thmod = (nmod > iam ? 1:0 );
      nsize=nt*nz/thmax ;
      ntzs=nsize*iam + thmod*iam + nmod *(1-thmod);
      ntze=ntzs+nsize+thmod;
      PF_t= ntzs / nz;
      PF_iyf=nxs+nxs*ny*(ntzs%nz)+nxs*ny*nz*PF_t;
      PF_i0=nxs*ny*(ntzs%nz)+nxs*ny*nz*PF_t;
      const scs_t *in_PF_iyf = in + PF_iyf;
      scs_t *out_PF_i0 = out+PF_i0;
      const pglus_t gy_PF_i0 = gy+PF_i0;
      int te= ntze / nz;
      int ze= ntze % nz;
      int i0e=  nxs*ny*ze + nxs*ny*nz*te;
#else
      int PF_iyf,PF_t,PF_i0;
      ntzs=0; ntze=nt*nz; 
      PF_iyf=nxs;
      PF_i0=0;
      int i0e=  nxs*ny*nz*nt;
#endif
      int i0 = (ntzs/nz)*nxs*ny*nz + (ntzs%nz)*nxs*ny;
      for(int tz = ntzs; tz < ntze; ++tz) {
        for(int y = 0; y < ny; ++y) {

          for(int x = 0; x < nxs; ++x, ++i0) {

            int ixf = i0+1;
            int ixb = i0-1;
            int pfx;
            if(x+_PFI < nxs){
              pfx=x+_PFI;
            }else{
              pfx=0;
            }

#if VLENS == 16
            if ( i0+_PFI != i0e) {
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[383]),1,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[288]),0,1);
              __builtin_prefetch(&(ff[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff[0][3][0][0]),1,1);
              __builtin_prefetch(&(ff[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff[1][3][0][0]),1,1);
              __builtin_prefetch(&(ff[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_in[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_g[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][1][VLENS*2-1]),1,1);
              if (pfx != nxs-1){
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[383]),0,1);
              }
              if (pfx != 0){
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[287]),0,1);
              }
            } else{
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[383]),1,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[287]),0,1);
            }
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
            if ( i0+_PFI != i0e) {
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[96]) ,0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[96]) ,1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[96]) ,0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[143]),0,1);
              __builtin_prefetch(&(ff[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff[0][3][0][0]),1,1);
              __builtin_prefetch(&(ff[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff[1][3][0][0]),1,1);
              __builtin_prefetch(&(ff[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_in[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_g[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][1][VLENS*2-1]),1,1);
              if (pfx != nxs-1){
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[96]) ,0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[191]),0,1);
              }
              if (pfx != 0){
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[96]) ,0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[96]) ,0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[143]),0,1);
              }
            } else{
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[96]) ,0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[96]) ,1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[143]),0,1);
            }
#else
            if ( i0+_PFI != i0e) {
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[143]),0,1);
              __builtin_prefetch(&(ff[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_in[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_g[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][1][VLENS*2-1]),1,1);
              if (pfx != nxs-1){
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[191]),0,1);
              }
              if (pfx != 0){
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[64]) ,0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[143]),0,1);
              }
            } else{
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[64]) ,0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]) ,1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[143]),0,1);
            }
#endif
#elif VLENS == 4
            if ( i0+_PFI != i0e) {
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[32]) ,0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in+i0+_PFI)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[32]) ,1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[95]),1,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[32]) ,0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gx+i0+_PFI)->c_prefetch)[71]),0,1);
              __builtin_prefetch(&(ff[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_in[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[1][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_in[2][3][1][VLENS*2-1]),1,1);
              __builtin_prefetch(&(ff_g[0][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[0][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[1][1][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][0][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][0][0]),1,1);
              __builtin_prefetch(&(ff_g[2][2][1][VLENS*2-1]),1,1);
              if (pfx != nxs-1){
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[32]) ,0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+ixf+_PFI)->c_prefetch)[95]),0,1);
              }
              if (pfx != 0){
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[32]) ,0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+ixb+_PFI)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[32]) ,0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gx+ixb+_PFI)->c_prefetch)[71]),0,1);
              }
            } else{
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[32]) ,0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_iyf)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]) ,1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[95]),1,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gy_PF_i0 )->c_prefetch)[71]),0,1);
            }
#endif

            //
            // X-forward
            //
            {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              if ( x != nxs-1 ) {

                for(int c = 0; c < 3; ++c) {
                  for(int s = 0; s < 4; ++s) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j]; 
                        ff[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+ixf)->c))[c][s][ri][j];
                      }
                    }
                  }
                }

              } else {

                for(int c = 0; c < 3; ++c) {
                  for(int s = 0; s < 4; ++s) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j]; 
                        ff[c][s][ri][j+VLENS] = 0.0f;
                      }
                    }
                  }
                }

              }

              __mult_x_forw_pre_2_(a,ff); // with forward simd-site-shift
              const __restrict__ pglus_t gx_i0 = gx+i0;
              __mult_u_y_(ua,a,(*(gx_i0)));
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_x_forw_pst_((*out_i0),ua);
            }

            //
            // X-backward
            //
            {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

              if ( x != 0 ) {

                for(int c = 0; c < 3; ++c) {
                  for(int s = 0; s < 4; ++s) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff_in[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+ixb)->c))[c][s][ri][j];
                        ff_in[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j];
                      }
                    }
                  }
                }

                //
                // load link U.  WITHOUT taking Hermitian-conjugate
                //
                for (int c2 = 0; c2 < 3; ++c2) {
                  for (int c1 = 0; c1 < 3; ++c1) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff_g[c2][c1][ri][j      ] = ((float(*)[3][2][VLENS])((gx+ixb)->c))[c2][c1][ri][j];
                        ff_g[c2][c1][ri][j+VLENS] = ((float(*)[3][2][VLENS])((gx+i0)->c))[c2][c1][ri][j];
                      }
                    }
                  }
                }

              } else {

                for (int c = 0; c < 3; ++c) {
                  for (int s = 0; s < 4; ++s) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff_in[c][s][ri][j      ] = 0.0f;
                        ff_in[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j];
                      }
                    }
                  }
                }

                //
                // load link U.  WITHOUT taking Hermitian-conjugate
                //
                for (int c2 = 0; c2 < 3; ++c2) {
                  for (int c1 = 0; c1 < 3; ++c1) {
                    for(int ri = 0; ri < 2; ++ri) { 
                      for(int j = 0; j < VLENS; ++j) { 
                        ff_g[c2][c1][ri][j      ] = 0.0f;
                        ff_g[c2][c1][ri][j+VLENS] = ((float(*)[3][2][VLENS])((gx+i0)->c))[c2][c1][ri][j];
                      }
                    }
                  }
                }

              }

              __mult_x_back_pre_2_(a,ff_in);// with backward simd-site-shift
              __mult_udag_y_2_(ua,a,ff_g);  // with backward simd-site-shift
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_x_back_pst_((*out_i0),ua);

            }
          }
        }
      }
      //
      // Y
      //
#ifdef _OPENMP
      nmod = (nt*ny)%thmax;
      thmod = (nmod > iam ? 1:0 );
      nsize=nt*ny/thmax ;
      int ntys=nsize*iam + thmod*iam + nmod *(1-thmod);
      int ntye=ntys+nsize+thmod;
      PF_t=ntys/ny;
      int PF_y=ntys%ny;
      int PF_izf = nxs*PF_y + nxs*ny*nz*PF_t + nxs*ny;
      PF_i0= nxs*PF_y + nxs*ny*nz*PF_t;
      out_PF_i0=out+PF_i0;
      const scs_t *in_PF_izf = in + PF_izf;
      const pglus_t gz_PF_i0 = gz+PF_i0;
#else
      int ntys=0; 
      int ntye=nt*ny;
#endif
      i0 = (ntzs/nz)*nxs*ny*nz + (ntzs%nz)*nxs*ny;
      for(int tz = ntzs; tz < ntze; ++tz) {
        for(int y = 0; y < ny; ++y) {

          //      for(int x = 0; x < nxs; ++x) {
          for(int x = 0; x < nxs; ++x, ++i0) {

            //        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
            int iyf = i0 + nxs;
            int iyb = i0 - nxs;
            int pfy;
            if(x+_PFI < nxs){
              pfy=y;
            }else if (y != ny-1){
              pfy=y+1;
            }else{
              pfy=0;
            }
#if VLENS == 16
            if ( i0+_PFI != i0e) {
              if (pfy != ny-1){
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[287]),0,1);
              }
              if (pfy != 0){
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[287]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[383]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[287]),0,1);
            }
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
            if ( i0+_PFI != i0e) {
              if (pfy != ny-1){
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[143]),0,1);
              }
              if (pfy != 0){
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[143]),0,1);
            }
#else
            if ( i0+_PFI != i0e) {
              if (pfy != ny-1){
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[143]),0,1);
              }
              if (pfy != 0){
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[143]),0,1);
            }
#endif
#elif VLENS == 4
            if ( i0+_PFI != i0e) {
              if (pfy != ny-1){
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyf+_PFI)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+i0+_PFI )->c_prefetch)[71]),0,1);
              }
              if (pfy != 0){
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+iyb+_PFI)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gy+iyb+_PFI)->c_prefetch)[71]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFI)->c_prefetch)[95]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_izf)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gz_PF_i0 )->c_prefetch)[71]),0,1);
            }
#endif
            //
            // Y-forward
            //
            if ( y != ny-1 ) {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gy_i0 = gy+i0;
              const scs_t* __restrict__ in_iyf = in + iyf;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_y_forw_pre_(a,(*(in_iyf)));
              __mult_u_y_(ua,a,(*(gy_i0)));
              __mult_y_forw_pst_((*out_i0),ua);
            } 
            //
            // Y-backward
            //
            if ( y != 0 ) {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gy_iyb = gy+iyb;
              const scs_t* __restrict__ in_iyb = in + iyb;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_y_back_pre_(a,(*(in_iyb)));
              __mult_udag_y_(ua,a,(*(gy_iyb)));
              __mult_y_back_pst_((*out_i0),ua);
            } 
          }
        }
      }

#pragma omp barrier
#if VLENS == 16
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[192]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[256]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[320]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[383]),1,1);
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[96]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[160]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
#else
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
#endif
#elif VLENS == 4
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[95]),1,1);
#endif
      //
      // Z
      //
#ifdef _OPENMP
      nmod = (nz*ny)%thmax;
      thmod = (nmod > iam ? 1:0 );
      nsize=nz*ny/thmax ;
      int nzys=nsize*iam + thmod*iam + nmod *(1-thmod);
      int nzye=nzys+nsize+thmod;
      int PF_z= nzys / ny;
      int PF_itf=nxs*ny*PF_z + nxs*ny*nz + nxs*(nzys%ny);
      PF_i0=nxs*ny*PF_z + nxs*(nzys%ny);
      te= ntye / ny; 
      int ye= ntye % ny;
      i0e=  nxs*ye + nxs*ny*nz*te;
      out_PF_i0=out+PF_i0;
      const scs_t *in_PF_itf = in + PF_itf;
      const pglus_t gt_PF_i0 = gt+PF_i0;
#else
      int nzys=0; 
      int nzye=nz*ny;
#endif
      i0 = (ntys/ny)*nxs*ny*nz + (ntys%ny)*nxs;
      for(int ty = ntys; ty < ntye; ++ty) {
        for(int z = 0; z < nz; ++z) {
          for(int x = 0; x < nxs; ++x, ++i0) {
            //        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
            int izf = i0 + nxs*ny;
            int izb = i0 - nxs*ny;
            int _PFIZ;
            int pfz;
            if(x+_PFI < nxs){
              _PFIZ=_PFI;
              pfz=z;
            }else if (z != nz-1){
              _PFIZ = nxs * ny -nxs+_PFI;
              pfz=z+1;
            }else{
              if(ty%ny != ny-1){
              _PFIZ = nxs-nxs*ny*(nz-1)-nxs+_PFI;
              }else{
//              _PFIZ = -nxs*ny*(nz-1)-nxs*(ny-1)+nxs*ny*nz-nxs+_PFI;
              _PFIZ = _PFI;
              }
              pfz=0;
            }
#if VLENS == 16
            if ( i0+_PFIZ != i0e) {
              if (pfz != nz-1){
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[287]),0,1);
              }
              if (pfz != 0){
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[287]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[383]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[287]),0,1);
            }
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
            if ( i0+_PFIZ != i0e) {
              if (pfz != nz-1){
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[143]),0,1);
              }
              if (pfz != 0){
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[143]),0,1);
            }
#else
            if ( i0+_PFIZ != i0e) {
              if (pfz != nz-1){
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[143]),0,1);
              }
              if (pfz != 0){
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[143]),0,1);
            }
#endif
#elif VLENS == 4
            if ( i0+_PFIZ != i0e) {
              if (pfz != nz-1){
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izf+_PFIZ)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+i0+_PFIZ )->c_prefetch)[71]),0,1);
              }
              if (pfz != 0){
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+izb+_PFIZ)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gz+izb+_PFIZ)->c_prefetch)[71]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIZ)->c_prefetch)[95]),1,1);
            } else{
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_itf)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gt_PF_i0 )->c_prefetch)[71]),0,1);
            }
#endif
            //
            // Z-forward
            //
            if ( z != nz-1 ) {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gz_i0 = gz+i0;
              const scs_t* __restrict__ in_izf = in + izf;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_z_forw_pre_(a,(*(in_izf)));
              __mult_u_y_(ua,a,(*(gz_i0)));
              __mult_z_forw_pst_((*out_i0),ua);
            }
            //
            // Z-backward
            //
            if ( z != 0 ) {
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gz_izb = gz+izb;
              const scs_t* __restrict__ in_izb = in + izb;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_z_back_pre_(a,(*(in_izb)));
              __mult_udag_y_(ua,a,(*(gz_izb)));
              __mult_z_back_pst_((*out_i0),ua);
            } 
          }
          i0=i0+nxs*ny-nxs;
        }
        //i0=i0+nxs-nxs*ny*nz;
        if(ty%ny != ny-1){
        i0=i0+nxs-nxs*ny*nz;
        }else{
        i0=i0-nxs*ny+nxs;
        }
      }

#pragma omp barrier
#if VLENS == 16
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[192]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[256]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[320]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[383]),1,1);
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[96]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[160]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
#else
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[128]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[191]),1,1);
#endif
#elif VLENS == 4
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[0]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[32]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[64]),1,1);
      __builtin_prefetch(&(((out_PF_i0 )->c_prefetch)[95]),1,1);
#endif

      //
      // T
      //
      ze= nzye / ny; 
      ye= nzye % ny;
      i0e=  nxs*ye + nxs*ny*ze;
      i0 = (nzys/ny)*nxs*ny + (nzys%ny)*nxs;
      const scs_t *in_PF_i0 = in + PF_i0;
      const pclvs_t cl_PF_i0 = cl+PF_i0;
      for(int zy = nzys; zy < nzye; ++zy) {
        for(int t = 0; t < nt; ++t) {
          for(int x = 0; x < nxs; ++x, ++i0) {
            int itf = i0 + nxs*ny*nz;
            int itb = i0 - nxs*ny*nz;
            int _PFIT;
            int pft;
            if(x+_PFI < nxs){
              _PFIT=_PFI;
              pft = t;
            }else if (t != nt-1){
              _PFIT =  nxs * ny * nz -  nxs + _PFI;
              pft = t+1;
            }else{
              _PFIT = _PFI - nxs * ny * nz* (nt-1) ;
              pft = 0;
            }
#if VLENS == 16
            if ( i0+_PFIT != i0e) {
              if (pft != nt-1){
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[287]),0,1);
              }
              if (pft != 0){
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[320]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[383]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[192]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[256]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[287]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[383]),1,1);
            } else{
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[383]),1,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[576]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[640]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[704]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[768]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[832]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[896]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[960]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[1024]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[1088]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[1151]),0,1);
            }
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
            if ( i0+_PFIT != i0e) {
              if (pft != nt-1){
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[143]),0,1);
              }
              if (pft != 0){
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[160]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[96]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[224]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[288]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[352]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[416]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[480]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[544]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[575]),0,1);
            }
#else
            if ( i0+_PFIT != i0e) {
              if (pft != nt-1){
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[143]),0,1);
              }
              if (pft != 0){
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[191]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[128]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[143]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[191]),1,1);
            } else{
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[575]),0,1);
            }
#endif
#elif VLENS == 4
            if ( i0+_PFIT != i0e) {
              if (pft != nt-1){
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itf+_PFIT)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+i0+_PFIT )->c_prefetch)[71]),0,1);
              }
              if (pft != 0){
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( in+itb+_PFIT)->c_prefetch)[95]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[0]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[32]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[64]),0,1);
                __builtin_prefetch(&((( gt+itb+_PFIT)->c_prefetch)[71]),0,1);
              }
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[95]),1,1);
            } else{
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[95]),1,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in_PF_i0)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[224]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl_PF_i0)->c_prefetch)[287]),0,1);
            }
#endif
            //
            // T-forward
            //
            if ( t != nt-1 ) { 
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gt_i0 = gt+i0;
              const scs_t* __restrict__ in_itf = in + itf;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_t_forw_pre_(a,(*(in_itf)));
              __mult_u_y_(ua,a,(*(gt_i0)));
              __mult_t_forw_pst_((*out_i0),ua);
            } 
            //
            // T-backward
            //
            if ( t != 0 ) { 
              __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
              const __restrict__ pglus_t gt_itb = gt + itb;
              const scs_t* __restrict__ in_itb = in + itb;
              scs_t* __restrict__ out_i0 = out+i0;
              __mult_t_back_pre_(a,(*(in_itb)));
              __mult_udag_y_(ua,a,(*(gt_itb)));
              __mult_t_back_pst_((*out_i0),ua);
            }
          }
          i0 = i0 + nxs*ny*nz - nxs;
        }
        i0 = i0 + nxs - nxs*ny*nz*nt;
      }
      //
      // CLV
      //
      PF_t= ntzs / nz;
      PF_i0=nxs*ny*(ntzs%nz) + nxs*ny*nz*PF_t;
      i0 = (nzys/ny)*nxs*ny + (nzys%ny)*nxs;
      out_PF_i0=out+PF_i0;
      const scs_t *in_PF_ixf = in + PF_i0;
      const pglus_t gx_PF_i0 = gx+PF_i0;
      float fmkappa = mkappa;
      for(int zy = nzys; zy < nzye; ++zy) {
        for(int t = 0; t < nt; ++t) {
          for(int x = 0; x < nxs; ++x, ++i0) {
            int _PFIT;
            if(x+_PFI < nxs){
              _PFIT=_PFI;
            }else if (t != nt-1){
              _PFIT =  nxs * ny * nz -  nxs + _PFI;
            }else{
              _PFIT = _PFI - nxs * ny * nz* (nt-1) ;
            }
#if VLENS == 16
            if ( i0+_PFIT != i0e) {
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[383]),1,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[576]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[640]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[704]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[768]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[832]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[896]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[960]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[1024]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[1088]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[1151]),0,1);
            } else{
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[383]),0,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[192]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[256]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[320]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[383]),1,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[287]),0,1);
            }
#elif VLENS == 8
#ifdef DS_TO_DOUBLE
            if ( i0+_PFIT != i0e) {
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[224]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[288]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[352]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[416]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[480]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[544]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[575]),0,1);
            } else{
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[96]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[160]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[96]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[143]),0,1);
            }
#else
            if ( i0+_PFIT != i0e) {
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[192]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[320]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[384]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[448]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[512]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[575]),0,1);
            } else{
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[128]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[191]),1,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[143]),0,1);
            }
#endif
#elif VLENS == 4
            if ( i0+_PFIT != i0e) {
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out+i0+_PFIT)->c_prefetch)[95]),1,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((in+i0+_PFIT)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[128]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[160]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[191]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[224]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[256]),0,1);
              __builtin_prefetch(&(((cl+i0+_PFIT)->c_prefetch)[287]),0,1);
            } else{
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( in_PF_ixf)->c_prefetch)[95]),0,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[0]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[32]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[64]),1,1);
              __builtin_prefetch(&(((out_PF_i0)->c_prefetch)[95]),1,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[0]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[32]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[64]),0,1);
              __builtin_prefetch(&((( gx_PF_i0)->c_prefetch)[71]),0,1);
            }
#endif

            __mult_clvs(out[i0].cv, cl[i0].cv);

#pragma loop norecurrence
            for (int c = 0; c < 3; ++c){
              for (int s = 0; s < 4; ++s){
                for (int ri = 0; ri < 2; ++ri){
                  for (int j = 0; j < VLENS; ++j){
                    ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] = ((float(*)[4][2][VLENS])((in + i0)->c))[c][s][ri][j] + ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] * fmkappa;
                  }
                }
              }
            }
          }
          i0 = i0 + nxs*ny*nz - nxs;
        }
        i0 = i0 + nxs - nxs*ny*nz*nt;
      }
#pragma omp barrier
  }
  void ddd_in_s_(scs_t* __restrict__ out, const scs_t* __restrict__ in, const int* DEO) {
	  ddd_in_s_noprl(out, in, DEO);
  }
#ifdef __cplusplus
}
#endif
