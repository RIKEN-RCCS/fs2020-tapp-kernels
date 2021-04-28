#include "qws.h"
#include "wilson_s.h"
#include "clover_s.h"
#include "profiler.h"
#include "init.h"
#include "omp.h"

#define min(a,b) (a)>(b)?(b):(a) 

#ifdef __cplusplus
extern "C"{
#endif

  __attribute__((aligned(64))) pglus_t glus;
  __attribute__((aligned(64))) pclvs_t clvs;

  projscs1_t *xfs_send, *xfs_recv;
  projscs1_t *xbs_send, *xbs_recv;
  projscs_t *yfs_send, *yfs_recv;
  projscs_t *ybs_send, *ybs_recv;
  projscs_t *zfs_send, *zfs_recv;
  projscs_t *zbs_send, *zbs_recv;
  projscs_t *tfs_send, *tfs_recv;
  projscs_t *tbs_send, *tbs_recv;

  double kappa, kappa2, mkappa;

  //  nh=nx/2, nd=nx/2/8, nxs=nx/2/16
  int nx, ny, nz, nt, nxh, nxd, nxs;
  int vold;// = nxd * ny * nz * nt;
  int vols;// = nxs * ny * nz * nt;
  int npe[4];// number of proccesses
  double fbc[4][2];// fermion boundary condiion [mu][0=fwd, 1=bwd]
  int rank, size, px, py, pz, pt;// process coordinate
  int pxf, pyf, pzf, ptf;
  int pxb, pyb, pzb, ptb;
  int domain_e, domain_o;
  int thmax;
//#pragma omp threadprivate(iam)

  int *pce, *pco;//even-odd precondition notation

  block_map_t* block_map;
  int num_blocks;

  //---------------------------------------------------------------------------------------- init
  void qws_init_(int* lx,  int* ly, int* lz, int* lt, 
		 int* npe_f, int* fbc_f, int* pce_f, int* pco_f, int* block_size){
    int i;
    pce = pce_f;
    pco = pco_f;

    for (i=0;i<4;i++){
      npe[i]=npe_f[i];
    }

    nx=*lx;  ny=*ly;  nz=*lz;  nt=*lt;
    nxh=nx/2;
    nxd=nx/2/VLEND;
    nxs=nx/2/VLENS;
    vold = nxd * ny * nz * nt;
    vols = nxs * ny * nz * nt;
#pragma omp parallel 
    {
      int iam = omp_get_thread_num();
      if(iam == 0) {
        thmax = omp_get_num_threads();
      }
#pragma omp barrier
    }

    static glus_t     glus_l     [VOLS*NDIM*NEO] ;
    static clvs_t     clvs_l     [VOLS*NEO]      ;
    fill_rand((float*)glus_l, (float*)(glus_l)+sizeof(glus_l)/sizeof(float));
    fill_rand((float*)clvs_l, (float*)(clvs_l)+sizeof(clvs_l)/sizeof(float));

    // static projscs1_t xfs_send_l [NY*NZ*NT]      ;
    // static projscs1_t xbs_send_l [NY*NZ*NT]      ;
    // static projscs_t  yfs_send_l [NXS*NZ*NT]     ;
    // static projscs_t  ybs_send_l [NXS*NZ*NT]     ;
    // static projscs_t  zfs_send_l [NXS*NY*NT]     ;
    // static projscs_t  zbs_send_l [NXS*NY*NT]     ;
    // static projscs_t  tfs_send_l [NXS*NY*NZ]     ;
    // static projscs_t  tbs_send_l [NXS*NY*NZ]     ;
    // static projscs1_t xfs_recv_l [NY*NZ*NT]      ;
    // static projscs1_t xbs_recv_l [NY*NZ*NT]      ;
    // static projscs_t  yfs_recv_l [NXS*NZ*NT]     ;
    // static projscs_t  ybs_recv_l [NXS*NZ*NT]     ;
    // static projscs_t  zfs_recv_l [NXS*NY*NT]     ;
    // static projscs_t  zbs_recv_l [NXS*NY*NT]     ;
    // static projscs_t  tfs_recv_l [NXS*NY*NZ]     ;
    // static projscs_t  tbs_recv_l [NXS*NY*NZ]     ;

    glus     = glus_l     ;
    clvs     = clvs_l     ;
    // xfs_send = xfs_send_l ;
    // xbs_send = xbs_send_l ;
    // yfs_send = yfs_send_l ;
    // ybs_send = ybs_send_l ;
    // zfs_send = zfs_send_l ;
    // zbs_send = zbs_send_l ;
    // tfs_send = tfs_send_l ;
    // tbs_send = tbs_send_l ;
    // xfs_recv = xfs_recv_l ;
    // xbs_recv = xbs_recv_l ;
    // yfs_recv = yfs_recv_l ;
    // ybs_recv = ybs_recv_l ;
    // zfs_recv = zfs_recv_l ;
    // zbs_recv = zbs_recv_l ;
    // tfs_recv = tfs_recv_l ;
    // tbs_recv = tbs_recv_l ;

    rank = 0;
    size = 1;

    if (size == 1){
      px=0;    py=0;    pz=0;    pt=0;
      domain_e = ( py + pz + pt)%2;
      domain_o = 1 - domain_e;
    }
    if (pt == npe[3]-1 && pt != 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)(fbc_f[3]*2);
      fbc[3][1]=(double)2;
    } else if (pt != npe[3]-1 && pt == 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)2;
      fbc[3][1]=(double)(fbc_f[3]*2);
    } else if (pt == npe[3]-1 && pt == 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)(fbc_f[3]*2);
      fbc[3][1]=(double)(fbc_f[3]*2);
    } else {
      fbc[3][0]=(double)2;
      fbc[3][1]=(double)2;
    }
    //exit(0);
  }

  //---------------------------------------------------------------------------------------- for DD solver
  void ddd_out_pre_s_(scs_t* in, int* domain);
  void ddd_out_pos_s_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_in_s_(scs_t* out, scs_t* in, int* domain, int iam);
  int domain0 = 0;
  int domain1 = 1;
  void jinv_ddd_in_s_(scs_t* x, scs_t* b, int *domain, int* maxiter);
  //----------------------------------------------------------------------------------------

  extern float* result_arr;
  extern int result_arr_len;
  void prec_ddd_s_(scs_t* out, scs_t* in, int* nsap, int* nm){
    __attribute__((aligned(256))) static scs_t tmp0 [VOLS*2] ;
    __attribute__((aligned(256))) static scs_t tmp1 [VOLS*2] ;
    __attribute__((aligned(256))) static scs_t q [VOLS*2] ;

    // changeed //
    //static scs_t *tmp0, *tmp1, *q;
    //int ret = 0;
    //if( tmp0==0) ret = posix_memalign((void **)&tmp0,256,sizeof(scs_t)*VOLS*2);
    //if( tmp1==0) ret = posix_memalign((void **)&tmp1,256,sizeof(scs_t)*VOLS*2);
    //if( q == 0) ret = posix_memalign((void **)&q, 256, sizeof(scs_t)*VOLS*2);


    static bool first_call = true;

    if(first_call) {
      fill_rand((float*)tmp0, (float*)(tmp0)+sizeof(tmp0)/sizeof(float));
      fill_rand((float*)tmp1, (float*)(tmp1)+sizeof(tmp1)/sizeof(float));
      fill_rand((float*)q,    (float*)(q)   +sizeof(q)   /sizeof(float));
      result_arr = (float*)q;
      result_arr_len = sizeof(q)/sizeof(float);
      first_call = false;
    }

#pragma omp parallel
    {
	int iam2 = omp_get_thread_num();
      for (int isap=0; isap < *nsap; isap++) {
        ddd_in_s_(    q, &tmp1[vols*domain_e], &domain_e, iam2);
#pragma omp barrier
        ddd_in_s_(    q, &tmp1[vols*domain_o], &domain_o, iam2);
#pragma omp barrier
      }
      ddd_in_s_(        &out[vols*domain_e], &tmp1[vols*domain_e], &domain_e, iam2);
#pragma omp barrier
      ddd_in_s_(       q, &tmp1[vols*domain_o], &domain_o, iam2);
    }
  }

#ifdef __cplusplus
}
#endif
