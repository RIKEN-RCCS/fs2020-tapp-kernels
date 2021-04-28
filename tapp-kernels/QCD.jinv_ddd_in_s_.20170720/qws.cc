#include "qws.h"
#define _PFI1 2
#define _PFI2 1
#define _PFI3 1
#define _PFI4 2
#include "wilson_d.h"
#include "wilson_s.h"
#include "clover_d.h"
#include "clover_s.h"
#include "clover_def.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define min(a,b) (a)>(b)?(b):(a) 

#include "init.h"
#include "time.h"
#include "tools.h"

#ifdef __cplusplus
extern "C"{
#endif

  __attribute__((aligned(64))) pglud_t glud;
  __attribute__((aligned(64))) pglus_t glus;
  __attribute__((aligned(64))) pclvd_t clvd;
  __attribute__((aligned(64))) pclvs_t clvs;
  //projscd_t *xfd_buff;
  //projscd_t *xbd_buff;
  projscd1_t *xfd_send, *xfd_recv;
  projscd1_t *xbd_send, *xbd_recv;
  projscd_t *yfd_send, *yfd_recv;
  projscd_t *ybd_send, *ybd_recv;
  projscd_t *zfd_send, *zfd_recv;
  projscd_t *zbd_send, *zbd_recv;
  projscd_t *tfd_send, *tfd_recv;
  projscd_t *tbd_send, *tbd_recv;

  projscs1_t *xfs_send, *xfs_recv;
  projscs1_t *xbs_send, *xbs_recv;
  projscs_t *yfs_send, *yfs_recv;
  projscs_t *ybs_send, *ybs_recv;
  projscs_t *zfs_send, *zfs_recv;
  projscs_t *zbs_send, *zbs_recv;
  projscs_t *tfs_send, *tfs_recv;
  projscs_t *tbs_send, *tbs_recv;
  scs_t *qs;
#ifdef _MPI_
  MPI_Request sd_req[8];
  MPI_Request rd_req[8];
  MPI_Request ss_req[8];
  MPI_Request rs_req[8];
#endif

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
  int thmax,volsize;
  int iam;
#pragma omp threadprivate(iam)
  int volse;
#pragma omp threadprivate(volse)

  int *pce, *pco;//even-odd precondition notation

  block_map_t* block_map;
  int num_blocks;

  int prof_flag;
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
    if (nxd == 0) {printf("ERROR at %s : nx/2/VLEND must be >=1\n", __func__); exit(1);}
    if (nxs == 0) {printf("ERROR at %s : nx/2/VLENS must be >=1\n", __func__); exit(1);}
    vold = nxd * ny * nz * nt;
    vols = nxs * ny * nz * nt;
#ifdef _OPENMP
#pragma omp parallel 
  {
    iam = omp_get_thread_num();
#pragma omp master
   {
    thmax = omp_get_num_threads();
    volsize = vols / thmax;

   }
#pragma omp barrier
    volse = volsize * (iam+1);
  }
#else
   iam=0;
   thmax=1;
   volsize=vols
   volse=vols;

#endif

    static glus_t     glus_l     [VOLS*NDIM*NEO] ;
    static clvs_t     clvs_l     [VOLS*NEO]      ;
    fill_rand((float*)glus_l, (float*)(glus_l)+sizeof(glus_l)/sizeof(float));
    fill_rand((float*)clvs_l, (float*)(clvs_l)+sizeof(clvs_l)/sizeof(float));
    
    //printf("vold  %d \n", vold);
    //printf("vols  %d \n", vols);


    // int bx = block_size[0];
    // int by = block_size[1];
    // int bz = block_size[2];
    // int bt = block_size[3];
    // 
    // int nby = ((ny+by-1)/by);
    // int nbz = ((nz+bz-1)/bz);
    // int nbt = ((nt+bt-1)/bt);
    // int num_blocks = nby*nbz*nbt;
    // 
    // block_map = (block_map_t*)malloc(sizeof(block_map_t)* num_blocks);
    // //printf("block_num %d %d %d\n",nby,nbz, nbt);
    // //printf("block_size %d %d %d %d\n", bx, by, bz, bt);
    // i =0;
    // for (int t=0; t<nt; t+=bt){
    //   for (int z=0; z<nz; z+=bz){
    // 	for (int y=0;y<ny; y+=by){
    // 	  block_map[i].sy = y;
    // 	  block_map[i].sz = z;
    // 	  block_map[i].st = t;
    // 	  block_map[i].ey = min(y+by, ny);
    // 	  block_map[i].ez = min(z+bz, nz);
    // 	  block_map[i].et = min(t+bt, nt);
    // 	  i++;
    // 	}
    //   }
    // }

    //    g33d = (pg33d_t)malloc( sizeof(g33d_t) * vold*NDIM*NEO);
    glud = (pglud_t)malloc( sizeof(glud_t) * vold*NDIM*NEO); /////////////////
//  glus = (pglus_t)malloc( sizeof(glus_t) * vols*NDIM*NEO);
    clvd = (pclvd_t)malloc( sizeof(clvd_t) * vold*NEO);
//  clvs = (pclvs_t)malloc( sizeof(clvs_t) * vols*NEO);

    glus     = glus_l     ;
    clvs     = clvs_l     ;
    xfd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    xfd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);

//    xfs_send = (projscs1_t*)malloc( sizeof(projscs1_t) * ny *nz*nt);
//    xbs_send = (projscs1_t*)malloc( sizeof(projscs1_t) * ny *nz*nt);
//    yfs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*nz*nt);
//    ybs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*nz*nt);
//    zfs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nt);
//    zbs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nt);
//    tfs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nz);
//    tbs_send = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nz);
    posix_memalign((void **)&xfs_send,256,sizeof(projscs1_t)*ny*nz*nt);
    posix_memalign((void **)&xbs_send,256,sizeof(projscs1_t)*ny*nz*nt);
    posix_memalign((void **)&yfs_send,256,sizeof(projscs_t)*nxs*nz*nt);
    posix_memalign((void **)&ybs_send,256,sizeof(projscs_t)*nxs*nz*nt);
    posix_memalign((void **)&zfs_send,256,sizeof(projscs_t)*nxs*ny*nt);
    posix_memalign((void **)&zbs_send,256,sizeof(projscs_t)*nxs*ny*nt);
    posix_memalign((void **)&tfs_send,256,sizeof(projscs_t)*nxs*ny*nz);
    posix_memalign((void **)&tbs_send,256,sizeof(projscs_t)*nxs*ny*nz);
//    xfs_recv = (projscs1_t*)malloc( sizeof(projscs1_t) * ny *nz*nt);
//    xbs_recv = (projscs1_t*)malloc( sizeof(projscs1_t) * ny *nz*nt);
//    yfs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*nz*nt);
//    ybs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*nz*nt);
//    zfs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nt);
//    zbs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nt);
//    tfs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nz);
//    tbs_recv = (projscs_t*)malloc( sizeof(projscs_t) * nxs*ny*nz);
    posix_memalign((void **)&xfs_recv,256,sizeof(projscs1_t)*ny*nz*nt);
    posix_memalign((void **)&xbs_recv,256,sizeof(projscs1_t)*ny*nz*nt);
    posix_memalign((void **)&yfs_recv,256,sizeof(projscs_t)*nxs*nz*nt);
    posix_memalign((void **)&ybs_recv,256,sizeof(projscs_t)*nxs*nz*nt);
    posix_memalign((void **)&zfs_recv,256,sizeof(projscs_t)*nxs*ny*nt);
    posix_memalign((void **)&zbs_recv,256,sizeof(projscs_t)*nxs*ny*nt);
    posix_memalign((void **)&tfs_recv,256,sizeof(projscs_t)*nxs*ny*nz);
    posix_memalign((void **)&tbs_recv,256,sizeof(projscs_t)*nxs*ny*nz);
    posix_memalign((void **)&qs,256,sizeof(scs_t)*vols);
    for(i=0; i<vols; i++){
      for(int j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  qs[i].ccs[j].v[v]  = 0.0f;
	}
      }
    }

#ifdef _MPI_
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif

    if (size == 1){
      printf("single %d\n", rank);
      px=0;    py=0;    pz=0;    pt=0;
      domain_e = ( py + pz + pt)%2;
      domain_o = 1 - domain_e;
    } else {
      px = rank;
      //printf("MPI %d\n", px);
      pt = px / (npe[0]*npe[1]*npe[2]);
      px-= npe[0]*npe[1]*npe[2]*pt;
      pz = px / (npe[0]*npe[1]);
      px-= npe[0]*npe[1]*pz;
      py = px / npe[0];
      px-= npe[0]*py;
      pxf = ((       px+1)%npe[0]) + npe[0]*py + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
      pxb = ((npe[0]+px-1)%npe[0]) + npe[0]*py + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
      pyf = px + ((       py+1)%npe[1])*npe[0] + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
      pyb = px + ((npe[1]+py-1)%npe[1])*npe[0] + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
      pzf = px + npe[0]*py + ((       pz+1)%npe[2])*npe[0]*npe[1] + npe[0]*npe[1]*npe[2]*pt;
      pzb = px + npe[0]*py + ((npe[2]+pz-1)%npe[2])*npe[0]*npe[1] + npe[0]*npe[1]*npe[2]*pt;
      ptf = px + npe[0]*py + npe[0]*npe[1]*pz + ((       pt+1)%npe[3])*npe[0]*npe[1]*npe[2];
      ptb = px + npe[0]*py + npe[0]*npe[1]*pz + ((npe[3]+pt-1)%npe[3])*npe[0]*npe[1]*npe[2];

      //printf("MPI rank=%d px=%d py=%d pz=%d pt=%d\n", rank, px, py, pz, pt);
      //printf("MPI rank=%d pxf=%d pxb=%d\n", rank, pxf, pxb);
      //printf("MPI rank=%d pyf=%d pyb=%d\n", rank, pyf, pyb);
      //printf("MPI rank=%d pzf=%d pzb=%d\n", rank, pzf, pzb);
      //printf("MPI rank=%d ptf=%d ptb=%d\n", rank, ptf, ptb);

      domain_e = ( py + pz + pt)%2;
      domain_o = 1 - domain_e;
#ifdef _MPI_
      // tag may be changed
      MPI_Send_init(xfd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 0, MPI_COMM_WORLD, &sd_req[0]);
      MPI_Send_init(xbd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 1, MPI_COMM_WORLD, &sd_req[1]);
      MPI_Send_init(yfd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 2, MPI_COMM_WORLD, &sd_req[2]);
      MPI_Send_init(ybd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 3, MPI_COMM_WORLD, &sd_req[3]);
      MPI_Send_init(zfd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 4, MPI_COMM_WORLD, &sd_req[4]);
      MPI_Send_init(zbd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 5, MPI_COMM_WORLD, &sd_req[5]);
      MPI_Send_init(tfd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 6, MPI_COMM_WORLD, &sd_req[6]);
      MPI_Send_init(tbd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 7, MPI_COMM_WORLD, &sd_req[7]);
      MPI_Recv_init(xfd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 0, MPI_COMM_WORLD, &rd_req[0]);
      MPI_Recv_init(xbd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 1, MPI_COMM_WORLD, &rd_req[1]);
      MPI_Recv_init(yfd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 2, MPI_COMM_WORLD, &rd_req[2]);
      MPI_Recv_init(ybd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 3, MPI_COMM_WORLD, &rd_req[3]);
      MPI_Recv_init(zfd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 4, MPI_COMM_WORLD, &rd_req[4]);
      MPI_Recv_init(zbd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 5, MPI_COMM_WORLD, &rd_req[5]);
      MPI_Recv_init(tfd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 6, MPI_COMM_WORLD, &rd_req[6]);
      MPI_Recv_init(tbd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 7, MPI_COMM_WORLD, &rd_req[7]);
      
      MPI_Send_init(xfs_send, 12*ny *nz*nt, MPI_REAL, pxb, 0, MPI_COMM_WORLD, &ss_req[0]);
      MPI_Send_init(xbs_send, 12*ny *nz*nt, MPI_REAL, pxf, 1, MPI_COMM_WORLD, &ss_req[1]);
      MPI_Send_init(yfs_send, 12*nxh*nz*nt, MPI_REAL, pyb, 2, MPI_COMM_WORLD, &ss_req[2]);
      MPI_Send_init(ybs_send, 12*nxh*nz*nt, MPI_REAL, pyf, 3, MPI_COMM_WORLD, &ss_req[3]);
      MPI_Send_init(zfs_send, 12*nxh*ny*nt, MPI_REAL, pzb, 4, MPI_COMM_WORLD, &ss_req[4]);
      MPI_Send_init(zbs_send, 12*nxh*ny*nt, MPI_REAL, pzf, 5, MPI_COMM_WORLD, &ss_req[5]);
      MPI_Send_init(tfs_send, 12*nxh*ny*nz, MPI_REAL, ptb, 6, MPI_COMM_WORLD, &ss_req[6]);
      MPI_Send_init(tbs_send, 12*nxh*ny*nz, MPI_REAL, ptf, 7, MPI_COMM_WORLD, &ss_req[7]);
      MPI_Recv_init(xfs_recv, 12*ny *nz*nt, MPI_REAL, pxf, 0, MPI_COMM_WORLD, &rs_req[0]);
      MPI_Recv_init(xbs_recv, 12*ny *nz*nt, MPI_REAL, pxb, 1, MPI_COMM_WORLD, &rs_req[1]);
      MPI_Recv_init(yfs_recv, 12*nxh*nz*nt, MPI_REAL, pyf, 2, MPI_COMM_WORLD, &rs_req[2]);
      MPI_Recv_init(ybs_recv, 12*nxh*nz*nt, MPI_REAL, pyb, 3, MPI_COMM_WORLD, &rs_req[3]);
      MPI_Recv_init(zfs_recv, 12*nxh*ny*nt, MPI_REAL, pzf, 4, MPI_COMM_WORLD, &rs_req[4]);
      MPI_Recv_init(zbs_recv, 12*nxh*ny*nt, MPI_REAL, pzb, 5, MPI_COMM_WORLD, &rs_req[5]);
      MPI_Recv_init(tfs_recv, 12*nxh*ny*nz, MPI_REAL, ptf, 6, MPI_COMM_WORLD, &rs_req[6]);
      MPI_Recv_init(tbs_recv, 12*nxh*ny*nz, MPI_REAL, ptb, 7, MPI_COMM_WORLD, &rs_req[7]);
#endif
      //printf("qws_init end: MPI rank=%d px=%d py=%d pz=%d pt=%d\n", rank, px, py, pz, pt);
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


  //---------------------------------------------------------------------------------------- COMM
  void xbound(int req, int prec) {
    if (npe[req/2] != 1) {
#ifdef _MPI_
      if (prec == 8) {
	//printf("xbound %d %d\n",prec,req);
	MPI_Start(&rd_req[req]);
	MPI_Start(&sd_req[req]);
      } else {
	MPI_Start(&rs_req[req]);
	MPI_Start(&ss_req[req]);
      }
#endif
    } else {
      switch (req) {
      case 0 :
	if (prec == 8) {
	  memcpy(xfd_recv, xfd_send, sizeof(double)*12*ny*nz*nt);
	} else {
	  memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
	}
	break;
      case 1 :
	if (prec == 8) {
	  memcpy(xbd_recv, xbd_send, sizeof(double)*12*ny*nz*nt);
	} else {
	  memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
	}
	break;
      case 2 :
	if (prec == 8) {
	  memcpy(yfd_recv, yfd_send, sizeof(projscd_t)*nxd*nz*nt);
	} else {
	  memcpy(yfs_recv, yfs_send, sizeof(projscs_t)*nxs*nz*nt);
	}
	break;
      case 3 :
	if (prec == 8) {
	  memcpy(ybd_recv, ybd_send, sizeof(projscd_t)*nxd*nz*nt);
	} else {
	  memcpy(ybs_recv, ybs_send, sizeof(projscs_t)*nxs*nz*nt);
	}
	break;
      case 4 :
	if (prec == 8) {
	  memcpy(zfd_recv, zfd_send, sizeof(projscd_t)*nxd*ny*nt);
	} else {
	  memcpy(zfs_recv, zfs_send, sizeof(projscs_t)*nxs*ny*nt);
	}
	break;
      case 5 :
	if (prec == 8) {
	  memcpy(zbd_recv, zbd_send, sizeof(projscd_t)*nxd*ny*nt);
	} else {
	  memcpy(zbs_recv, zbs_send, sizeof(projscs_t)*nxs*ny*nt);
	}
	break;
      case 6 :
	if (prec == 8) {
	  memcpy(tfd_recv, tfd_send, sizeof(projscd_t)*nxd*ny*nz);
	} else {
	  memcpy(tfs_recv, tfs_send, sizeof(projscs_t)*nxs*ny*nz);
	}
	break;
      case 7 :
	if (prec == 8) {
	  memcpy(tbd_recv, tbd_send, sizeof(projscd_t)*nxd*ny*nz);
	} else {
	  memcpy(tbs_recv, tbs_send, sizeof(projscs_t)*nxs*ny*nz);
	}
	break;
      }
    }
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_wait(int req, int prec) {
#if _MPI_
    MPI_Status status;
    if (npe[req/2] != 1) {
      if (prec == 8) {
	MPI_Wait(&rd_req[req], &status);
      } else {
	MPI_Wait(&rs_req[req], &status);
      }
    }
#else
    return;
#endif
  }
  //---------------------------------------------------------------------------------------- send
  void xbound_send_waitall(int prec) {
#if _MPI_
    int i;
    MPI_Status status;
    if (prec == 8) {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&sd_req[0+2*i], &status);
	  MPI_Wait(&sd_req[1+2*i], &status);
	}
      }
    } else {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&ss_req[0+2*i], &status);
	  MPI_Wait(&ss_req[1+2*i], &status);
	}
      }
    }
#else
    return;
#endif
  }
  //---------------------------------------------------------------------------------------- recv
  void xbound_recv_waitall(int prec) {
#if _MPI_
    int i;
    MPI_Status status;
    if (prec == 8) {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&rd_req[0+2*i], &status);
	  MPI_Wait(&rd_req[1+2*i], &status);
	}
      }
    } else {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&rs_req[0+2*i], &status);
	  MPI_Wait(&rs_req[1+2*i], &status);
	}
      }
    }
#else
    return;
#endif
  }
  //---------------------------------------------------------------------------------------- for testing, not using
#ifndef RDC
  void deo_test_(int* pe, int* po, scd_t* out, scd_t* in){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;

    pglud_t gxf = &glud[vold*0 + NDIM*vold*(*pe)];
    pglud_t gxb = &glud[vold*0 + NDIM*vold*(*po)];
    pglud_t gyf = &glud[vold*1 + NDIM*vold*(*pe)];
    pglud_t gyb = &glud[vold*1 + NDIM*vold*(*po)];
    pglud_t gzf = &glud[vold*2 + NDIM*vold*(*pe)];
    pglud_t gzb = &glud[vold*2 + NDIM*vold*(*po)];
    pglud_t gtf = &glud[vold*3 + NDIM*vold*(*pe)];
    pglud_t gtb = &glud[vold*3 + NDIM*vold*(*po)];

    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick=  ( (*pe) + y     + z + t)%2;
	  //printf("ick e y z t %d %d %d %d %d\n",ick,e,y,z,t);
	  for(x=0; x<nxd; x++){
	    //printf("iiiii x=%d y=%d z=%d t=%d iyf=%d\n",x, y, z, t, iyf);
	    __zero_sc(tmp.c);
	    int i0  =      x        + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixf = (    x+1)%nxd + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixb = (nxd+x-1)%nxd + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int iyf = x + ((   y+1)%ny)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int iyb = x + ((ny+y-1)%ny)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int izf = x + nxd*y + ((   z+1)%nz)*nxd*ny + nxd*ny*nz*t;
	    int izb = x + nxd*y + ((nz+z-1)%nz)*nxd*ny + nxd*ny*nz*t;
	    int itf = x + nxd*y + nxd*ny*z + ((   t+1)%nt)*nxd*ny*nz;
	    int itb = x + nxd*y + nxd*ny*z + ((nt+t-1)%nt)*nxd*ny*nz;
	    if (ick == 0) {
	      __mult_xfwd0(tmp.c, gxf, in, i0);
	      __mult_xbwd0(tmp.c, gxb, in, i0, ixb);
	    } else {
	      __mult_xfwd1(tmp.c, gxf, in, i0, ixf);
	      __mult_xbwd1(tmp.c, gxb, in, i0);
	    }
	    __mult_yfwd(tmp.c, gyf, in, i0, iyf);
	    __mult_ybwd(tmp.c, gyb, in,     iyb);
	    __mult_zfwd(tmp.c, gzf, in, i0, izf);
	    __mult_zbwd(tmp.c, gzb, in,     izb);
	    __mult_tfwd(tmp.c, gtf, in, i0, itf);
	    __mult_tbwd(tmp.c, gtb, in,     itb);
	    __store_sc(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }

  //----------------------------------------------------------------------------------------
  void clvd_vm_(int* pe, scd_t* inout){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxd; x++){
	    int i0  = x + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    __load_sc(tmp.c, inout[i0].c);
	    //__mult_clvd( tmp.cv, clv[i0 + vold*(*pe)].cv);
	    __mult_clvd_def( tmp.cv, clvd[i0 + vold*(*pe)].cv);
	    __store_sc(inout[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  void clv_s_(int* pe, scs_t* inout){
    __attribute__((aligned(64))) scs_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxs; x++){
	    int i0  = x + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    __load_sc_s(tmp.c, inout[i0].c);
	    __mult_clvs( tmp.cv, clvs[i0 + vols*(*pe)].cv);
	    __store_sc_s(inout[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  void clvd_vm2_(int* pe, scd_t* out, scd_t* in){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxd; x++){
	    int i0  = x + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    __load_sc(tmp.c, in[i0].c);
	    __mult_clvd( tmp.cv, clvd[i0 + vold*(*pe)].cv);
	    __store_sc(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void zero_scd_field(scd_t* inout, int size){
#pragma omp parallel for
    for (int i=0;i<size;i++) {
      for (int j=0;j<24;j++) {
	for (int iv=0;iv<VLEND;iv++) {
	  inout[i].ccs[j].v[iv]=0;
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void zero_scs_field(scs_t* inout, int size){
#pragma omp parallel for
    for (int i=0;i<size;i++) {
      for (int j=0;j<24;j++) {
	for (int iv=0;iv<VLENS;iv++) {
	  inout[i].ccs[j].v[iv]=0;
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void print_scdnorm2(char* a, scd_t* in, int size){
    rvecd_t rvd0;
    double rtmp0 = 0;
#pragma omp parallel for private(rvd0) reduction(+:rtmp0)
    for(int i=0; i<size; i++){
      for(int j=0; j<24; j++){
	rvd0 = fmul_d(in[i].ccs[j],in[i].ccs[j]);
	rtmp0 += fsum_d(rvd0);
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE, &rtmp0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    if(rank==0)printf("%30s : %24.14e\n", a, rtmp0);
  }
  //----------------------------------------------------------------------------
  void print_scsnorm2(char* a, scs_t* in, int size){
    rvecs_t rvd0;
    float rtmp0 = 0;
#pragma omp parallel for private(rvd0) reduction(+:rtmp0)
    for(int i=0; i<size; i++){
      for(int j=0; j<24; j++){
	rvd0 = fmul_s(in[i].ccs[j],in[i].ccs[j]);
	rtmp0 += fsum_s(rvd0);
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE, &rtmp0,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    if(rank==0)printf("%30s : %24.14e\n", a, rtmp0);
  }
#endif
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
#ifndef RDC
  void deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void dee_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void one_minus_dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
  void one_minus_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
  void deo_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_dag_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void dee_deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void dee_deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  //---------------------------------------------------------------------------------------- double precision
  void deo_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_out_pre_vm_(pe, po, out, in);
    deo_in_vm_(pe, po, out, in);
    deo_out_pos_vm_(pe, po, out, in, 1);
  }
  void dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_out_pre_vm_(pe, po, out, in);
    dee_deo_in_vm_(pe, po, out, in);
    dee_deo_out_pos_vm_(pe, po, out, in, 1);
  }
  void one_minus_dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0){
    double fac = -kappa2;
    deo_out_pre_vm_(pe, po, out, in);
    one_minus_dee_deo_in_vm_(pe, po, out, in, in0, fac);
    dee_deo_out_pos_vm_(pe, po, out, in, fac );
  }
  void deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_dag_out_pre_vm_(pe, po, out, in);
    deo_dag_in_vm_(pe, po, out, in);
    deo_dag_out_pos_vm_(pe, po, out, in, 1);
  }
  void dee_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_dag_out_pre_vm_(pe, po, out, in);
    dee_deo_dag_in_vm_(pe, po, out, in);
    dee_deo_dag_out_pos_vm_(pe, po, out, in, 1);
  }
  void one_minus_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0){
    double fac = -kappa2;
    deo_dag_out_pre_vm_(pe, po, out, in);
    one_minus_deo_dag_in_vm_(pe, po, out, in, in0, fac);
    deo_dag_out_pos_vm_(pe, po, out, in, fac );
  }
  //----------------------------------------------------------------------------------------
  void mtilde_vm_(scd_t* out, scd_t* in){
    scd_t* tmp = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    dee_deo_vm_(pco, pce, tmp, in);
    one_minus_dee_deo_vm_(pce, pco, out, tmp, in);
    free(tmp);
  }
  void mtilde_dag_vm_(scd_t* out, scd_t* in){
    scd_t* tmp0 = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    scd_t* tmp1 = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    clvd_vm2_(pce, tmp0, in);
    dee_deo_dag_vm_(pco, pce, tmp1, tmp0);
    one_minus_deo_dag_vm_(pce, pco, out, tmp1, in);
    free(tmp0);
    free(tmp1);
  }
  void mtdagmt_vm_(scd_t* out, scd_t* in){
    scd_t* tmp = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    mtilde_vm_(tmp, in);
    mtilde_dag_vm_(out, tmp);
    free(tmp);
  }

  //---------------------------------------------------------------------------------------- single precision
  void deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
  void dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
  void one_minus_dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0, float factor);
  void deo_out_pre_s_(const int* pe, const int* po, const scs_t* out, const scs_t* in);
  void deo_out_pos_s_(const int* pe, const int* po, const scs_t* out, const scs_t* in, float factor);
  void dee_deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor);
  //----------------------------------------------------------------------------------------
  void deo_s_(int* pe, int* po, scs_t* out, scs_t* in){
    float fac = 1;
    deo_out_pre_s_(pe, po, out, in);
    deo_in_s_(pe, po, out, in);
    deo_out_pos_s_(pe, po, out, in, fac);
  }
  void dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in){
    deo_out_pre_s_(pe, po, out, in);
    dee_deo_in_s_(pe, po, out, in);
    dee_deo_out_pos_s_(pe, po, out, in, 1);
  }
  void one_minus_dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0){
    float fac = -(float)kappa2;
    deo_out_pre_s_(pe, po, out, in);
    one_minus_dee_deo_in_s_(pe, po, out, in, in0, fac);
    dee_deo_out_pos_s_(pe, po, out, in, fac );
  }
  //----------------------------------------------------------------------------------------
  void mtilde_s_(scs_t* out, scs_t* in){
    scs_t* tmp = (scs_t*)malloc( sizeof(scs_t) * nxs*ny*nz*nt);
    dee_deo_s_(pco, pce, tmp, in);
    one_minus_dee_deo_s_(pce, pco, out, tmp, in);
    free(tmp);
  }
#endif

  //---------------------------------------------------------------------------------------- for DD solver
#ifndef RDC
  void ddd_out_pre_d_(scd_t* in, int* domain);
  void ddd_out_pos_d_(scd_t* out, scd_t* in, int* domain);
  void ddd_in_d_(scd_t* out, scd_t* in, int* domain);
#endif

#ifndef RDC
  void ddd_out_pre_s_(scs_t* in, int* domain);
#endif
  void ddd_out_pre_s_noprl_(scs_t* in, int* domain);
#ifndef RDC
  void ddd_out_pre_s_no_timer_(scs_t* in, int* domain);
  void ddd_out_pre_s_noprl_no_timer_(scs_t* in, int* domain);
#endif
#ifndef RDC
  void ddd_out_pos_s_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_out_pos_s_no_timer_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_in_s_(scs_t* out, scs_t* in, int* domain);
#endif
  void ddd_out_pos_s_noprl_(scs_t* out, scs_t* in, int* domain, float factor);
#ifndef RDC
  void ddd_out_pos_s_noprl_no_timer_(scs_t* out, scs_t* in, int* domain, float factor);
#endif
  void ddd_in_s_noprl(scs_t* out, scs_t* in, int* domain);
  int domain0 = 0;
  int domain1 = 1;
  //----------------------------------------------------------------------------------------
#ifndef RDC
  void ddd_d_(scd_t* out, scd_t* in){
    ddd_out_pre_d_(in, &domain0);
    ddd_in_d_(     &out[vold*0], &in[vold*0], &domain0);
    ddd_out_pos_d_(&out[vold*0], &in[vold*0], &domain0);
    ddd_out_pre_d_(in, &domain1);
    ddd_in_d_(     &out[vold*1], &in[vold*1], &domain1);
    ddd_out_pos_d_(&out[vold*1], &in[vold*0], &domain1);
  }
  //----------------------------------------------------------------------------------------
  void ddd_eo_s_(scs_t* out, scs_t* in, int* domain){
    
    ddd_out_pre_s_no_timer_(in, domain );
    ddd_in_s_(     &out[vols*(*domain)], &in[vols*(*domain)], domain);
    ddd_out_pos_s_no_timer_(&out[vols*(*domain)], in, domain, (float)mkappa);
  }
  //----------------------------------------------------------------------------------------
  void ddd_s_(scs_t* out, scs_t* in){
    ddd_eo_s_(out, in, &domain_e);
    ddd_eo_s_(out, in, &domain_o);
  }
  void jinv_ddd_in_s_(scs_t* x, scs_t* b, int *domain, int* maxiter);
#endif
  void jinv_ddd_in_s_noprl(scs_t* x, scs_t* b, int *domain, int* maxiter);
  //----------------------------------------------------------------------------------------
  void prec_ddd_s_(scs_t* out, scs_t* in, int* nsap, int* nm){
    static scs_t *tmp0, *tmp1;
    int ret;
    ret=0;
    if( tmp0==0) ret = posix_memalign((void **)&tmp0,256,sizeof(scs_t)*vols*2);
    if( tmp1==0) ret = posix_memalign((void **)&tmp1,256,sizeof(scs_t)*vols*2);
    int i, j;

    if (*nsap > 10) { printf("prec_ddd_s_ nsap > 10\n"); exit(1); }
    if (*nm   > 10) { printf("prec_ddd_s_ nn   > 10\n"); exit(1); }

    _PREC_DDD_S_TIC_;

#pragma omp parallel 
  {
#pragma omp for private(i,j)
    for(i=0; i<vols*2; i++){
      if (i+_PFI1 < volse*2) {
        __builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[0]),0,1);
        __builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[64]),0,1);
        __builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[128]),0,1);
        __builtin_prefetch(&(((tmp0+i+_PFI1 )->c_prefetch)[0]),1,1);
        __builtin_prefetch(&(((tmp0+i+_PFI1 )->c_prefetch)[64]),1,1);
        __builtin_prefetch(&(((tmp0+i+_PFI1 )->c_prefetch)[128]),1,1);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  tmp0[i].ccs[j].v[v]  = in[i].ccs[j].v[v];
	}
      }
    }
    _OHTER_CALC_TOC_;

    for (int isap=0; isap < *nsap; isap++) {
      //
      //
      _JINV_DDD_IN_S_TIC_;
      jinv_ddd_in_s_noprl(&tmp1[vols*domain_e], &tmp0[vols*domain_e], &domain_e, nm);
      _JINV_DDD_IN_S_TOC_;

      _DDD_OUT_PRE_S_TIC_;
      ddd_out_pre_s_noprl_(   tmp1, &domain_o);
      _DDD_OUT_PRE_S_TOC_;

      _DDD_IN_S_TIC_;
      ddd_in_s_noprl( &out[vols*domain_e], &tmp1[vols*domain_e], &domain_e);
      _DDD_IN_S_TOC_;

      _OHTER_CALC_TIC_;
#pragma omp for private(i,j)
      for(i=0; i<vols; i++){
        if (i+_PFI2 < volse) {
           __builtin_prefetch(&(((  in+i+vols*domain_e+_PFI2 )->c_prefetch)[0]),0,1);
           __builtin_prefetch(&(((  in+i+vols*domain_e+_PFI2 )->c_prefetch)[64]),0,1);
           __builtin_prefetch(&(((  in+i+vols*domain_e+_PFI2 )->c_prefetch)[128]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_e+_PFI2 )->c_prefetch)[0]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_e+_PFI2 )->c_prefetch)[64]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_e+_PFI2 )->c_prefetch)[128]),0,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_e+_PFI2 )->c_prefetch)[0]),1,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_e+_PFI2 )->c_prefetch)[64]),1,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_e+_PFI2 )->c_prefetch)[128]),1,1);
        }
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    tmp0[i+vols*domain_e].ccs[j].v[v] += in[i+vols*domain_e].ccs[j].v[v] - out[i+vols*domain_e].ccs[j].v[v];
	  }
	}
      }
      _OHTER_CALC_TOC_;

      _DDD_OUT_POS_S_TIC_;
      ddd_out_pos_s_noprl_(&tmp0[vols*domain_o], tmp1, &domain_o, (float)kappa);
      _DDD_OUT_POS_S_TOC_;
      //
      //
      _JINV_DDD_IN_S_TIC_;
      jinv_ddd_in_s_noprl(&tmp1[vols*domain_o], &tmp0[vols*domain_o], &domain_o, nm);
      _JINV_DDD_IN_S_TOC_;

      _DDD_OUT_PRE_S_TIC_;
      ddd_out_pre_s_noprl_(   tmp1, &domain_e);
      _DDD_OUT_PRE_S_TOC_;

      _DDD_IN_S_TIC_;
      ddd_in_s_noprl( &out[vols*domain_o], &tmp1[vols*domain_o], &domain_o);
      _DDD_IN_S_TOC_;

      _OHTER_CALC_TIC_;
#pragma omp for private(i,j)
      for(i=0; i<vols; i++){
        if (i+_PFI2 < volse) {
           __builtin_prefetch(&(((  in+i+vols*domain_o+_PFI2 )->c_prefetch)[0]),0,1);
           __builtin_prefetch(&(((  in+i+vols*domain_o+_PFI2 )->c_prefetch)[64]),0,1);
           __builtin_prefetch(&(((  in+i+vols*domain_o+_PFI2 )->c_prefetch)[128]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_o+_PFI2 )->c_prefetch)[0]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_o+_PFI2 )->c_prefetch)[64]),0,1);
           __builtin_prefetch(&((( out+i+vols*domain_o+_PFI2 )->c_prefetch)[128]),0,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI2 )->c_prefetch)[0]),1,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI2 )->c_prefetch)[64]),1,1);
           __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI2 )->c_prefetch)[128]),1,1);
        }
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    tmp0[i+vols*domain_o].ccs[j].v[v]  += in[i+vols*domain_o].ccs[j].v[v] - out[i+vols*domain_o].ccs[j].v[v];
	  }
	}
      }
      _OHTER_CALC_TOC_;
      _DDD_OUT_POS_S_TIC_;
      ddd_out_pos_s_noprl_(&tmp0[vols*domain_e], tmp1, &domain_e, (float)kappa);
      _DDD_OUT_POS_S_TOC_;
    }

    //
    //
   _JINV_DDD_IN_S_TIC_;
    jinv_ddd_in_s_noprl(&tmp1[vols*domain_e], &tmp0[vols*domain_e], &domain_e, nm);
   _JINV_DDD_IN_S_TOC_;

    _DDD_OUT_PRE_S_TIC_;
    ddd_out_pre_s_noprl_(   tmp1, &domain_o);
    _DDD_OUT_PRE_S_TOC_;

    _DDD_IN_S_TIC_;
    ddd_in_s_noprl(        &out[vols*domain_e], &tmp1[vols*domain_e], &domain_e);
    _DDD_IN_S_TOC_;

    _OHTER_CALC_TIC_;
#pragma omp  for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI4 < volse) {
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI4 )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI4 )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI4 )->c_prefetch)[128]),1,1);
      }
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          out[i+vols*domain_o].ccs[j].v[v] = 0;
        }
      }
    }
    _OHTER_CALC_TOC_;

    _DDD_OUT_POS_S_TIC_;
    ddd_out_pos_s_noprl_(   &out[vols*domain_o], tmp1,  &domain_o, (float)mkappa);
    _DDD_OUT_POS_S_TOC_;

    _OHTER_CALC_TIC_;
#pragma omp for private(i, j)
    for(i=0; i<vols; i++){
      if (i+_PFI3 < volse) {
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[128]),0,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[128]),1,1);
      }
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          tmp0[i+vols*domain_o].ccs[j].v[v] -=  out[i+vols*domain_o].ccs[j].v[v];
        }
      }
    }
    _OHTER_CALC_TOC_;

    //
    _JINV_DDD_IN_S_TIC_;
    jinv_ddd_in_s_noprl(&tmp1[vols*domain_o], &tmp0[vols*domain_o], &domain_o, nm);
    _JINV_DDD_IN_S_TOC_;

    _DDD_OUT_PRE_S_TIC_;
    ddd_out_pre_s_noprl_(  tmp1 , &domain_e);
    _DDD_OUT_PRE_S_TOC_;
  
    _DDD_IN_S_TIC_;
    ddd_in_s_noprl( &tmp0[vols*domain_o], &tmp1[vols*domain_o], &domain_o);
    _DDD_IN_S_TOC_;

    _OHTER_CALC_TIC_;
#pragma omp for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI3 < volse) {
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&((( out+i+vols*domain_o+_PFI3 )->c_prefetch)[128]),1,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((tmp0+i+vols*domain_o+_PFI3 )->c_prefetch)[128]),0,1);
      }
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          out[i+vols*domain_o].ccs[j].v[v] += tmp0[i+vols*domain_o].ccs[j].v[v];
        }
      }
    }
    _OHTER_CALC_TOC_;
    _DDD_OUT_POS_S_TIC_;
    ddd_out_pos_s_noprl_(&out[vols*domain_e], tmp1, &domain_e, (float)mkappa);
    _DDD_OUT_POS_S_TOC_;

    _OHTER_CALC_TIC_;
  } //end parallel
    _PREC_DDD_S_TOC_;

  }

  //----------------------------------------------------------------------------------------
#ifndef RDC
  void prec_s_(scs_t* out, scs_t* in, int* nsap, int* nm){
    static scs_t *tmp0, *tmp1, *q;
    if(tmp0==0) tmp0 = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if(tmp1==0) tmp1 = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if(   q==0)    q = (scs_t*)malloc( sizeof(scs_t) * vols);
    int i, j;

    if (*nsap > 10) { printf("prec_ddd_s_ nsap > 10\n"); exit(1); }
    if (*nm   > 10) { printf("prec_ddd_s_ nn   > 10\n"); exit(1); }
    //printf("MPI rank=%d px=%d py=%d pz=%d pt=%d domain_e=%d, domain_o=%d\n", rank, px, py, pz, pt, domain_e, domain_o);
    //printf("prec_ddd_s_ start %d %d\n",*nsap,*nm);

    _PREC_S_TIC_;

#pragma omp parallel for private(i, j)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  tmp0[i].ccs[j].v[v]  = in[i].ccs[j].v[v];
	}
      }
    }

#pragma omp parallel for private(i,j)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v]  = 0;
	}
      }
    }

    for (int isap=0; isap < *nsap; isap++) {
      //
      //
      jinv_ddd_in_s_(&tmp1[vols*domain_e], &tmp0[vols*domain_e], &domain_e, nm);
      ddd_out_pre_s_no_timer_(   tmp1, &domain_o);
      ddd_in_s_(    q, &tmp1[vols*domain_e], &domain_e);
#pragma omp parallel for private(i, j)
      for(i=0; i<vols; i++){
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    tmp0[i+vols*domain_e].ccs[j].v[v]  = in[i+vols*domain_e].ccs[j].v[v] - q[i].ccs[j].v[v] + tmp0[i+vols*domain_e].ccs[j].v[v];
	  }
	}
      }
      ddd_out_pos_s_no_timer_(&tmp0[vols*domain_o], tmp1, &domain_o, (float)kappa);
      //
      //
      jinv_ddd_in_s_(&tmp1[vols*domain_o], &tmp0[vols*domain_o], &domain_o, nm);
      ddd_out_pre_s_no_timer_(   tmp1, &domain_e);
      ddd_in_s_(    q, &tmp1[vols*domain_o], &domain_o);
#pragma omp parallel for private(i,j)
      for(i=0; i<vols; i++){
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    tmp0[i+vols*domain_o].ccs[j].v[v]  = in[i+vols*domain_o].ccs[j].v[v] - q[i].ccs[j].v[v] + tmp0[i+vols*domain_o].ccs[j].v[v];
	  }
	}
      }
      ddd_out_pos_s_no_timer_(&tmp0[vols*domain_e], tmp1, &domain_e, (float)kappa);
    }


    jinv_ddd_in_s_(&out[vols*domain_e], &tmp0[vols*domain_e], &domain_e, nm);
    ddd_out_pre_s_no_timer_(   out, &domain_o);
#pragma omp parallel for private(i, j)
    for(i=0; i<vols; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
	  q[i].ccs[j].v[v]  = tmp0[i+vols*domain_o].ccs[j].v[v];
	}
      }
    }
    ddd_out_pos_s_no_timer_(   q, tmp1,          &domain_o, (float)kappa);
    jinv_ddd_in_s_(&out[vols*domain_o], q, &domain_o, nm);

    //free(tmp0);
    //free(tmp1);
    //free(q);

    _PREC_S_TOC_;
  }
#endif


#ifdef __cplusplus
}
#endif
