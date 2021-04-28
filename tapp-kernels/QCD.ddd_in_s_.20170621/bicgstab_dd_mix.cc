#include "qws.h"
#include "qwsintrin.h"
#include "profiler.h"

#ifdef __cplusplus
extern "C"{
#endif

  using namespace std;
  extern int rank, vold, vols;
  extern void bicgstab_precdd_s_(scs_t* x, scs_t* b, double* tol, int* iter, int* maxiter, int* nsap, int* nm);
  //----------------------------------------------------------------------------------------
  void approx_inv_dirac_op(scd_t* x, scd_t* b, double *tol, int* iter, int* maxiter, int* nsap, int* nm){
    __attribute__((aligned(256))) static scs_t b_s[VOLS*2];
    __attribute__((aligned(256))) static scs_t t_s[VOLS*2];

    // change //
    //__attribute__((aligned(64))) static scs_t *b_s;
    //static scs_t* t_s;
    //static bool is_first = true;
    //if(is_first){
    //  b_s = (scs_t*)malloc(sizeof(scs_t) * VOLS*2);
    //  posix_memalign((void **)&t_s, 256, sizeof(scs_t)*VOLS*2);
    //}

    bicgstab_precdd_s_(t_s, b_s, tol, iter, maxiter, nsap, nm);
  }
  //----------------------------------------------------------------------------------------
  void bicgstab_dd_mix_(scd_t* x, scd_t* b, double *tol, int* conviter, int* maxiter, double *tol_s, int* maxiter_s, int* nsap, int* nm){
    int iter;
    static int iter_s_l;
    int *iter_s = &iter_s_l;

#ifndef RDC
    PROF_START("ddd_in_s_");
#endif
    for (iter=1; iter<(*maxiter);iter++){
#ifdef RDC
      if(iter==2) PROF_START("ddd_in_s_"); // iter==1: warm up
#endif
      // u = Mp
      // q = Au
      approx_inv_dirac_op(0, 0, tol_s, iter_s, maxiter_s, nsap, nm);

      // u = Mr
      // t = Au
      approx_inv_dirac_op(0, 0, tol_s, iter_s, maxiter_s, nsap, nm);
    }//iter
    PROF_STOP("ddd_in_s_");
  }//bicgstab


#ifdef __cplusplus
}
#endif
