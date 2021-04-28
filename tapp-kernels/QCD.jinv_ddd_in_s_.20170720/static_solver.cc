#include "qws.h"
#define _PFI 1
#include "qwsintrin.h"
#include <math.h>


#ifdef __cplusplus
extern "C"{
#endif

  extern int vols;
  extern scs_t * qs;
  extern int iam,volse,volsize;
#pragma omp threadprivate(volse)
  extern void ddd_in_s_noprl(scs_t* out, scs_t* in, int *DEO);
#include "time.h"

  void jinv_ddd_in_s_noprl(scs_t* x, scs_t* b, int *DEO, int* maxiter){
    int i, j, iter;

    ddd_in_s_noprl(qs, b, DEO);
#pragma omp  for private(i, j) 
    for(i=0; i<vols; i++){
      if (i+_PFI < volse) {
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[128]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[128]),0,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[128]),1,1);

      }else{
       __builtin_prefetch(&(((qs+i+_PFI-volsize )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&(((qs+i+_PFI-volsize )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&(((qs+i+_PFI-volsize )->c_prefetch)[128]),1,1);
       __builtin_prefetch(&(((x+i+_PFI-volsize )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((x+i+_PFI-volsize )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((x+i+_PFI-volsize )->c_prefetch)[128]),0,1);
      }
      for(j=0; j<24; j++){
        for(int v=0; v < VLENS; v++) {
          x[i].ccs[j].v[v] = 2*b[i].ccs[j].v[v] - qs[i].ccs[j].v[v];
        }
      }
    }
    for (iter=1; iter<(*maxiter);iter++){
      // q = Ax
      ddd_in_s_noprl(qs, x, DEO);
#pragma omp  for private(i, j) 
      for(i=0; i<vols; i++){
      if (i+_PFI < volse) {
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((b+i+_PFI )->c_prefetch)[128]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[0]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[64]),0,1);
       __builtin_prefetch(&(((qs+i+_PFI )->c_prefetch)[128]),0,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[0]),1,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[64]),1,1);
       __builtin_prefetch(&(((x+i+_PFI )->c_prefetch)[128]),1,1);

      }
	for(j=0; j<24; j++){
	  for(int v=0; v < VLENS; v++) {
	    x[i].ccs[j].v[v] += b[i].ccs[j].v[v] - qs[i].ccs[j].v[v];
	  }
	}
      }
    }//iter
    // free(q);
  }//mr
  void jinv_ddd_in_s_(scs_t* x, scs_t* b, int *DEO, int* maxiter){
#pragma omp parallel
  {
    jinv_ddd_in_s_noprl(x, b, DEO, maxiter);
  }
  }
#ifdef __cplusplus
}
#endif
