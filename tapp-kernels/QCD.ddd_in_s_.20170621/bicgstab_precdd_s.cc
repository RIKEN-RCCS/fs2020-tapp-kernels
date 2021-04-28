#include "qws.h"
#include "qwsintrin.h"

#ifdef __cplusplus
extern "C"{
#endif

  using namespace std;
  extern int vols;
  extern void prec_ddd_s_(scs_t* out, scs_t* in, int* nsap, int* nm);

  void bicgstab_precdd_s_(scs_t* x, scs_t* b, double* tol, int* conviter, int* maxiter, int* nsap, int* nm){
    __attribute__((aligned(256))) static scs_t q [VOLS*2];
    __attribute__((aligned(256))) static scs_t r [VOLS*2];
    __attribute__((aligned(256))) static scs_t p [VOLS*2];
    __attribute__((aligned(256))) static scs_t t [VOLS*2];

    // change //
    //__attribute__((aligned(64))) static scs_t q [VOLS*2];
    //__attribute__((aligned(64))) static scs_t r [VOLS*2];
    //__attribute__((aligned(64))) static scs_t p [VOLS*2];
    //__attribute__((aligned(64))) static scs_t t [VOLS*2];
    //int ret = 0;
    //if( q==0) ret = posix_memalign((void **)&q,256,sizeof(scs_t)*VOLS*2);
    //if( r==0) ret = posix_memalign((void **)&r,256,sizeof(scs_t)*VOLS*2);
    //if( p==0) ret = posix_memalign((void **)&p,256,sizeof(scs_t)*VOLS*2);
    //if( t==0) ret = posix_memalign((void **)&t,256,sizeof(scs_t)*VOLS*2);
    
    
    int iter;

    for (iter=0; iter<(*maxiter);iter++){
      // q = Ap
      prec_ddd_s_(q, p, nsap, nm);

      // t = Ar
      prec_ddd_s_(t, r, nsap, nm);
    }//iter
  }//bicgstab


#ifdef __cplusplus
}
#endif
