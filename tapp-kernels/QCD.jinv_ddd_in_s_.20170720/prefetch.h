#ifdef PREFETCH

#define __prefetch_su3(a,i){\
  __builtin_prefetch((void*)&((*(a+i)).c_prefetch[0  ]), 0, 1);\
  __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64 ]), 0, 1);\
  __builtin_prefetch((void*)&((*(a+i)).c_prefetch[128]), 0, 1);\
  }
#define __prefetch_inp(a,i){\
  __builtin_prefetch(&(*(a+i)).c_prefetch[0  ], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[64 ], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[128], 0, 1);\
  }
#define __prefetch_out(a,i){\
  __builtin_prefetch(&(*(a+i)).c_prefetch[0  ], 1, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[64 ], 1, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[128], 1, 1);\
  }
#define __prefetch_clv(a,i){\
  __builtin_prefetch(&(*(a+i)).c_prefetch[0  ], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[64 ], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[128], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[192], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[256], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[320], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[384], 0, 1);\
  __builtin_prefetch(&(*(a+i)).c_prefetch[512], 0, 1);\
  }
#else

#define __prefetch_su3(a, i)
#define __prefetch_inp(a, i)
#define __prefetch_out(a, i)
#define __prefetch_clv(a, i)

#endif

