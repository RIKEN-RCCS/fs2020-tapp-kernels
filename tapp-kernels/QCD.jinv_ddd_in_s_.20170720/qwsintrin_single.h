#ifndef QWS_INTRINSICS_SINGLE
#define QWS_INTRINSICS_SINGLE

//=0;
inline rvecs_t fzero_s(void) __attribute__((always_inline));
inline rvecs_t fzero_s(void){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = 0;}
  return tmp;
}
//=a;
inline rvecs_t fcopy_s(const rvecs_t &a) __attribute__((always_inline));
inline rvecs_t fcopy_s(const rvecs_t &a){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i];}
  return tmp;
}
//=-a;
inline rvecs_t mfcopy_s(const rvecs_t &a) __attribute__((always_inline));
inline rvecs_t mfcopy_s(const rvecs_t &a){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] =-a.v[i];}
  return tmp;
}

//a[v]=a;
inline rvecs_t fload1_s(const float &a) __attribute__((always_inline));
inline rvecs_t fload1_s(const float &a){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a;}
  return tmp;
}
// a+b;
inline rvecs_t fadd_s(const rvecs_t &a, const rvecs_t &b) __attribute__((always_inline));
inline rvecs_t fadd_s(const rvecs_t &a, const rvecs_t &b){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i]+b.v[i];}
  return tmp;
}
// a-b;
inline rvecs_t fsub_s(const rvecs_t &a, const rvecs_t &b) __attribute__((always_inline));
inline rvecs_t fsub_s(const rvecs_t &a, const rvecs_t &b){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i]-b.v[i];}
  return tmp;
}

// a*b;
inline rvecs_t fmul_s(const rvecs_t &a, const rvecs_t &b) __attribute__((always_inline));
inline rvecs_t fmul_s(const rvecs_t &a, const rvecs_t &b){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i]*b.v[i];}
  return tmp;
}

// a*b+c;
inline rvecs_t fmadd_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c) __attribute__((always_inline));
inline rvecs_t fmadd_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i]*b.v[i]+c.v[i];}
  return tmp;
}

//-a*b+c;
inline rvecs_t fnmadd_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c) __attribute__((always_inline));
inline rvecs_t fnmadd_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] =-a.v[i]*b.v[i]+c.v[i];}
  return tmp;
}

// a*b-c;
inline rvecs_t fmsub_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c) __attribute__((always_inline));
inline rvecs_t fmsub_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i]*b.v[i]-c.v[i];}
  return tmp;
}

//-a*b-c;
inline rvecs_t fnmsub_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c) __attribute__((always_inline));
inline rvecs_t fnmsub_s(const rvecs_t &a, const rvecs_t &b,const rvecs_t &c){
  int i;
  rvecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] =-a.v[i]*b.v[i]-c.v[i];}
  return tmp;
}

//=sum(a);
inline float fsum_s(const rvecs_t &a) __attribute__((always_inline));
inline float fsum_s(const rvecs_t &a){
  int i;
  float tmp=0;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp += a.v[i];}
  return tmp;
}
#endif
