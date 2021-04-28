#ifndef QWS_INTRINSICS_DOUBLE
#define QWS_INTRINSICS_DOUBLE

//=0;
inline rvecd_t fzero_d(void) __attribute__((always_inline));
inline rvecd_t fzero_d(void){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = 0;}
  return tmp;
}
//=a;
inline rvecd_t fcopy_d(const rvecd_t &a) __attribute__((always_inline));
inline rvecd_t fcopy_d(const rvecd_t &a){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i];}
  return tmp;
}
//=-a;
inline rvecd_t mfcopy_d(const rvecd_t &a) __attribute__((always_inline));
inline rvecd_t mfcopy_d(const rvecd_t &a){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] =-a.v[i];}
  return tmp;
}

//a[v]=a;
inline rvecd_t fload1_d(const double &a) __attribute__((always_inline));
inline rvecd_t fload1_d(const double &a){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a;}
  return tmp;
}
// a+b;
inline rvecd_t fadd_d(const rvecd_t &a, const rvecd_t &b) __attribute__((always_inline));
inline rvecd_t fadd_d(const rvecd_t &a, const rvecd_t &b){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i]+b.v[i];}
  return tmp;
}
// a-b;
inline rvecd_t fsub_d(const rvecd_t &a, const rvecd_t &b) __attribute__((always_inline));
inline rvecd_t fsub_d(const rvecd_t &a, const rvecd_t &b){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i]-b.v[i];}
  return tmp;
}

// a*b;
inline rvecd_t fmul_d(const rvecd_t &a, const rvecd_t &b) __attribute__((always_inline));
inline rvecd_t fmul_d(const rvecd_t &a, const rvecd_t &b){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i]*b.v[i];}
  return tmp;
}

// a*b+c;
inline rvecd_t fmadd_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c) __attribute__((always_inline));
inline rvecd_t fmadd_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i]*b.v[i]+c.v[i];}
  return tmp;
}

//-a*b+c;
inline rvecd_t fnmadd_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c) __attribute__((always_inline));
inline rvecd_t fnmadd_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] =-a.v[i]*b.v[i]+c.v[i];}
  return tmp;
}

// a*b-c;
inline rvecd_t fmsub_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c) __attribute__((always_inline));
inline rvecd_t fmsub_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] = a.v[i]*b.v[i]-c.v[i];}
  return tmp;
}

//-a*b-c;
inline rvecd_t fnmsub_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c) __attribute__((always_inline));
inline rvecd_t fnmsub_d(const rvecd_t &a, const rvecd_t &b,const rvecd_t &c){
  int i;
  rvecd_t tmp;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp.v[i] =-a.v[i]*b.v[i]-c.v[i];}
  return tmp;
}

//=sum(a);
inline double fsum_d(const rvecd_t &a) __attribute__((always_inline));
inline double fsum_d(const rvecd_t &a){
  int i;
  double tmp=0;
#pragma loop nounroll
  for (i=0;i<VLEND;i++){tmp += a.v[i];}
  return tmp;
}
#endif
