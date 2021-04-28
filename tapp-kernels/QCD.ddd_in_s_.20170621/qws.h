#ifndef QWS_H
#define QWS_H

#define VLEND 8
//#ifndef V512
//#define VLENS 8
//#else
//#define VLENS 16
//#endif
#define VLENS 16

#ifdef DS_TO_DOUBLE
//#ifndef V512
//#define VLENS 4
//#else
//#define VLENS 8
//#endif
#define VLENS 8
#define float double
#endif

#define SU3_RECONSTRUCT_D 18
#define SU3_RECONSTRUCT_S 18

#define NEO 2
#define NDIM 4
#define NCOL 3

#ifndef DS_TO_DOUBLE
#define NX 32
#else
#define NX 16
#endif
//#define NY 3
//#define NZ 6
//#define NT 4

//#define NX 16
#define NY 6
#define NZ 6
#define NT 2


#define NXH (NX/2)
#define NXS (NXH/VLENS)
#define NXD (NXH/VLEND)
#define VOLS (NXS*NY*NZ*NT)

// double precision
typedef struct{
  double v[VLEND];
} rvecd_t;

typedef union {
  double   c[3][3][2][VLEND];
  rvecd_t cv[3][3][2];
} g33d_t, *pg33d_t;

typedef union {
#if SU3_RECONSTRUCT_D == 18
  double   c[3][3][2][VLEND];
#elif SU3_RECONSTRUCT_D == 12
  double   c[2][3][2][VLEND];
#endif
} glud_t, *pglud_t;

typedef union {
  double   c[2][36][VLEND];
  rvecd_t cv[2][36];
} clvd_t, *pclvd_t;

typedef union {
  double   c[3][4][2][VLEND];
  rvecd_t cv[3][4][2];
  rvecd_t cs[12][2];
  rvecd_t ccs[24];
} scd_t;

typedef union {
  double   c[3][2][2][VLEND];
  rvecd_t cv[3][2][2];
} projscd_t;

typedef union {
  double   c[3][2][2];
} projscd1_t;

// single precision
typedef struct{
  float  v[VLENS];
} rvecs_t;

typedef union {
  float   c[3][3][2][VLENS];
  rvecs_t cv[3][3][2];
} g33s_t, *pg33s_t;

typedef union {
#if SU3_RECONSTRUCT_S == 18
  float   c[3][3][2][VLENS];
  float   c_prefetch[18*VLENS];
#elif SU3_RECONSTRUCT_S == 12
  float   c[2][3][2][VLENS];
  float   c_prefetch[12*VLENS];
#endif
} glus_t, *pglus_t;

typedef union {
  float   c[2][36][VLENS];
  float   c_prefetch[2*36*VLENS];
  rvecs_t cv[2][36];
} clvs_t, *pclvs_t;

typedef union {
  float   c[3][4][2][VLENS];
  float   c_prefetch[24*VLENS];
  rvecs_t cv[3][4][2];
  rvecs_t cs[12][2];
  rvecs_t ccs[24];
} scs_t;

typedef union {
  float   c[3][2][2][VLENS];
  rvecs_t cv[3][2][2];
} projscs_t;

typedef union {
  float   c[3][2][2];
} projscs1_t;

// omp block
typedef struct{
  int sy, sz, st;
  int ey, ez, et;
} block_map_t;

#endif
