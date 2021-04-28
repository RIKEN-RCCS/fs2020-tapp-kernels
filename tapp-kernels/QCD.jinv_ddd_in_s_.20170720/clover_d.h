#ifndef CLOVER_D_H
#define CLOVER_D_H

#include "qws.h"
#include "qwsintrin.h"
static inline void __mult_clvd(rvecd_t sc[3][4][2], rvecd_t a[2][36]) {
  int c, ri, i;
  rvecd_t x[2][6][2];
  rvecd_t y[2][2];

  for (c=0;c<3;c++){
    for (ri=0;ri<2;ri++){
      x[0][0+c][ri] = fadd_d(sc[c][0][ri], sc[c][2][ri]);
      x[0][3+c][ri] = fadd_d(sc[c][1][ri], sc[c][3][ri]);
      x[1][0+c][ri] = fsub_d(sc[c][0][ri], sc[c][2][ri]);
      x[1][3+c][ri] = fsub_d(sc[c][1][ri], sc[c][3][ri]);
    }
  }

  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][0], x[i][0][0]);
    y[i][1] = fmul_d(a[i][0], x[i][0][1]);
    y[i][0] = fnmadd_d(a[i][15], x[i][5][1], fmadd_d(a[i][14], x[i][5][0], 
              fnmadd_d(a[i][13], x[i][4][1], fmadd_d(a[i][12], x[i][4][0], 
              fnmadd_d(a[i][11], x[i][3][1], fmadd_d(a[i][10], x[i][3][0], 
              fnmadd_d(a[i][ 9], x[i][2][1], fmadd_d(a[i][ 8], x[i][2][0], 
              fnmadd_d(a[i][ 7], x[i][1][1], fmadd_d(a[i][ 6], x[i][1][0], y[i][0]))))))))));
    y[i][1] =  fmadd_d(a[i][15], x[i][5][0], fmadd_d(a[i][14], x[i][5][1], 
               fmadd_d(a[i][13], x[i][4][0], fmadd_d(a[i][12], x[i][4][1], 
               fmadd_d(a[i][11], x[i][3][0], fmadd_d(a[i][10], x[i][3][1], 
               fmadd_d(a[i][ 9], x[i][2][0], fmadd_d(a[i][ 8], x[i][2][1], 
               fmadd_d(a[i][ 7], x[i][1][0], fmadd_d(a[i][ 6], x[i][1][1], y[i][1]))))))))));
  }

  sc[0][0][0] = fadd_d(y[0][0], y[1][0]);
  sc[0][0][1] = fadd_d(y[0][1], y[1][1]);
  sc[0][2][0] = fsub_d(y[0][0], y[1][0]);
  sc[0][2][1] = fsub_d(y[0][1], y[1][1]);

  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][1], x[i][1][0]);
    y[i][1] = fmul_d(a[i][1], x[i][1][1]);
    y[i][0] = fnmadd_d(a[i][23], x[i][5][1], fmadd_d(a[i][22], x[i][5][0], 
              fnmadd_d(a[i][21], x[i][4][1], fmadd_d(a[i][20], x[i][4][0],
              fnmadd_d(a[i][19], x[i][3][1], fmadd_d(a[i][18], x[i][3][0], 
              fnmadd_d(a[i][17], x[i][2][1], fmadd_d(a[i][16], x[i][2][0],
               fmadd_d(a[i][ 7], x[i][0][1], fmadd_d(a[i][ 6], x[i][0][0], y[i][0]))))))))));
    y[i][1] =  fmadd_d(a[i][23], x[i][5][0], fmadd_d(a[i][22], x[i][5][1], 
               fmadd_d(a[i][21], x[i][4][0], fmadd_d(a[i][20], x[i][4][1],
	       fmadd_d(a[i][19], x[i][3][0], fmadd_d(a[i][18], x[i][3][1], 
               fmadd_d(a[i][17], x[i][2][0], fmadd_d(a[i][16], x[i][2][1],
	      fnmadd_d(a[i][ 7], x[i][0][0], fmadd_d(a[i][ 6], x[i][0][1], y[i][1]))))))))));
  }

  sc[1][0][0] = fadd_d(y[0][0], y[1][0]);
  sc[1][0][1] = fadd_d(y[0][1], y[1][1]);
  sc[1][2][0] = fsub_d(y[0][0], y[1][0]);
  sc[1][2][1] = fsub_d(y[0][1], y[1][1]);

  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][2], x[i][2][0]);
    y[i][1] = fmul_d(a[i][2], x[i][2][1]);
    y[i][0] = fnmadd_d(a[i][29], x[i][5][1], fmadd_d(a[i][28], x[i][5][0], 
              fnmadd_d(a[i][27], x[i][4][1], fmadd_d(a[i][26], x[i][4][0], 
              fnmadd_d(a[i][25], x[i][3][1], fmadd_d(a[i][24], x[i][3][0], 
               fmadd_d(a[i][17], x[i][1][1], fmadd_d(a[i][16], x[i][1][0], 
               fmadd_d(a[i][ 9], x[i][0][1], fmadd_d(a[i][ 8], x[i][0][0], y[i][0]))))))))));
    y[i][1] =  fmadd_d(a[i][29], x[i][5][0], fmadd_d(a[i][28], x[i][5][1], 
               fmadd_d(a[i][27], x[i][4][0], fmadd_d(a[i][26], x[i][4][1], 
               fmadd_d(a[i][25], x[i][3][0], fmadd_d(a[i][24], x[i][3][1], 
              fnmadd_d(a[i][17], x[i][1][0], fmadd_d(a[i][16], x[i][1][1], 
              fnmadd_d(a[i][ 9], x[i][0][0], fmadd_d(a[i][ 8], x[i][0][1], y[i][1]))))))))));
  }

  sc[2][0][0] = fadd_d(y[0][0], y[1][0]);
  sc[2][0][1] = fadd_d(y[0][1], y[1][1]);
  sc[2][2][0] = fsub_d(y[0][0], y[1][0]);
  sc[2][2][1] = fsub_d(y[0][1], y[1][1]);

  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][3], x[i][3][0]);
    y[i][1] = fmul_d(a[i][3], x[i][3][1]);
    y[i][0] = fnmadd_d(a[i][33], x[i][5][1], fmadd_d(a[i][32], x[i][5][0], 
              fnmadd_d(a[i][31], x[i][4][1], fmadd_d(a[i][30], x[i][4][0],
               fmadd_d(a[i][25], x[i][2][1], fmadd_d(a[i][24], x[i][2][0], 
               fmadd_d(a[i][19], x[i][1][1], fmadd_d(a[i][18], x[i][1][0], 
               fmadd_d(a[i][11], x[i][0][1], fmadd_d(a[i][10], x[i][0][0], y[i][0]))))))))));
    y[i][1] =  fmadd_d(a[i][33], x[i][5][0], fmadd_d(a[i][32], x[i][5][1], 
               fmadd_d(a[i][31], x[i][4][0], fmadd_d(a[i][30], x[i][4][1], 
              fnmadd_d(a[i][25], x[i][2][0], fmadd_d(a[i][24], x[i][2][1], 
              fnmadd_d(a[i][19], x[i][1][0], fmadd_d(a[i][18], x[i][1][1],
              fnmadd_d(a[i][11], x[i][0][0], fmadd_d(a[i][10], x[i][0][1], y[i][1]))))))))));
  }

  sc[0][1][0] = fadd_d(y[0][0], y[1][0]);
  sc[0][1][1] = fadd_d(y[0][1], y[1][1]);
  sc[0][3][0] = fsub_d(y[0][0], y[1][0]);
  sc[0][3][1] = fsub_d(y[0][1], y[1][1]);


  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][4], x[i][4][0]);
    y[i][1] = fmul_d(a[i][4], x[i][4][1]);
    y[i][0] = fnmadd_d(a[i][35], x[i][5][1], fmadd_d(a[i][34], x[i][5][0], 
               fmadd_d(a[i][31], x[i][3][1], fmadd_d(a[i][30], x[i][3][0],
               fmadd_d(a[i][27], x[i][2][1], fmadd_d(a[i][26], x[i][2][0], 
               fmadd_d(a[i][21], x[i][1][1], fmadd_d(a[i][20], x[i][1][0], 
               fmadd_d(a[i][13], x[i][0][1], fmadd_d(a[i][12], x[i][0][0], y[i][0]))))))))));
    y[i][1] =  fmadd_d(a[i][35], x[i][5][0], fmadd_d(a[i][34], x[i][5][1], 
              fnmadd_d(a[i][31], x[i][3][0], fmadd_d(a[i][30], x[i][3][1], 
              fnmadd_d(a[i][27], x[i][2][0], fmadd_d(a[i][26], x[i][2][1], 
              fnmadd_d(a[i][21], x[i][1][0], fmadd_d(a[i][20], x[i][1][1], 
              fnmadd_d(a[i][13], x[i][0][0], fmadd_d(a[i][12], x[i][0][1], y[i][1]))))))))));
  }

  sc[1][1][0] = fadd_d(y[0][0], y[1][0]);
  sc[1][1][1] = fadd_d(y[0][1], y[1][1]);
  sc[1][3][0] = fsub_d(y[0][0], y[1][0]);
  sc[1][3][1] = fsub_d(y[0][1], y[1][1]);

  for (i=0;i<2;i++){
    y[i][0] = fmul_d(a[i][5], x[i][5][0]);
    y[i][1] = fmul_d(a[i][5], x[i][5][1]);
    y[i][0] =  fmadd_d(a[i][35], x[i][4][1], fmadd_d(a[i][34], x[i][4][0], 
               fmadd_d(a[i][33], x[i][3][1], fmadd_d(a[i][32], x[i][3][0],
	       fmadd_d(a[i][29], x[i][2][1], fmadd_d(a[i][28], x[i][2][0], 
               fmadd_d(a[i][23], x[i][1][1], fmadd_d(a[i][22], x[i][1][0],
	       fmadd_d(a[i][15], x[i][0][1], fmadd_d(a[i][14], x[i][0][0], y[i][0]))))))))));
    y[i][1] = fnmadd_d(a[i][35], x[i][4][0], fmadd_d(a[i][34], x[i][4][1], 
              fnmadd_d(a[i][33], x[i][3][0], fmadd_d(a[i][32], x[i][3][1],
	      fnmadd_d(a[i][29], x[i][2][0], fmadd_d(a[i][28], x[i][2][1], 
              fnmadd_d(a[i][23], x[i][1][0], fmadd_d(a[i][22], x[i][1][1],
	      fnmadd_d(a[i][15], x[i][0][0], fmadd_d(a[i][14], x[i][0][1], y[i][1]))))))))));
  }

  sc[2][1][0] = fadd_d(y[0][0], y[1][0]);
  sc[2][1][1] = fadd_d(y[0][1], y[1][1]);
  sc[2][3][0] = fsub_d(y[0][0], y[1][0]);
  sc[2][3][1] = fsub_d(y[0][1], y[1][1]);

}
#endif
