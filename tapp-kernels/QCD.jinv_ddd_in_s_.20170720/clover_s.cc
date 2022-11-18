#include "qws.h"
#include "qwsintrin.h"
int true1 = 1, true2 = 1, true3 = 1;
void __mult_clvs(rvecs_t sc[3][4][2], rvecs_t a[2][36]) {
#pragma procedure noeval
  rvecs_t x[2][6][2];

  float z0_000[VLENS];
  float z0_001[VLENS];
  float z0_010[VLENS];
  float z0_011[VLENS];
  float z0_100[VLENS];
  float z0_101[VLENS];
  float z0_110[VLENS];
  float z0_111[VLENS];

  float z0_200[VLENS];
  float z0_201[VLENS];
  float z0_210[VLENS];
  float z0_211[VLENS];
  float z0_300[VLENS];
  float z0_301[VLENS];
  float z0_310[VLENS];
  float z0_311[VLENS];

  float z0_400[VLENS];
  float z0_401[VLENS];
  float z0_410[VLENS];
  float z0_411[VLENS];
  float z0_500[VLENS];
  float z0_501[VLENS];
  float z0_510[VLENS];
  float z0_511[VLENS];

  float z1_000[VLENS];
  float z1_001[VLENS];
  float z1_010[VLENS];
  float z1_011[VLENS];
  float z1_100[VLENS];
  float z1_101[VLENS];
  float z1_110[VLENS];
  float z1_111[VLENS];

  float z1_200[VLENS];
  float z1_201[VLENS];
  float z1_210[VLENS];
  float z1_211[VLENS];
  float z1_300[VLENS];
  float z1_301[VLENS];
  float z1_310[VLENS];
  float z1_311[VLENS];

  float z1_400[VLENS];
  float z1_401[VLENS];
  float z1_410[VLENS];
  float z1_411[VLENS];
  float z1_500[VLENS];
  float z1_501[VLENS];
  float z1_510[VLENS];
  float z1_511[VLENS];

  for (int c=0;c<3;c++)
    for (int ri=0;ri<2;ri++)
      for(int v=0; v < VLENS; v++) {
        x[0][0+c][ri].v[v] = sc[c][0][ri].v[v] + sc[c][2][ri].v[v];
        x[0][3+c][ri].v[v] = sc[c][1][ri].v[v] + sc[c][3][ri].v[v];
        x[1][0+c][ri].v[v] = sc[c][0][ri].v[v] - sc[c][2][ri].v[v];
        x[1][3+c][ri].v[v] = sc[c][1][ri].v[v] - sc[c][3][ri].v[v];
      }

   if (true1) {
    for (int v=0; v < VLENS; v++) {
      z0_000[v]  = a[0][ 0].v[v] * x[0][0][0].v[v];
      z0_001[v]  = a[0][ 0].v[v] * x[0][0][1].v[v];
      z0_010[v]  = a[1][ 0].v[v] * x[1][0][0].v[v];
      z0_011[v]  = a[1][ 0].v[v] * x[1][0][1].v[v];
      z0_100[v]  = a[0][ 1].v[v] * x[0][1][0].v[v];
      z0_101[v]  = a[0][ 1].v[v] * x[0][1][1].v[v];
      z0_110[v]  = a[1][ 1].v[v] * x[1][1][0].v[v];
      z0_111[v]  = a[1][ 1].v[v] * x[1][1][1].v[v];
      z1_000[v]  = 0;
      z1_001[v]  = 0;
      z1_010[v]  = 0;
      z1_011[v]  = 0;
      z1_100[v]  = 0;
      z1_101[v]  = 0;
      z1_110[v]  = 0;
      z1_111[v]  = 0;
    }
    for (int v=0; v < VLENS; v++) {
      z0_000[v] += a[0][ 6].v[v] * x[0][1][0].v[v];
      z1_000[v] -= a[0][ 7].v[v] * x[0][1][1].v[v];
      z0_001[v] += a[0][ 6].v[v] * x[0][1][1].v[v];
      z1_001[v] += a[0][ 7].v[v] * x[0][1][0].v[v];
      z0_010[v] += a[1][ 6].v[v] * x[1][1][0].v[v];
      z1_010[v] -= a[1][ 7].v[v] * x[1][1][1].v[v];
      z0_011[v] += a[1][ 6].v[v] * x[1][1][1].v[v];
      z1_011[v] += a[1][ 7].v[v] * x[1][1][0].v[v];
      z0_100[v] += a[0][ 6].v[v] * x[0][0][0].v[v];
      z1_100[v] += a[0][ 7].v[v] * x[0][0][1].v[v];
      z0_101[v] += a[0][ 6].v[v] * x[0][0][1].v[v];
      z1_101[v] -= a[0][ 7].v[v] * x[0][0][0].v[v];
      z0_110[v] += a[1][ 6].v[v] * x[1][0][0].v[v];
      z1_110[v] += a[1][ 7].v[v] * x[1][0][1].v[v];
      z0_111[v] += a[1][ 6].v[v] * x[1][0][1].v[v];
      z1_111[v] -= a[1][ 7].v[v] * x[1][0][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_000[v] += a[0][ 8].v[v] * x[0][2][0].v[v];
      z1_000[v] -= a[0][ 9].v[v] * x[0][2][1].v[v];
      z0_001[v] += a[0][ 8].v[v] * x[0][2][1].v[v];
      z1_001[v] += a[0][ 9].v[v] * x[0][2][0].v[v];
      z0_010[v] += a[1][ 8].v[v] * x[1][2][0].v[v];
      z1_010[v] -= a[1][ 9].v[v] * x[1][2][1].v[v];
      z0_011[v] += a[1][ 8].v[v] * x[1][2][1].v[v];
      z1_011[v] += a[1][ 9].v[v] * x[1][2][0].v[v];
      z0_100[v] += a[0][16].v[v] * x[0][2][0].v[v];
      z1_100[v] -= a[0][17].v[v] * x[0][2][1].v[v];
      z0_101[v] += a[0][16].v[v] * x[0][2][1].v[v];
      z1_101[v] += a[0][17].v[v] * x[0][2][0].v[v];
      z0_110[v] += a[1][16].v[v] * x[1][2][0].v[v];
      z1_110[v] -= a[1][17].v[v] * x[1][2][1].v[v];
      z0_111[v] += a[1][16].v[v] * x[1][2][1].v[v];
      z1_111[v] += a[1][17].v[v] * x[1][2][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_000[v] += a[0][10].v[v] * x[0][3][0].v[v];
      z1_000[v] -= a[0][11].v[v] * x[0][3][1].v[v];
      z0_001[v] += a[0][10].v[v] * x[0][3][1].v[v];
      z1_001[v] += a[0][11].v[v] * x[0][3][0].v[v];
      z0_010[v] += a[1][10].v[v] * x[1][3][0].v[v];
      z1_010[v] -= a[1][11].v[v] * x[1][3][1].v[v];
      z0_011[v] += a[1][10].v[v] * x[1][3][1].v[v];
      z1_011[v] += a[1][11].v[v] * x[1][3][0].v[v];
      z0_100[v] += a[0][18].v[v] * x[0][3][0].v[v];
      z1_100[v] -= a[0][19].v[v] * x[0][3][1].v[v];
      z0_101[v] += a[0][18].v[v] * x[0][3][1].v[v];
      z1_101[v] += a[0][19].v[v] * x[0][3][0].v[v];
      z0_110[v] += a[1][18].v[v] * x[1][3][0].v[v];
      z1_110[v] -= a[1][19].v[v] * x[1][3][1].v[v];
      z0_111[v] += a[1][18].v[v] * x[1][3][1].v[v];
      z1_111[v] += a[1][19].v[v] * x[1][3][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_000[v] += a[0][12].v[v] * x[0][4][0].v[v];
      z1_000[v] -= a[0][13].v[v] * x[0][4][1].v[v];
      z0_001[v] += a[0][12].v[v] * x[0][4][1].v[v];
      z1_001[v] += a[0][13].v[v] * x[0][4][0].v[v];
      z0_010[v] += a[1][12].v[v] * x[1][4][0].v[v];
      z1_010[v] -= a[1][13].v[v] * x[1][4][1].v[v];
      z0_011[v] += a[1][12].v[v] * x[1][4][1].v[v];
      z1_011[v] += a[1][13].v[v] * x[1][4][0].v[v];
      z0_100[v] += a[0][20].v[v] * x[0][4][0].v[v];
      z1_100[v] -= a[0][21].v[v] * x[0][4][1].v[v];
      z0_101[v] += a[0][20].v[v] * x[0][4][1].v[v];
      z1_101[v] += a[0][21].v[v] * x[0][4][0].v[v];
      z0_110[v] += a[1][20].v[v] * x[1][4][0].v[v];
      z1_110[v] -= a[1][21].v[v] * x[1][4][1].v[v];
      z0_111[v] += a[1][20].v[v] * x[1][4][1].v[v];
      z1_111[v] += a[1][21].v[v] * x[1][4][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_000[v] += a[0][14].v[v] * x[0][5][0].v[v];
      z1_000[v] -= a[0][15].v[v] * x[0][5][1].v[v];
      z0_001[v] += a[0][14].v[v] * x[0][5][1].v[v];
      z1_001[v] += a[0][15].v[v] * x[0][5][0].v[v];
      z0_010[v] += a[1][14].v[v] * x[1][5][0].v[v];
      z1_010[v] -= a[1][15].v[v] * x[1][5][1].v[v];
      z0_011[v] += a[1][14].v[v] * x[1][5][1].v[v];
      z1_011[v] += a[1][15].v[v] * x[1][5][0].v[v];
      z0_100[v] += a[0][22].v[v] * x[0][5][0].v[v];
      z1_100[v] -= a[0][23].v[v] * x[0][5][1].v[v];
      z0_101[v] += a[0][22].v[v] * x[0][5][1].v[v];
      z1_101[v] += a[0][23].v[v] * x[0][5][0].v[v];
      z0_110[v] += a[1][22].v[v] * x[1][5][0].v[v];
      z1_110[v] -= a[1][23].v[v] * x[1][5][1].v[v];
      z0_111[v] += a[1][22].v[v] * x[1][5][1].v[v];
      z1_111[v] += a[1][23].v[v] * x[1][5][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      sc[0][0][0].v[v] = (z0_000[v]+z1_000[v]) + (z0_010[v]+z1_010[v]);
      sc[0][0][1].v[v] = (z0_001[v]+z1_001[v]) + (z0_011[v]+z1_011[v]);
      sc[0][2][0].v[v] = (z0_000[v]+z1_000[v]) - (z0_010[v]+z1_010[v]);
      sc[0][2][1].v[v] = (z0_001[v]+z1_001[v]) - (z0_011[v]+z1_011[v]);
      sc[1][0][0].v[v] = (z0_100[v]+z1_100[v]) + (z0_110[v]+z1_110[v]);
      sc[1][0][1].v[v] = (z0_101[v]+z1_101[v]) + (z0_111[v]+z1_111[v]);
      sc[1][2][0].v[v] = (z0_100[v]+z1_100[v]) - (z0_110[v]+z1_110[v]);
      sc[1][2][1].v[v] = (z0_101[v]+z1_101[v]) - (z0_111[v]+z1_111[v]);
    }
   }

   if (true2) {
    for (int v=0; v < VLENS; v++) {
      z0_200[v]  = a[0][ 2].v[v] *  x[0][2][0].v[v];
      z0_201[v]  = a[0][ 2].v[v] *  x[0][2][1].v[v];
      z0_210[v]  = a[1][ 2].v[v] *  x[1][2][0].v[v];
      z0_211[v]  = a[1][ 2].v[v] *  x[1][2][1].v[v];
      z0_300[v]  = a[0][ 3].v[v] *  x[0][3][0].v[v];
      z0_301[v]  = a[0][ 3].v[v] *  x[0][3][1].v[v];
      z0_310[v]  = a[1][ 3].v[v] *  x[1][3][0].v[v];
      z0_311[v]  = a[1][ 3].v[v] *  x[1][3][1].v[v];
      z1_200[v]  = 0;
      z1_201[v]  = 0;
      z1_210[v]  = 0;
      z1_211[v]  = 0;
      z1_300[v]  = 0;
      z1_301[v]  = 0;
      z1_310[v]  = 0;
      z1_311[v]  = 0;
    }
    for (int v=0; v < VLENS; v++) {
      z0_200[v] += a[0][ 8].v[v] *  x[0][0][0].v[v];
      z1_200[v] += a[0][ 9].v[v] *  x[0][0][1].v[v];
      z0_201[v] += a[0][ 8].v[v] *  x[0][0][1].v[v];
      z1_201[v] -= a[0][ 9].v[v] *  x[0][0][0].v[v];
      z0_210[v] += a[1][ 8].v[v] *  x[1][0][0].v[v];
      z1_210[v] += a[1][ 9].v[v] *  x[1][0][1].v[v];
      z0_211[v] += a[1][ 8].v[v] *  x[1][0][1].v[v];
      z1_211[v] -= a[1][ 9].v[v] *  x[1][0][0].v[v];
      z0_300[v] += a[0][10].v[v] *  x[0][0][0].v[v];
      z1_300[v] += a[0][11].v[v] *  x[0][0][1].v[v];
      z0_301[v] += a[0][10].v[v] *  x[0][0][1].v[v];
      z1_301[v] -= a[0][11].v[v] *  x[0][0][0].v[v];
      z0_310[v] += a[1][10].v[v] *  x[1][0][0].v[v];
      z1_310[v] += a[1][11].v[v] *  x[1][0][1].v[v];
      z0_311[v] += a[1][10].v[v] *  x[1][0][1].v[v];
      z1_311[v] -= a[1][11].v[v] *  x[1][0][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_200[v] += a[0][16].v[v] *  x[0][1][0].v[v];
      z1_200[v] += a[0][17].v[v] *  x[0][1][1].v[v];
      z0_201[v] += a[0][16].v[v] *  x[0][1][1].v[v];
      z1_201[v] -= a[0][17].v[v] *  x[0][1][0].v[v];
      z0_210[v] += a[1][16].v[v] *  x[1][1][0].v[v];
      z1_210[v] += a[1][17].v[v] *  x[1][1][1].v[v];
      z0_211[v] += a[1][16].v[v] *  x[1][1][1].v[v];
      z1_211[v] -= a[1][17].v[v] *  x[1][1][0].v[v];
      z0_300[v] += a[0][18].v[v] *  x[0][1][0].v[v];
      z1_300[v] += a[0][19].v[v] *  x[0][1][1].v[v];
      z0_301[v] += a[0][18].v[v] *  x[0][1][1].v[v];
      z1_301[v] -= a[0][19].v[v] *  x[0][1][0].v[v];
      z0_310[v] += a[1][18].v[v] *  x[1][1][0].v[v];
      z1_310[v] += a[1][19].v[v] *  x[1][1][1].v[v];
      z0_311[v] += a[1][18].v[v] *  x[1][1][1].v[v];
      z1_311[v] -= a[1][19].v[v] *  x[1][1][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_200[v] += a[0][24].v[v] *  x[0][3][0].v[v];
      z1_200[v] -= a[0][25].v[v] *  x[0][3][1].v[v];
      z0_201[v] += a[0][24].v[v] *  x[0][3][1].v[v];
      z1_201[v] += a[0][25].v[v] *  x[0][3][0].v[v];
      z0_210[v] += a[1][24].v[v] *  x[1][3][0].v[v];
      z1_210[v] -= a[1][25].v[v] *  x[1][3][1].v[v];
      z0_211[v] += a[1][24].v[v] *  x[1][3][1].v[v];
      z1_211[v] += a[1][25].v[v] *  x[1][3][0].v[v];
      z0_300[v] += a[0][24].v[v] *  x[0][2][0].v[v];
      z1_300[v] += a[0][25].v[v] *  x[0][2][1].v[v];
      z0_301[v] += a[0][24].v[v] *  x[0][2][1].v[v];
      z1_301[v] -= a[0][25].v[v] *  x[0][2][0].v[v];
      z0_310[v] += a[1][24].v[v] *  x[1][2][0].v[v];
      z1_310[v] += a[1][25].v[v] *  x[1][2][1].v[v];
      z0_311[v] += a[1][24].v[v] *  x[1][2][1].v[v];
      z1_311[v] -= a[1][25].v[v] *  x[1][2][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_200[v] += a[0][26].v[v] *  x[0][4][0].v[v];
      z1_200[v] -= a[0][27].v[v] *  x[0][4][1].v[v];
      z0_201[v] += a[0][26].v[v] *  x[0][4][1].v[v];
      z1_201[v] += a[0][27].v[v] *  x[0][4][0].v[v];
      z0_210[v] += a[1][26].v[v] *  x[1][4][0].v[v];
      z1_210[v] -= a[1][27].v[v] *  x[1][4][1].v[v];
      z0_211[v] += a[1][26].v[v] *  x[1][4][1].v[v];
      z1_211[v] += a[1][27].v[v] *  x[1][4][0].v[v];
      z0_300[v] += a[0][30].v[v] *  x[0][4][0].v[v];
      z1_300[v] -= a[0][31].v[v] *  x[0][4][1].v[v];
      z0_301[v] += a[0][30].v[v] *  x[0][4][1].v[v];
      z1_301[v] += a[0][31].v[v] *  x[0][4][0].v[v];
      z0_310[v] += a[1][30].v[v] *  x[1][4][0].v[v];
      z1_310[v] -= a[1][31].v[v] *  x[1][4][1].v[v];
      z0_311[v] += a[1][30].v[v] *  x[1][4][1].v[v];
      z1_311[v] += a[1][31].v[v] *  x[1][4][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_200[v] += a[0][28].v[v] *  x[0][5][0].v[v];
      z1_200[v] -= a[0][29].v[v] *  x[0][5][1].v[v];
      z0_201[v] += a[0][28].v[v] *  x[0][5][1].v[v];
      z1_201[v] += a[0][29].v[v] *  x[0][5][0].v[v];
      z0_210[v] += a[1][28].v[v] *  x[1][5][0].v[v];
      z1_210[v] -= a[1][29].v[v] *  x[1][5][1].v[v];
      z0_211[v] += a[1][28].v[v] *  x[1][5][1].v[v];
      z1_211[v] += a[1][29].v[v] *  x[1][5][0].v[v];
      z0_300[v] += a[0][32].v[v] *  x[0][5][0].v[v];
      z1_300[v] -= a[0][33].v[v] *  x[0][5][1].v[v];
      z0_301[v] += a[0][32].v[v] *  x[0][5][1].v[v];
      z1_301[v] += a[0][33].v[v] *  x[0][5][0].v[v];
      z0_310[v] += a[1][32].v[v] *  x[1][5][0].v[v];
      z1_310[v] -= a[1][33].v[v] *  x[1][5][1].v[v];
      z0_311[v] += a[1][32].v[v] *  x[1][5][1].v[v];
      z1_311[v] += a[1][33].v[v] *  x[1][5][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      sc[2][0][0].v[v] = (z0_200[v]+z1_200[v]) +  (z0_210[v]+z1_210[v]);
      sc[2][0][1].v[v] = (z0_201[v]+z1_201[v]) +  (z0_211[v]+z1_211[v]);
      sc[2][2][0].v[v] = (z0_200[v]+z1_200[v]) -  (z0_210[v]+z1_210[v]);
      sc[2][2][1].v[v] = (z0_201[v]+z1_201[v]) -  (z0_211[v]+z1_211[v]);
      sc[0][1][0].v[v] = (z0_300[v]+z1_300[v]) +  (z0_310[v]+z1_310[v]);
      sc[0][1][1].v[v] = (z0_301[v]+z1_301[v]) +  (z0_311[v]+z1_311[v]);
      sc[0][3][0].v[v] = (z0_300[v]+z1_300[v]) -  (z0_310[v]+z1_310[v]);
      sc[0][3][1].v[v] = (z0_301[v]+z1_301[v]) -  (z0_311[v]+z1_311[v]);
    }
   }

   if (true3) {
    for (int v=0; v < VLENS; v++) {
      z0_400[v]  = a[0][ 4].v[v] *  x[0][4][0].v[v];
      z0_401[v]  = a[0][ 4].v[v] *  x[0][4][1].v[v];
      z0_410[v]  = a[1][ 4].v[v] *  x[1][4][0].v[v];
      z0_411[v]  = a[1][ 4].v[v] *  x[1][4][1].v[v];
      z0_500[v]  = a[0][ 5].v[v] *  x[0][5][0].v[v];
      z0_501[v]  = a[0][ 5].v[v] *  x[0][5][1].v[v];
      z0_510[v]  = a[1][ 5].v[v] *  x[1][5][0].v[v];
      z0_511[v]  = a[1][ 5].v[v] *  x[1][5][1].v[v];
      z1_400[v]  = 0;
      z1_401[v]  = 0;
      z1_410[v]  = 0;
      z1_411[v]  = 0;
      z1_500[v]  = 0;
      z1_501[v]  = 0;
      z1_510[v]  = 0;
      z1_511[v]  = 0;
    }
    for (int v=0; v < VLENS; v++) {
      z0_400[v] += a[0][12].v[v] *  x[0][0][0].v[v];
      z1_400[v] += a[0][13].v[v] *  x[0][0][1].v[v];
      z0_401[v] += a[0][12].v[v] *  x[0][0][1].v[v];
      z1_401[v] -= a[0][13].v[v] *  x[0][0][0].v[v];
      z0_410[v] += a[1][12].v[v] *  x[1][0][0].v[v];
      z1_410[v] += a[1][13].v[v] *  x[1][0][1].v[v];
      z0_411[v] += a[1][12].v[v] *  x[1][0][1].v[v];
      z1_411[v] -= a[1][13].v[v] *  x[1][0][0].v[v];
      z0_500[v] += a[0][14].v[v] *  x[0][0][0].v[v];
      z1_500[v] += a[0][15].v[v] *  x[0][0][1].v[v];
      z0_501[v] += a[0][14].v[v] *  x[0][0][1].v[v];
      z1_501[v] -= a[0][15].v[v] *  x[0][0][0].v[v];
      z0_510[v] += a[1][14].v[v] *  x[1][0][0].v[v];
      z1_510[v] += a[1][15].v[v] *  x[1][0][1].v[v];
      z0_511[v] += a[1][14].v[v] *  x[1][0][1].v[v];
      z1_511[v] -= a[1][15].v[v] *  x[1][0][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_400[v] += a[0][20].v[v] *  x[0][1][0].v[v];
      z1_400[v] += a[0][21].v[v] *  x[0][1][1].v[v];
      z0_401[v] += a[0][20].v[v] *  x[0][1][1].v[v];
      z1_401[v] -= a[0][21].v[v] *  x[0][1][0].v[v];
      z0_410[v] += a[1][20].v[v] *  x[1][1][0].v[v];
      z1_410[v] += a[1][21].v[v] *  x[1][1][1].v[v];
      z0_411[v] += a[1][20].v[v] *  x[1][1][1].v[v];
      z1_411[v] -= a[1][21].v[v] *  x[1][1][0].v[v];
      z0_500[v] += a[0][22].v[v] *  x[0][1][0].v[v];
      z1_500[v] += a[0][23].v[v] *  x[0][1][1].v[v];
      z0_501[v] += a[0][22].v[v] *  x[0][1][1].v[v];
      z1_501[v] -= a[0][23].v[v] *  x[0][1][0].v[v];
      z0_510[v] += a[1][22].v[v] *  x[1][1][0].v[v];
      z1_510[v] += a[1][23].v[v] *  x[1][1][1].v[v];
      z0_511[v] += a[1][22].v[v] *  x[1][1][1].v[v];
      z1_511[v] -= a[1][23].v[v] *  x[1][1][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_400[v] += a[0][26].v[v] *  x[0][2][0].v[v];
      z1_400[v] += a[0][27].v[v] *  x[0][2][1].v[v];
      z0_401[v] += a[0][26].v[v] *  x[0][2][1].v[v];
      z1_401[v] -= a[0][27].v[v] *  x[0][2][0].v[v];
      z0_410[v] += a[1][26].v[v] *  x[1][2][0].v[v];
      z1_410[v] += a[1][27].v[v] *  x[1][2][1].v[v];
      z0_411[v] += a[1][26].v[v] *  x[1][2][1].v[v];
      z1_411[v] -= a[1][27].v[v] *  x[1][2][0].v[v];
      z0_500[v] += a[0][28].v[v] *  x[0][2][0].v[v];
      z1_500[v] += a[0][29].v[v] *  x[0][2][1].v[v];
      z0_501[v] += a[0][28].v[v] *  x[0][2][1].v[v];
      z1_501[v] -= a[0][29].v[v] *  x[0][2][0].v[v];
      z0_510[v] += a[1][28].v[v] *  x[1][2][0].v[v];
      z1_510[v] += a[1][29].v[v] *  x[1][2][1].v[v];
      z0_511[v] += a[1][28].v[v] *  x[1][2][1].v[v];
      z1_511[v] -= a[1][29].v[v] *  x[1][2][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_400[v] += a[0][30].v[v] *  x[0][3][0].v[v];
      z1_400[v] += a[0][31].v[v] *  x[0][3][1].v[v];
      z0_401[v] += a[0][30].v[v] *  x[0][3][1].v[v];
      z1_401[v] -= a[0][31].v[v] *  x[0][3][0].v[v];
      z0_410[v] += a[1][30].v[v] *  x[1][3][0].v[v];
      z1_410[v] += a[1][31].v[v] *  x[1][3][1].v[v];
      z0_411[v] += a[1][30].v[v] *  x[1][3][1].v[v];
      z1_411[v] -= a[1][31].v[v] *  x[1][3][0].v[v];
      z0_500[v] += a[0][32].v[v] *  x[0][3][0].v[v];
      z1_500[v] += a[0][33].v[v] *  x[0][3][1].v[v];
      z0_501[v] += a[0][32].v[v] *  x[0][3][1].v[v];
      z1_501[v] -= a[0][33].v[v] *  x[0][3][0].v[v];
      z0_510[v] += a[1][32].v[v] *  x[1][3][0].v[v];
      z1_510[v] += a[1][33].v[v] *  x[1][3][1].v[v];
      z0_511[v] += a[1][32].v[v] *  x[1][3][1].v[v];
      z1_511[v] -= a[1][33].v[v] *  x[1][3][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      z0_400[v] += a[0][34].v[v] *  x[0][5][0].v[v];
      z1_400[v] -= a[0][35].v[v] *  x[0][5][1].v[v];
      z0_401[v] += a[0][34].v[v] *  x[0][5][1].v[v];
      z1_401[v] += a[0][35].v[v] *  x[0][5][0].v[v];
      z0_410[v] += a[1][34].v[v] *  x[1][5][0].v[v];
      z1_410[v] -= a[1][35].v[v] *  x[1][5][1].v[v];
      z0_411[v] += a[1][34].v[v] *  x[1][5][1].v[v];
      z1_411[v] += a[1][35].v[v] *  x[1][5][0].v[v];
      z0_500[v] += a[0][34].v[v] *  x[0][4][0].v[v];
      z1_500[v] += a[0][35].v[v] *  x[0][4][1].v[v];
      z0_501[v] += a[0][34].v[v] *  x[0][4][1].v[v];
      z1_501[v] -= a[0][35].v[v] *  x[0][4][0].v[v];
      z0_510[v] += a[1][34].v[v] *  x[1][4][0].v[v];
      z1_510[v] += a[1][35].v[v] *  x[1][4][1].v[v];
      z0_511[v] += a[1][34].v[v] *  x[1][4][1].v[v];
      z1_511[v] -= a[1][35].v[v] *  x[1][4][0].v[v];
    }
    for (int v=0; v < VLENS; v++) {
      sc[1][1][0].v[v] = (z0_400[v]+z1_400[v]) +  (z0_410[v]+z1_410[v]);
      sc[1][1][1].v[v] = (z0_401[v]+z1_401[v]) +  (z0_411[v]+z1_411[v]);
      sc[1][3][0].v[v] = (z0_400[v]+z1_400[v]) -  (z0_410[v]+z1_410[v]);
      sc[1][3][1].v[v] = (z0_401[v]+z1_401[v]) -  (z0_411[v]+z1_411[v]);
      sc[2][1][0].v[v] = (z0_500[v]+z1_500[v]) +  (z0_510[v]+z1_510[v]);
      sc[2][1][1].v[v] = (z0_501[v]+z1_501[v]) +  (z0_511[v]+z1_511[v]);
      sc[2][3][0].v[v] = (z0_500[v]+z1_500[v]) -  (z0_510[v]+z1_510[v]);
      sc[2][3][1].v[v] = (z0_501[v]+z1_501[v]) -  (z0_511[v]+z1_511[v]);
    }
   }
}