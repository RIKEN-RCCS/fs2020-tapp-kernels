#ifndef _MULT_ALL_H
#define _MULT_ALL_H

#define _ALIGN_SIZE  (4)

//
// link multiplication on two-spinor
// upy1 = u * y1
// upy2 = u * y2
//
#define __mult_u_y_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u).c[0][c][0][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u).c[0][c][0][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[1][c][0][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[1][c][0][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[2][c][0][j] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[2][c][0][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] -= (u).c[0][c][1][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[0][c][1][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] -= (u).c[1][c][1][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[1][c][1][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] -= (u).c[2][c][1][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[2][c][1][j] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

//
// link multiplication on two-spinor on a site with simd-site last site link
// upy1 = u * y1
// upy2 = u * y2
//
#define __mult_u_y_3_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    (upy).c[c][s][0]  = (u).c[0][c][0][VLENS-1] * (py).c[0][s][0];\
    (upy).c[c][s][1]  = (u).c[0][c][0][VLENS-1] * (py).c[0][s][1];\
    (upy).c[c][s][0] += (u).c[1][c][0][VLENS-1] * (py).c[1][s][0];\
    (upy).c[c][s][1] += (u).c[1][c][0][VLENS-1] * (py).c[1][s][1];\
    (upy).c[c][s][0] += (u).c[2][c][0][VLENS-1] * (py).c[2][s][0];\
    (upy).c[c][s][1] += (u).c[2][c][0][VLENS-1] * (py).c[2][s][1];\
    (upy).c[c][s][0] -= (u).c[0][c][1][VLENS-1] * (py).c[0][s][1];\
    (upy).c[c][s][1] += (u).c[0][c][1][VLENS-1] * (py).c[0][s][0];\
    (upy).c[c][s][0] -= (u).c[1][c][1][VLENS-1] * (py).c[1][s][1];\
    (upy).c[c][s][1] += (u).c[1][c][1][VLENS-1] * (py).c[1][s][0];\
    (upy).c[c][s][0] -= (u).c[2][c][1][VLENS-1] * (py).c[2][s][1];\
    (upy).c[c][s][1] += (u).c[2][c][1][VLENS-1] * (py).c[2][s][0];\
  }\
  }\
}



//
// link multiplication on two-spinor
// upy1 = udag * y1
// upy2 = udag * y2
//
#define __mult_udag_y_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u).c[c][0][0][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u).c[c][0][0][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][1][0][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[c][1][0][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][2][0][j] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[c][2][0][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][0][1][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][0][1][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] += (u).c[c][1][1][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][1][1][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] += (u).c[c][2][1][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][2][1][j] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

//
// link multiplication on two-spinor with simd-site-shifted gauge link
// upy1(j) = udag(j+VLEN-1) * y1(j)
// upy2(j) = udag(j+VLEN-1) * y2(j)
//
#define __mult_udag_y_2_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u)[c][0][0][j+VLENS-1] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u)[c][0][0][j+VLENS-1] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][1][0][j+VLENS-1] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u)[c][1][0][j+VLENS-1] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][2][0][j+VLENS-1] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u)[c][2][0][j+VLENS-1] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][0][1][j+VLENS-1] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][0][1][j+VLENS-1] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] += (u)[c][1][1][j+VLENS-1] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][1][1][j+VLENS-1] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] += (u)[c][2][1][j+VLENS-1] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][2][1][j+VLENS-1] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

//
// link multiplication on two-spinor on a site with the last simd-site gauge
// upy1 = udag(VLEN-1) * y1
// upy2 = udag(VLEN-1) * y2
//
#define __mult_udag_y_3_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
      (upy).c[c][s][0]  = (u).c[c][0][0][VLENS-1] * (py).c[0][s][0];\
      (upy).c[c][s][1]  = (u).c[c][0][0][VLENS-1] * (py).c[0][s][1];\
      (upy).c[c][s][0] += (u).c[c][1][0][VLENS-1] * (py).c[1][s][0];\
      (upy).c[c][s][1] += (u).c[c][1][0][VLENS-1] * (py).c[1][s][1];\
      (upy).c[c][s][0] += (u).c[c][2][0][VLENS-1] * (py).c[2][s][0];\
      (upy).c[c][s][1] += (u).c[c][2][0][VLENS-1] * (py).c[2][s][1];\
      (upy).c[c][s][0] += (u).c[c][0][1][VLENS-1] * (py).c[0][s][1];\
      (upy).c[c][s][1] -= (u).c[c][0][1][VLENS-1] * (py).c[0][s][0];\
      (upy).c[c][s][0] += (u).c[c][1][1][VLENS-1] * (py).c[1][s][1];\
      (upy).c[c][s][1] -= (u).c[c][1][1][VLENS-1] * (py).c[1][s][0];\
      (upy).c[c][s][0] += (u).c[c][2][1][VLENS-1] * (py).c[2][s][1];\
      (upy).c[c][s][1] -= (u).c[c][2][1][VLENS-1] * (py).c[2][s][0];\
  }\
  }\
}

//
// X-forward spin pre projection
// py1 = y1 - i y4
// py2 = y2 - i y3
//
#define __mult_x_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][3][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][3][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][2][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][2][0][j];\
    }\
  }\
}

//
// X-forward spin pre projection with forward simd-site-shifting
// py1(j) = y1(j+1) - i y4(j+1)
// py2(j) = y2(j+1) - i y3(j+1)
//
#define __mult_x_forw_pre_2_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y)[c][0][0][j+1] + (y)[c][3][1][j+1];\
      (py).c[c][0][1][j] = (y)[c][0][1][j+1] - (y)[c][3][0][j+1];\
      (py).c[c][1][0][j] = (y)[c][1][0][j+1] + (y)[c][2][1][j+1];\
      (py).c[c][1][1][j] = (y)[c][1][1][j+1] - (y)[c][2][0][j+1];\
    }\
  }\
}

//
// X-forward spin pre projection with first simd-site spionr
// py1 = y1(0) - i y4(0)
// py2 = y2(0) - i y3(0)
//
#define __mult_x_forw_pre_3_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    (py).c[c][0][0] = (y).c[c][0][0][0] + (y).c[c][3][1][0];\
    (py).c[c][0][1] = (y).c[c][0][1][0] - (y).c[c][3][0][0];\
    (py).c[c][1][0] = (y).c[c][1][0][0] + (y).c[c][2][1][0];\
    (py).c[c][1][1] = (y).c[c][1][1][0] - (y).c[c][2][0][0];\
  }\
}

//
// X-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 + i upy2
// my4 = my4 + i upy1
//
#define __mult_x_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] = (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] = (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] = (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] = (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] = -(upy).c[c][1][1][j];\
      (my).c[c][2][1][j] = (upy).c[c][1][0][j];\
      (my).c[c][3][0][j] = -(upy).c[c][0][1][j];\
      (my).c[c][3][1][j] = (upy).c[c][0][0][j];\
    }\
  }\
}

//
// X-forward spin post reconstruction at simd-site-shift last site
// my1(VLEN-1) = my1(VLEN-1) +   upy1
// my2(VLEN-1) = my2(VLEN-1) +   upy2
// my3(VLEN-1) = my3(VLEN-1) + i upy2
// my4(VLEN-1) = my4(VLEN-1) + i upy1
//
#define __mult_x_forw_pst_3_(my,upy) {\
  for (int c = 0; c < 3; ++c) {\
    (my).c[c][0][0][VLENS-1] += (upy).c[c][0][0];\
    (my).c[c][0][1][VLENS-1] += (upy).c[c][0][1];\
    (my).c[c][1][0][VLENS-1] += (upy).c[c][1][0];\
    (my).c[c][1][1][VLENS-1] += (upy).c[c][1][1];\
    (my).c[c][2][0][VLENS-1] -= (upy).c[c][1][1];\
    (my).c[c][2][1][VLENS-1] += (upy).c[c][1][0];\
    (my).c[c][3][0][VLENS-1] -= (upy).c[c][0][1];\
    (my).c[c][3][1][VLENS-1] += (upy).c[c][0][0];\
  }\
}


//
// X-backward spin pre projection
// py1 = y1 + i y4
// py2 = y2 + i y3
//
#define __mult_x_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][3][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][3][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][2][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][2][0][j];\
    }\
  }\
}

//
// X-backward spin pre projection with backward simd-site-shifti(ng
// py1(j) = y1(j+VLEN-1) + i y4(j+VLEN-1)
// py2(j) = y2(j+VLEN-1) + i y3(j+VLEN-1)
//
#define __mult_x_back_pre_2_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y)[c][0][0][j+VLENS-1] - (y)[c][3][1][j+VLENS-1];\
      (py).c[c][0][1][j] = (y)[c][0][1][j+VLENS-1] + (y)[c][3][0][j+VLENS-1];\
      (py).c[c][1][0][j] = (y)[c][1][0][j+VLENS-1] - (y)[c][2][1][j+VLENS-1];\
      (py).c[c][1][1][j] = (y)[c][1][1][j+VLENS-1] + (y)[c][2][0][j+VLENS-1];\
    }\
  }\
}

//
// X-backward spin pre projection with last simd-site spinor
// py1 = y1(VLEN-1) + i y4(VLEN-1)
// py2 = y2(VLEN-1) + i y3(VLEN-1)
//
#define __mult_x_back_pre_3_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    (py).c[c][0][0] = (y).c[c][0][0][VLENS-1] - (y).c[c][3][1][VLENS-1];\
    (py).c[c][0][1] = (y).c[c][0][1][VLENS-1] + (y).c[c][3][0][VLENS-1];\
    (py).c[c][1][0] = (y).c[c][1][0][VLENS-1] - (y).c[c][2][1][VLENS-1];\
    (py).c[c][1][1] = (y).c[c][1][1][VLENS-1] + (y).c[c][2][0][VLENS-1];\
  }\
}

//
// X-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 - i upy2
// my4 = my4 - i upy1
//
#define __mult_x_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][1][j] -= (upy).c[c][1][0][j];\
      (my).c[c][3][0][j] += (upy).c[c][0][1][j];\
      (my).c[c][3][1][j] -= (upy).c[c][0][0][j];\
    }\
  }\
}

//
// X-backward spin post reconstruction at simd-site-sift first site
// my1(0) = my1(0) +   upy1
// my2(0) = my2(0) +   upy2
// my3(0) = my3(0) - i upy2
// my4(0) = my4(0) - i upy1
//
#define __mult_x_back_pst_3_(my,upy) {\
  for (int c = 0; c < 3; ++c) {\
    (my).c[c][0][0][0] += (upy).c[c][0][0];\
    (my).c[c][0][1][0] += (upy).c[c][0][1];\
    (my).c[c][1][0][0] += (upy).c[c][1][0];\
    (my).c[c][1][1][0] += (upy).c[c][1][1];\
    (my).c[c][2][0][0] += (upy).c[c][1][1];\
    (my).c[c][2][1][0] -= (upy).c[c][1][0];\
    (my).c[c][3][0][0] += (upy).c[c][0][1];\
    (my).c[c][3][1][0] -= (upy).c[c][0][0];\
  }\
}

//
// Y-forward spin pre projection
// py1 = y1 - y4
// py2 = y2 + y3
//
#define __mult_y_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][3][0][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][3][1][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][2][0][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][2][1][j];\
    }\
  }\
}

//
// Y-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 +   upy2
// my4 = my4 -   upy1
//
#define __mult_y_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][2][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][3][0][j] -= (upy).c[c][0][0][j];\
      (my).c[c][3][1][j] -= (upy).c[c][0][1][j];\
    }\
  }\
}

//
// Y-backward spin pre projection
// py1 = y1 + y4
// py2 = y2 - y3
//
#define __mult_y_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][3][0][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][3][1][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][2][0][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][2][1][j];\
    }\
  }\
}

//
// Y-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 -   upy2
// my4 = my4 +   upy1
//
#define __mult_y_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] -= (upy).c[c][1][0][j];\
      (my).c[c][2][1][j] -= (upy).c[c][1][1][j];\
      (my).c[c][3][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][3][1][j] += (upy).c[c][0][1][j];\
    }\
  }\
}

//
// Z-forward spin pre projection
// py1 = y1 - i y3
// py2 = y2 + i y4
//
#define __mult_z_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][2][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][2][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][3][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][3][0][j];\
    }\
  }\
}

//
// Z-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 + i upy1
// my4 = my4 - i upy2
//
#define __mult_z_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] -= (upy).c[c][0][1][j];\
      (my).c[c][2][1][j] += (upy).c[c][0][0][j];\
      (my).c[c][3][0][j] += (upy).c[c][1][1][j];\
      (my).c[c][3][1][j] -= (upy).c[c][1][0][j];\
    }\
  }\
}

//
// Z-backward spin pre projection
// py1 = y1 + i y3
// py2 = y2 - i y4
//
#define __mult_z_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][2][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][2][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][3][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][3][0][j];\
    }\
  }\
}

//
// Z-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 - i upy1
// my4 = my4 + i upy2
//
#define __mult_z_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][0][1][j];\
      (my).c[c][2][1][j] -= (upy).c[c][0][0][j];\
      (my).c[c][3][0][j] -= (upy).c[c][1][1][j];\
      (my).c[c][3][1][j] += (upy).c[c][1][0][j];\
    }\
  }\
}

//
// T-forward spin pre projection
// py1 = 2*y3
// py2 = 2*y4
//
#define __mult_t_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = 2*(y).c[c][2][0][j];\
      (py).c[c][0][1][j] = 2*(y).c[c][2][1][j];\
      (py).c[c][1][0][j] = 2*(y).c[c][3][0][j];\
      (py).c[c][1][1][j] = 2*(y).c[c][3][1][j];\
    }\
  }\
}

//
// T-forward spin pre projection with boundary condition
// py1 = fbc *y3
// py2 = fbc *y4
//
// fbc = 2 * phase
//
#define __mult_t_forw_pre_bc_(py,y,fbc) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (fbc)*(y).c[c][2][0][j];\
      (py).c[c][0][1][j] = (fbc)*(y).c[c][2][1][j];\
      (py).c[c][1][0][j] = (fbc)*(y).c[c][3][0][j];\
      (py).c[c][1][1][j] = (fbc)*(y).c[c][3][1][j];\
    }\
  }\
}

//
// Y-forward spin post reconstruction
// my1 = my1
// my2 = my2
// my3 = my3 +   upy1
// my4 = my4 +   upy2
//
#define __mult_t_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][2][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][2][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][3][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][3][1][j] += (upy).c[c][1][1][j];\
    }\
  }\
}

//
// Y-forward spin post reconstruction with bounday-condition
// my1 = my1
// my2 = my2
// my3 = my3 +   upy1 * fbc
// my4 = my4 +   upy2 * fbc
//
#define __mult_t_forw_pst_bc_(my,upy,fbc) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][2][0][j] += (upy).c[c][0][0][j] * (fbc);\
      (my).c[c][2][1][j] += (upy).c[c][0][1][j] * (fbc);\
      (my).c[c][3][0][j] += (upy).c[c][1][0][j] * (fbc);\
      (my).c[c][3][1][j] += (upy).c[c][1][1][j] * (fbc);\
    }\
  }\
}


//
// T-backward spin pre projection
// py1 = 2*y1
// py2 = 2*y2
//
#define __mult_t_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = 2*(y).c[c][0][0][j];\
      (py).c[c][0][1][j] = 2*(y).c[c][0][1][j];\
      (py).c[c][1][0][j] = 2*(y).c[c][1][0][j];\
      (py).c[c][1][1][j] = 2*(y).c[c][1][1][j];\
    }\
  }\
}

//
// T-backward spin pre projection with boundary condition
// py1 = fbc*y1
// py2 = fbc*y2
//
// fbc = 2 * phase
//
#define __mult_t_back_pre_bc_(py,y,fbc) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (fbc)*(y).c[c][0][0][j];\
      (py).c[c][0][1][j] = (fbc)*(y).c[c][0][1][j];\
      (py).c[c][1][0][j] = (fbc)*(y).c[c][1][0][j];\
      (py).c[c][1][1][j] = (fbc)*(y).c[c][1][1][j];\
    }\
  }\
}

//
// Y-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3
// my4 = my4
//
#define __mult_t_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
    }\
  }\
}

//
// Y-backward spin post reconstruction with boundary-condition
// my1 = my1 +   upy1 * fbc
// my2 = my2 +   upy2 * fbc
// my3 = my3
// my4 = my4
//
#define __mult_t_back_pst_bc_(my,upy,fbc) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j] * (fbc);\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j] * (fbc);\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j] * (fbc);\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j] * (fbc);\
    }\
  }\
}

#endif
