
//
// X-forward (top end)
//
#define _mult_forw_x_recv_ {\
\
  int ix = y + ny*z + ny*nz*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs1_t ua;\
\
  __mult_u_y_3_(ua,(*(xfs_recv + ix)),(*(gx + i0)));\
  __mult_x_forw_pst_3_(tmp,ua);\
\
}

//
// X-backward (bottom end)
//
#define _mult_back_x_recv_ {\
\
  int ix = y + ny*z + ny*nz*t;\
\
  projscs1_t *xbs_recvi = xbs_recv + ix;\
\
  __mult_x_back_pst_3_(tmp,(*xbs_recvi));\
\
}

//
// Y-forward  (top end)
//
#define _mult_forw_y_recv_ {\
\
  int i2 = x +nxs*z +nxs*nz*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
\
  __mult_u_y_(ua,(*(yfs_recv + i2)),(*(gy + i0)));\
  __mult_y_forw_pst_(tmp,ua);\
\
}

//
// Y-backward (bottom end)
//
#define _mult_back_y_recv_ {\
\
  int i3 = x +nxs*z +nxs*nz*t;\
  projscs_t *ybs_recvi = ybs_recv + i3;\
\
  __mult_y_back_pst_(tmp,(*ybs_recvi));\
\
}

//
// Z-forward (top end)
//
#define _mult_forw_z_recv_ {\
\
  int i4 = x +nxs*y +nxs*ny*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
\
  __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + i0)));\
  __mult_z_forw_pst_(tmp,ua);\
\
}

//
// Z-backward (bottom end)
//
#define _mult_back_z_recv_ {\
\
  int i5 = x +nxs*y +nxs*ny*t;\
  projscs_t *zbs_recvi = zbs_recv + i5;\
\
  __mult_z_back_pst_(tmp,(*zbs_recvi));\
\
}

//
// T-forward (top end)
//
#define _mult_forw_t_recv_ {\
\
  int i6 = x +nxs*y +nxs*ny*z;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
  float tbc_fwd = ((float)fbc[3][0])*0.5f;\
\
  __mult_u_y_(ua,(*(tfs_recv + i6)),(*(gt + i0)));\
  __mult_t_forw_pst_bc_(tmp,ua,tbc_fwd);\
\
}

//
// T-backward (bottom end)
//
#define _mult_back_t_recv_ {\
\
  int i7 = x +nxs*y +nxs*ny*z;\
  projscs_t *tbs_recvi = tbs_recv + i7;\
  float tbc_bwd = ((float)fbc[3][1])*0.5f;\
\
  __mult_t_back_pst_bc_(tmp,(*tbs_recvi),tbc_bwd);\
\
}

//
// mult clover and accumulate 
//
#define _mult_clover_and_accumulate_recv_ {\
\
  __mult_clvs( tmp.cv, clvs[i0 + vols*(*idomain)].cv);\
\
  scs_t *outi = out + i0;\
  for (int c = 0; c < 3; c++) {\
  for (int s = 0; s < 4; s++) {\
    for (int j = 0; j < VLENS; j++) {\
      outi->c[c][s][0][j] += tmp.c[c][s][0][j] * factor;\
      outi->c[c][s][1][j] += tmp.c[c][s][1][j] * factor;\
    }\
  }\
  }\
}
