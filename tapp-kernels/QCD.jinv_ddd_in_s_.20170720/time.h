#ifndef TIMI_H
#define TIMI_H

#if defined(_FPCOLL)
///////////////////////////////
// Uses Fujitsu's monitor 
///////////////////////////////
#  include <fjcoll.h>
#  define _TIC_ start_collection(__func__)
#  define _TOC_  stop_collection(__func__)
#  define _BCG_ITER_TIC_ start_collection("bicgstab_iter_")
#  define _BCG_ITER_TOC_  stop_collection("bicgstab_iter_")

#  define _BCG_DDS_ITER_TIC_ start_collection("bicgstab_dd_s_iter_")
#  define _BCG_DDS_ITER_TOC_  stop_collection("bicgstab_dd_s_iter_")

#  define _BCG_PRECDDS_TIC_ start_collection("bicgstab_precdd_s_")
#  define _BCG_PRECDDS_TOC_  stop_collection("bicgstab_precdd_s_")
#  define _BCG_PRECDDS_ITER_TIC_ start_collection("bicgstab_precdd_s_iter_")
#  define _BCG_PRECDDS_ITER_TOC_  stop_collection("bicgstab_precdd_s_iter_")
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_ start_collection("bicgstab_precdd_s_iter_reduc1_")
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_  stop_collection("bicgstab_precdd_s_iter_reduc1_")
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_ start_collection("bicgstab_precdd_s_iter_reduc2_")
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_  stop_collection("bicgstab_precdd_s_iter_reduc2_")
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_ start_collection("bicgstab_precdd_s_iter_reduc3_")
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_  stop_collection("bicgstab_precdd_s_iter_reduc3_")
#  define _MCG_ITER_TIC_ start_collection("mcg_iter_")
#  define _MCG_ITER_TOC_  stop_collection("mcg_iter_")
#  define _CG_ITER_TIC_ start_collection("cg_iter_")
#  define _CG_ITER_TOC_  stop_collection("cg_iter_")

#  define _DEO_IN_TIC_ start_collection("deo_in_")
#  define _DEO_IN_TOC_  stop_collection("deo_in_")

#  define _DEE_DEO_IN_TIC_ start_collection("dee_deo_in_")
#  define _DEE_DEO_IN_TOC_  stop_collection("dee_deo_in_")

#  define _DEO_OUT_PRE_TIC_ start_collection("deo_out_pre_")
#  define _DEO_OUT_PRE_TOC_  stop_collection("deo_out_pre_")

#  define _DEE_DEO_OUT_POST_TIC_ start_collection("dee_deo_out_post_")
#  define _DEE_DEO_OUT_POST_TOC_  stop_collection("dee_deo_out_post_")

#  define _MTILDE_TIC_ start_collection("mtilde")
#  define _MTILDE_TOC_  stop_collection("mtilde")

#  define _PREC_S_TIC_ start_collection("prec_s_")
#  define _PREC_S_TOC_  stop_collection("prec_s_")
#  define _PREC_DDD_S_TIC_ start_collection("prec_ddd_s_")
#  define _PREC_DDD_S_TOC_  stop_collection("prec_ddd_s_")

#  define _DDD_IN_S_TIC_      start_collection("ddd_in_s_")
#  define _DDD_IN_S_TOC_       stop_collection("ddd_in_s_")
#  define _JINV_DDD_IN_S_TIC_ start_collection("jinv_ddd_in_s_")
#  define _JINV_DDD_IN_S_TOC_  stop_collection("jinv_ddd_in_s_")
#  define _DDD_OUT_PRE_S_TIC_ start_collection("ddd_out_pre_s_")
#  define _DDD_OUT_PRE_S_TOC_  stop_collection("ddd_out_pre_s_")
#  define _DDD_OUT_POS_S_TIC_ start_collection("ddd_out_pos_s_")
#  define _DDD_OUT_POS_S_TOC_  stop_collection("ddd_out_pos_s_")

#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_ start_collection("s_mult_wd_deo_out_recv_hpc_calc_")
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_ stop_collection("s_mult_wd_deo_out_recv_hpc_calc_")
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_ start_collection("s_mult_wd_deo_out_send_hpc_")
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_ stop_collection("s_mult_wd_deo_out_send_hpc_")

#  define _OHTER_CALC_TIC_
#  define _OHTER_CALC_TOC_

#  define _COMM_TIC_
#  define _COMM_TOC_

#elif defined(_CHECK_TIMING)
///////////////////////////////
// Uses costom timing monitor
///////////////////////////////
#  ifndef _TIMING_H_
#    define _TIMING_H_
void check_timing_ (char *);
#  endif  /* _TIMING_H_ */

#  define _TIC_ check_timing_(__func__)
#  define _TOC_ check_timing_(__func__)
#  define _BCG_ITER_TIC_ check_timing_("bicgstab_iter_")
#  define _BCG_ITER_TOC_ check_timing_("bicgstab_iter_")
#  define _BCG_DDS_ITER_TIC_ check_timing_("bicgstab_dd_s_iter_")
#  define _BCG_DDS_ITER_TOC_ check_timing_("bicgstab_dd_s_iter_")
#  define _BCG_PRECDDS_TIC_ check_timing_("bicgstab_precdd_s_")
#  define _BCG_PRECDDS_TOC_ check_timing_("bicgstab_precdd_s_")
#  define _BCG_PRECDDS_ITER_TIC_ check_timing_("bicgstab_precdd_s_iter_")
#  define _BCG_PRECDDS_ITER_TOC_ check_timing_("bicgstab_precdd_s_iter_")
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_ check_timing_("bicgstab_precdd_s_iter_reduc1_")
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_ check_timing_("bicgstab_precdd_s_iter_reduc1_")
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_ check_timing_("bicgstab_precdd_s_iter_reduc2_")
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_ check_timing_("bicgstab_precdd_s_iter_reduc2_")
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_ check_timing_("bicgstab_precdd_s_iter_reduc3_")
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_ check_timing_("bicgstab_precdd_s_iter_reduc3_")
#  define _MCG_ITER_TIC_ check_timing_("mcg_iter_")
#  define _MCG_ITER_TOC_ check_timing_("mcg_iter_")
#  define _CG_ITER_TIC_ check_timing_("cg_iter_")
#  define _CG_ITER_TOC_ check_timing_("cg_iter_")

#  define _DEO_IN_TIC_ check_timing_("deo_in_")
#  define _DEO_IN_TOC_ check_timing_("deo_in_")

#  define _DEE_DEO_IN_TIC_ check_timing_("dee_deo_in_")
#  define _DEE_DEO_IN_TOC_ check_timing_("dee_deo_in_")

#  define _DEO_OUT_PRE_TIC_ check_timing_("deo_out_pre_")
#  define _DEO_OUT_PRE_TOC_ check_timing_("deo_out_pre_")

#  define _DEE_DEO_OUT_POST_TIC_ check_timing_("dee_deo_out_post_")
#  define _DEE_DEO_OUT_POST_TOC_ check_timing_("dee_deo_out_post_")

#  define _MTILDE_TIC_ check_timing_("mtilde")
#  define _MTILDE_TOC_ check_timing_("mtilde")

#  define _PREC_S_TIC_ check_timing_("prec_s_")
#  define _PREC_S_TOC_ check_timing_("prec_s_")
#  define _PREC_DDD_S_TIC_ check_timing_("prec_ddd_s_")
#  define _PREC_DDD_S_TOC_ check_timing_("prec_ddd_s_")

#  define _DDD_IN_S_TIC_      check_timing_("ddd_in_s_")
#  define _DDD_IN_S_TOC_      check_timing_("ddd_in_s_")
#  define _JINV_DDD_IN_S_TIC_ check_timing_("jinv_ddd_in_s_")
#  define _JINV_DDD_IN_S_TOC_ check_timing_("jinv_ddd_in_s_")
#  define _DDD_OUT_PRE_S_TIC_ check_timing_("ddd_out_pre_s_")
#  define _DDD_OUT_PRE_S_TOC_ check_timing_("ddd_out_pre_s_")
#  define _DDD_OUT_POS_S_TIC_ check_timing_("ddd_out_pos_s_")
#  define _DDD_OUT_POS_S_TOC_ check_timing_("ddd_out_pos_s_")

#  define  _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_  check_timing_("s_mult_wd_deo_out_recv_hpc_calc_")
#  define  _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_  check_timing_("s_mult_wd_deo_out_recv_hpc_calc_")
#  define  _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_       check_timing_("s_mult_wd_deo_out_send_hpc_")
#  define  _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_       check_timing_("s_mult_wd_deo_out_send_hpc_")
#  define _OHTER_CALC_TIC_
#  define _OHTER_CALC_TOC_
#  define _COMM_TIC_
#  define _COMM_TOC_

#elif defined(_CHECK_PA)
#  define _TIC_
#  define _TOC_ 
#  define _BCG_ITER_TIC_
#  define _BCG_ITER_TOC_ 
#  define _BCG_DDS_ITER_TIC_
#  define _BCG_DDS_ITER_TOC_ 
#  define _BCG_PRECDDS_TIC_
#  define _BCG_PRECDDS_TOC_ 
#  define _BCG_PRECDDS_ITER_TIC_
#  define _BCG_PRECDDS_ITER_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_ 
#  define _MCG_ITER_TIC_
#  define _MCG_ITER_TOC_ 
#  define _CG_ITER_TIC_
#  define _CG_ITER_TOC_ 
#  define _DEO_IN_TIC_
#  define _DEO_IN_TOC_ 
#  define _DEE_DEO_IN_TIC_
#  define _DEE_DEO_IN_TOC_ 
#  define _DEO_OUT_PRE_TIC_
#  define _DEO_OUT_PRE_TOC_ 
#  define _DEE_DEO_OUT_POST_TIC_
#  define _DEE_DEO_OUT_POST_TOC_ 
#  define _MTILDE_TIC_
#  define _MTILDE_TOC_ 
#  define _PREC_S_TIC_
#  define _PREC_S_TOC_ 
#  define _PREC_DDD_S_TIC_
#  define _PREC_DDD_S_TOC_ 
#  define _JINV_DDD_IN_S_TIC_
#  define _JINV_DDD_IN_S_TOC_
#  define _DDD_IN_S_TIC_
#  define _DDD_IN_S_TOC_ 
#  define _DDD_OUT_PRE_S_TIC_
#  define _DDD_OUT_PRE_S_TOC_ 
#  define _DDD_OUT_POS_S_TIC_
#  define _DDD_OUT_POS_S_TOC_
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
#  define _COMM_TIC_
#  define _COMM_TOC_
#  define _OHTER_CALC_TIC_
#  define _OHTER_CALC_TOC_

#  include"profiler.h"

#  ifdef TARGET_JINV
#  define _JINV_DDD_IN_S_TIC_                   PROF_START_SRL
#  define _JINV_DDD_IN_S_TOC_                   PROF_STOP_SRL
#  elif defined(TARGET_IN)
#  define _DDD_IN_S_TIC_                        PROF_START_SRL
#  define _DDD_IN_S_TOC_                        PROF_STOP_SRL
#  elif defined(TARGET_PRE)
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_      PROF_START_SRL
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_      PROF_STOP_SRL
#  elif defined(TARGET_POS)
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_ PROF_START_SRL
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_ PROF_STOP_SRL
#  elif defined(TARGET_OTHER)
#  define _BCG_PRECDDS_ITER_TIC_                PROF_START
#  define _BCG_PRECDDS_ITER_TOC_                PROF_STOP
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_         PROF_STOP
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_         PROF_START
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_         PROF_STOP
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_         PROF_START
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_         PROF_STOP 
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_         PROF_START
#  define _OHTER_CALC_TIC_                      PROF_START_SRL
#  define _OHTER_CALC_TOC_                      PROF_STOP_SRL
#  endif

#elif defined(_CHECK_SIM)
#  define _TIC_
#  define _TOC_ 
#  define _BCG_ITER_TIC_
#  define _BCG_ITER_TOC_ 
#  define _BCG_DDS_ITER_TIC_
#  define _BCG_DDS_ITER_TOC_ 
#  define _BCG_PRECDDS_TIC_
#  define _BCG_PRECDDS_TOC_ 
#  define _MCG_ITER_TIC_
#  define _MCG_ITER_TOC_ 
#  define _CG_ITER_TIC_
#  define _CG_ITER_TOC_ 
#  define _DEO_IN_TIC_
#  define _DEO_IN_TOC_ 
#  define _DEE_DEO_IN_TIC_
#  define _DEE_DEO_IN_TOC_ 
#  define _DEO_OUT_PRE_TIC_
#  define _DEO_OUT_PRE_TOC_ 
#  define _DEE_DEO_OUT_POST_TIC_
#  define _DEE_DEO_OUT_POST_TOC_ 
#  define _MTILDE_TIC_
#  define _MTILDE_TOC_ 
#  define _PREC_S_TIC_
#  define _PREC_S_TOC_ 
#  define _PREC_DDD_S_TIC_
#  define _PREC_DDD_S_TOC_ 
#  define _DDD_OUT_PRE_S_TIC_
#  define _DDD_OUT_PRE_S_TOC_ 
#  define _DDD_OUT_POS_S_TIC_
#  define _DDD_OUT_POS_S_TOC_
#  define _COMM_TIC_
#  define _COMM_TOC_

#  include"profiler.h"

# ifdef TARGET_JINV
#  define _JINV_DDD_IN_S_TIC_                   if(prof_flag==1)PROF_START("name")
#  define _JINV_DDD_IN_S_TOC_                   if(prof_flag==1)PROF_STOP("name")
# else
#  define _JINV_DDD_IN_S_TIC_
#  define _JINV_DDD_IN_S_TOC_
# endif

# ifdef TARGET_IN
#  define _DDD_IN_S_TIC_                        if(prof_flag==1)PROF_START("name")
#  define _DDD_IN_S_TOC_                        if(prof_flag==1)PROF_STOP("name")
# else
#  define _DDD_IN_S_TIC_
#  define _DDD_IN_S_TOC_ 
# endif

# ifdef TARGET_PRE
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_      if(prof_flag==1)PROF_START("name")
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_      if(prof_flag==1)PROF_STOP("name")
# else
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
# endif

# ifdef TARGET_POS
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_ if(prof_flag==1)PROF_START("name")
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_ if(prof_flag==1)PROF_STOP("name")
# else
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_
# endif

# ifdef TARGET_OTHER
#  define _BCG_PRECDDS_ITER_TIC_                if(prof_flag==1)PROF_START("name")
#  define _BCG_PRECDDS_ITER_TOC_                if(prof_flag==1)PROF_STOP("name")
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_         if(prof_flag==1)PROF_START("name")
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_         if(prof_flag==1)PROF_STOP("name")
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_         if(prof_flag==1)PROF_START("name")
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_         if(prof_flag==1)PROF_STOP("name")
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_         if(prof_flag==1)PROF_START("name")
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_         if(prof_flag==1)PROF_STOP("name")
#  define _OHTER_CALC_TIC_                      if(prof_flag==1)PROF_START("name")
#  define _OHTER_CALC_TOC_                      if(prof_flag==1)PROF_STOP("name")
# else 
#  define _BCG_PRECDDS_ITER_TIC_
#  define _BCG_PRECDDS_ITER_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_ 
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_ 
#  define _OHTER_CALC_TIC_
#  define _OHTER_CALC_TOC_
# endif

#else

#  define _TIC_
#  define _TOC_
#  define _BCG_ITER_TIC_
#  define _BCG_ITER_TOC_
#  define _BCG_DDS_ITER_TIC_
#  define _BCG_DDS_ITER_TOC_
#  define _BCG_PRECDDS_TIC_
#  define _BCG_PRECDDS_TOC_
#  define _BCG_PRECDDS_ITER_TIC_
#  define _BCG_PRECDDS_ITER_TOC_
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_
#  define _MCG_ITER_TIC_
#  define _MCG_ITER_TOC_
#  define _CG_ITER_TIC_
#  define _CG_ITER_TOC_

#  define _DEO_IN_TIC_
#  define _DEO_IN_TOC_
#  define _DEE_DEO_IN_TIC_
#  define _DEE_DEO_IN_TOC_
#  define _DEO_OUT_PRE_TIC_
#  define _DEO_OUT_PRE_TOC_
#  define _DEE_DEO_OUT_POST_TIC_
#  define _DEE_DEO_OUT_POST_TOC_

#  define _MTILDE_TIC_
#  define _MTILDE_TOC_

#  define _PREC_S_TIC_
#  define _PREC_S_TOC_
#  define _PREC_DDD_S_TIC_
#  define _PREC_DDD_S_TOC_

#  define _DDD_IN_S_TIC_      
#  define _DDD_IN_S_TOC_      
#  define _JINV_DDD_IN_S_TIC_ 
#  define _JINV_DDD_IN_S_TOC_ 
#  define _DDD_OUT_PRE_S_TIC_ 
#  define _DDD_OUT_PRE_S_TOC_ 
#  define _DDD_OUT_POS_S_TIC_ 
#  define _DDD_OUT_POS_S_TOC_ 

#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
#  define _COMM_TIC_
#  define _COMM_TOC_
#  define _OHTER_CALC_TIC_
#  define _OHTER_CALC_TOC_

#endif

#endif
