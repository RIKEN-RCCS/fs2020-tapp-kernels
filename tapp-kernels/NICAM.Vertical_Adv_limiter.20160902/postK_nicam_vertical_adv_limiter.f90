
!-------------------------------------------------------------------------------
!> module PRECISION
!!
!! @par Description
!!          precision module
!!
!! @author NICAM developers
!!
!<
module mod_precision
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: SP = kind(0.E0) ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = kind(0.D0) ! Double Precision: kind(0.D0)

#ifdef SINGLE
  integer, public, parameter :: RP = SP
#else
  integer, public, parameter :: RP = DP
#endif

  !-----------------------------------------------------------------------------
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module mod_precision
!-------------------------------------------------------------------------------
module mod_adm
  implicit none

  !====== Basic definition & information ======
  !
  !------ Log file ID & Control file ID
  integer, public, save      :: ADM_LOG_FID = 6 ! default is STDOUT
  !--- Identifier of triangle element (i-axis-side or j-axis side)
  integer, public, parameter :: ADM_TI = 1
  integer, public, parameter :: ADM_TJ = 2
  !
  !--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  integer, public, parameter :: ADM_AI  = 1
  integer, public, parameter :: ADM_AIJ = 2
  integer, public, parameter :: ADM_AJ  = 3
  !
  !------ Identifier of 1 variable
  integer, public, parameter :: ADM_KNONE = 1

  !
  !====== Information for processes ======
  !
  !------ My process ID
  integer, public, save      :: ADM_prc_me
  !
  logical, public, save      :: ADM_have_pl =.true.

  !
  !====== Information for processes-region relationship ======
  !

  !
  !====== Information for regions ======
  !------ Maximum number of pole regions
  integer, public, parameter :: ADM_rgn_nmax_pl = 2
  !
  !------ Local region number
  integer, public, save      :: ADM_lall = 1
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  !

  !
  !====== Grid resolution informations  ======
  !
  !------ Horizontal grid numbers
  integer, public, save      :: ADM_gmin = 2
  integer, public, save      :: ADM_gmax = 129 
  integer, public, save      :: ADM_gall_1d = 130
#ifdef SINGLE_SIM
  integer, public, save      :: ADM_gall_1d_in = 65
#else
  integer, public, save      :: ADM_gall_1d_in = 130
#endif
  integer, public, save      :: ADM_gall = 16900
  !
  !----- grid number of inner region in the diamond
  integer, public, save      :: ADM_gall_in
  !
  !------ Identifiers of grid points around poles.
  integer, public, parameter :: ADM_gslf_pl = 1
  integer, public, parameter :: ADM_gmin_pl = 2
  integer, public, save      :: ADM_gmax_pl = 6    ! [mod] S.Iga 100607
  integer, public, save      :: ADM_gall_pl = 6    ! [mod] S.Iga 100607
  !
  !------ Vertica grid numbers
  integer, public, save      :: ADM_kmin = 2
  integer, public, save      :: ADM_kmax = 95
  integer, public, save      :: ADM_kall = 96

end module mod_adm

module mod_cnst
  use mod_precision
  implicit none
  private

  !====== Misc. constants ======
  !
  !------ Allowable minimum value
  real(RP), public :: CNST_EPS_ZERO = 1.E-16_RP
  !
  !------ Allowable maximum value
  real(RP), public, parameter :: CNST_MAX_REAL = 1.E+30_RP
  !
end module mod_cnst
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!>
!! Tracer advection module
!!
!! @par Description
!!         This module contains subroutines for tracer advection
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
!-------------------------------------------------------------------------------

module mod_src_tracer
  use mod_precision
  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE

    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gall_1d_in,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax

  implicit none
  private


  public :: vertical_limiter_thuburn

contains


  subroutine vertical_limiter_thuburn( &
       q_h_ij, q_h_pl, &
       q_ij,   q_pl,   &
       d_ij,   d_pl,   &
       ck_ij,  ck_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    use mod_adm, only: ADM_prc_me

!    use omp_lib
    implicit none
    external omp_get_num_threads
    integer omp_get_num_threads

    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
!
    real(RP), intent(inout) :: q_h_ij   (ADM_gall_1d_in, ADM_kall, ADM_gall_1d,  ADM_lall   )
    real(RP), intent(in)    :: q_ij     (ADM_gall_1d_in, ADM_kall, ADM_gall_1d,  ADM_lall   )
    real(RP), intent(in)    :: d_ij     (ADM_gall_1d_in, ADM_kall, ADM_gall_1d,  ADM_lall   )
    real(RP), intent(in)    :: ck_ij    (ADM_gall_1d_in, ADM_kall, ADM_gall_1d,  ADM_lall   ,2)

    real(RP) :: Qin_min_pl (ADM_gall_pl,ADM_kall,2)
    real(RP) :: Qin_max_pl (ADM_gall_pl,ADM_kall,2)

    real(RP) :: Qout_min   (ADM_gall_1d_in,ADM_gall_1d,   ADM_kall)
    real(RP) :: Qout_max   (ADM_gall_1d_in,ADM_gall_1d,   ADM_kall)
    real(RP) :: Qout_min_1d
    real(RP) :: Qout_max_1d
    real(RP) :: Qout_min_m1   (ADM_gall_1d_in)
    real(RP) :: Qout_max_m1   (ADM_gall_1d_in)

    real(RP) :: qnext_min, qnext_min_pl
    real(RP) :: qnext_max, qnext_max_pl
    real(RP) :: Cin,       Cin_pl
    real(RP) :: Cout,      Cout_pl
    real(RP) :: CQin_min,  CQin_min_pl
    real(RP) :: CQin_max,  CQin_max_pl

    real(RP) :: mask, mask1, mask2
    real(RP) :: mask3,qin_min1,qin_min2,qin_max1,qin_max2
    real(RP) :: zerosw

    integer :: n, k, l
    integer :: i,j
    integer :: num_threads,nn,nstart,nend,onestep,onethread,mod_num
    !---------------------------------------------------------------------------

!cx    if (mark_for_snapshots(ip_Vertical_Adv_limiter).ne.0) then
!cx    call snapshot_seq_bin_open ("Vertical_Adv_limiter", ADM_prc_me)
!cx    call V_Adv_limiter_snap_write (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Vertical_Adv_limiter) = 0
!cx    endif

!$omp parallel default(none),  &
!$omp private(nn,onestep,onethread,mod_num,nstart,nend,num_threads, &
!$omp         mask1,mask3,qin_min1,qin_max1,qin_min2,qin_max2,  &
!$omp         qnext_min,qnext_max,mask2,Cin,Cout,CQin_max,CQin_min,zerosw, &
!$omp         Qout_max_m1,Qout_min_m1,Qout_max_1d,Qout_min_1d, &
!$omp         k,i,j,l), &
!$omp shared(q_h_ij,q_ij,d_ij,ck_ij , &
!$omp        CNST_EPS_ZERO,ADM_gall_1d_in,ADM_gall_1d,ADM_kmax,ADM_kall,ADM_lall)
    num_threads = omp_get_num_threads()
!$omp do
    do nn = 1,num_threads
    nstart = 1
    nend   = ADM_gall_1d

    !nstart1,nend1 loop1 / ADM_gmin-1,ADM_gmin-1
    onestep = (nend - nstart + 1)
    onethread = onestep / num_threads
    mod_num = mod(onestep,num_threads)

    if(nn <= mod_num) then
      nstart = (onethread + 1) * (nn-1) + nstart
      nend = nstart + onethread
    else
      nstart = (onethread + 1) * mod_num + (onethread * ((nn-1) - mod_num)) + nstart
      nend = nstart + onethread - 1
    endif

    do l = 1, ADM_lall
!fj 2016       do n = 1, ADM_gall
       do j = nstart, nend
        k=1
        do i = 1, ADM_gall_1d_in
          mask1 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,1)-CNST_EPS_ZERO)
          mask3 = 0.5_RP - sign(0.5_RP,ck_ij(i,k+1,j,l,1)-CNST_EPS_ZERO)
          qin_min1 = CNST_MAX_REAL
          qin_min2 = mask3*CNST_MAX_REAL        +min(q_ij(i,k,j,l),q_ij(i,k+1,j,l))
          qin_max1 = -CNST_MAX_REAL
          qin_max2 = -mask3*CNST_MAX_REAL       +max(q_ij(i,k,j,l),q_ij(i,k+1,j,l))
          qnext_min = min(qin_min2,q_ij(i,k,j,l))
          qnext_max = max(qin_max2,q_ij(i,k,j,l))
          mask2 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,2)-CNST_EPS_ZERO)
          Cin  = (      mask1 ) * ck_ij(i,k,j,l,1) &
               + (      mask2 ) * ck_ij(i,k,j,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck_ij(i,k,j,l,1) &
               + ( 1.0_RP-mask2 ) * ck_ij(i,k,j,l,2)
          CQin_max = mask1 * ( ck_ij(i,k,j,l,1) * qin_max1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck_ij(i,k,j,l,1) * qin_min1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_min2 )
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1
!fj 2016          Qout_min(i,k) = ( q(i,j,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_min_m1(i) = ( q_ij(i,k,j,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw
!fj 2106          Qout_max(i,k) = ( q(i,j,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_max_m1(i) = ( q_ij(i,k,j,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw
        enddo
!fj>
        do k = 2, ADM_kall-1
!fj 2016       do n = 1, ADM_gall
        do i = 1, ADM_gall_1d_in
!fj>
          mask1 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,1)-CNST_EPS_ZERO)
          mask3 = 0.5_RP - sign(0.5_RP,ck_ij(i,k+1,j,l,1)-CNST_EPS_ZERO)
          qin_min1 = (1.0_RP-mask1)*CNST_MAX_REAL +min(q_ij(i,k,j,l),q_ij(i,k-1,j,l))
          qin_min2 = mask3*CNST_MAX_REAL        +min(q_ij(i,k,j,l),q_ij(i,k+1,j,l))
          qin_max1 = -(1.0_RP-mask1)*CNST_MAX_REAL+max(q_ij(i,k,j,l),q_ij(i,k-1,j,l))
          qin_max2 = -mask3*CNST_MAX_REAL       +max(q_ij(i,k,j,l),q_ij(i,k+1,j,l))
!fj<
          qnext_min = min(qin_min1,qin_min2,q_ij(i,k,j,l))
          qnext_max = max(qin_max1,qin_max2,q_ij(i,k,j,l))
!fj<

          mask2 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,2)-CNST_EPS_ZERO)

          Cin  = (      mask1 ) * ck_ij(i,k,j,l,1) &
               + (      mask2 ) * ck_ij(i,k,j,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck_ij(i,k,j,l,1) &
               + ( 1.0_RP-mask2 ) * ck_ij(i,k,j,l,2)

          CQin_max = mask1 * ( ck_ij(i,k,j,l,1) * qin_max1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck_ij(i,k,j,l,1) * qin_min1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_min2 )

          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1

!fj 2016          Qout_min(i,k) = ( q(i,j,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_min_1d = ( q_ij(i,k,j,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw
!fj 2016          Qout_max(i,k) = ( q(i,j,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_max_1d = ( q_ij(i,k,j,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw
          q_h_ij(i,k,j,l) = min( max( q_h_ij(i,k,j,l), min(q_ij(i,k,j,l),q_ij(i,k-1,j,l)) ), max(q_ij(i,k,j,l),q_ij(i,k-1,j,l)) ) 
          q_h_ij(i,k,j,l) = (      mask1 ) * max( min( q_h_ij(i,k,j,l), Qout_max_m1(i) ), Qout_min_m1(i) ) &
                     + ( 1.0_RP-mask1 ) * max( min( q_h_ij(i,k,j,l), Qout_max_1d), Qout_min_1d )
          Qout_max_m1(i) = Qout_max_1d
          Qout_min_m1(i) = Qout_min_1d
        enddo
       enddo
!fj>
!fj 2016       k=1
!fj 2016       do n = 1, ADM_gall
!fj 2016          mask1 = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)
!fj 2016          mask3 = 0.5_RP - sign(0.5_RP,ck(n,k+1,l,1)-CNST_EPS_ZERO)
!fj 2016          qin_min2 = mask3*CNST_MAX_REAL        +min(q(n,k,l),q(n,k+1,l))
!fj 2016          qin_max2 = -mask3*CNST_MAX_REAL       +max(q(n,k,l),q(n,k+1,l))
!fj 2016          qnext_min = min(qin_min2,q(n,k,l))
!fj 2016          qnext_max = max(qin_max2,q(n,k,l))
!fj 2016          mask2 = 0.5_RP - sign(0.5_RP,ck(n,k,l,2)-CNST_EPS_ZERO)
!fj 2016          Cin  = (      mask1 ) * ck(n,k,l,1) &
!fj 2016               + (      mask2 ) * ck(n,k,l,2)
!fj 2016          Cout = ( 1.0_RP-mask1 ) * ck(n,k,l,1) &
!fj 2016               + ( 1.0_RP-mask2 ) * ck(n,k,l,2)
!fj 2016          CQin_max = mask1 * ( ck(n,k,l,1) * qin_max1 ) &
!fj 2016                   + mask2 * ( ck(n,k,l,2) * qin_max2 )
!fj 2016          CQin_min = mask1 * ( ck(n,k,l,1) * qin_min1 ) &
!fj 2016                   + mask2 * ( ck(n,k,l,2) * qin_min2 )
!fj 2016          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1
!fj 2016          Qout_min(n,k) = ( q(n,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
!fj 2016                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
!fj 2016                        + q(n,k,l) * zerosw
!fj 2016          Qout_max(n,k) = ( q(n,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
!fj 2016                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
!fj 2016                        + q(n,k,l) * zerosw
!fj 2016       enddo
       k=ADM_kmax+1
!fj 2016       do n = 1, ADM_gall
        do i = 1, ADM_gall_1d_in
          mask1 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,1)-CNST_EPS_ZERO)
          qin_min1 = (1.0_RP-mask1)*CNST_MAX_REAL +min(q_ij(i,k,j,l),q_ij(i,k-1,j,l))
          qin_max1 = -(1.0_RP-mask1)*CNST_MAX_REAL+max(q_ij(i,k,j,l),q_ij(i,k-1,j,l))
          qnext_min = min(qin_min1,q_ij(i,k,j,l))
          qnext_max = max(qin_max1,q_ij(i,k,j,l))
          mask2 = 0.5_RP - sign(0.5_RP,ck_ij(i,k,j,l,2)-CNST_EPS_ZERO)
          Cin  = (      mask1 ) * ck_ij(i,k,j,l,1) &
               + (      mask2 ) * ck_ij(i,k,j,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck_ij(i,k,j,l,1) &
               + ( 1.0_RP-mask2 ) * ck_ij(i,k,j,l,2)
          CQin_max = mask1 * ( ck_ij(i,k,j,l,1) * qin_max1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck_ij(i,k,j,l,1) * qin_min1 ) &
                   + mask2 * ( ck_ij(i,k,j,l,2) * qin_min2 )
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1
!fj 2016          Qout_min(i,k) = ( q(i,j,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_min_1d = ( q_ij(i,k,j,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw
!fj 2016          Qout_max(i,k) = ( q(i,j,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(i,j,k,l)) ) &
          Qout_max_1d = ( q_ij(i,k,j,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d_ij(i,k,j,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q_ij(i,k,j,l) * zerosw

          q_h_ij(i,k,j,l) = min( max( q_h_ij(i,k,j,l), min(q_ij(i,k,j,l),q_ij(i,k-1,j,l)) ), max(q_ij(i,k,j,l),q_ij(i,k-1,j,l)) ) 
          q_h_ij(i,k,j,l) = (      mask1 ) * max( min( q_h_ij(i,k,j,l), Qout_max_m1(i) ), Qout_min_m1(i) ) &
                     + ( 1.0_RP-mask1 ) * max( min( q_h_ij(i,k,j,l), Qout_max_1d), Qout_min_1d )
        enddo
       enddo
!fj<

!fj 2016       do k = ADM_kmin, ADM_kmax+1
!fj 2016       do n = 1, ADM_gall
!fj 2016          mask = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)

!fj          q_h(n,k,l) = (      mask ) * min( max( q_h(n,k,l), Qin_min(n,k  ,1) ), Qin_max(n,k  ,1) ) &
!fj                     + ( 1._RP-mask ) * min( max( q_h(n,k,l), Qin_min(n,k-1,2) ), Qin_max(n,k-1,2) )
!fj 2016          q_h(n,k,l) = min( max( q_h(n,k,l), min(q(n,k,l),q(n,k-1,l)) ), max(q(n,k,l),q(n,k-1,l)) ) 

!fj 2016          q_h(n,k,l) = (      mask ) * max( min( q_h(n,k,l), Qout_max(n,k-1) ), Qout_min(n,k-1) ) &
!fj 2016                     + ( 1.0_RP-mask ) * max( min( q_h(n,k,l), Qout_max(n,k  ) ), Qout_min(n,k  ) )
!fj 2016       enddo
!fj 2016       enddo

    enddo
    enddo
!$omp end parallel


    return
  end subroutine vertical_limiter_thuburn

end module mod_src_tracer
