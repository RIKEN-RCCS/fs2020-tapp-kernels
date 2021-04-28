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
module mod_cnst
  use mod_precision
  implicit none
  private
  !++ Public procedure

  !++ Public parameters & variables
  real(RP), public, save :: CNST_ERADIUS = 6.37122E6_RP ! Radius of the Earth [m]
  real(RP), public, save :: CNST_EOHM    = 7.292E-5_RP   ! Angular velocity of the Earth [/s]
  real(RP), public, save :: CNST_EGRAV   = 9.80616_RP  ! Gravitational accerlaration of the Earth [m/s2]

  real(RP), public, save :: CNST_RAIR    =  287.0_RP   ! Gas constant of air
  real(RP), public, save :: CNST_RVAP    =  461.5_RP   ! Gas constant of vapor

  real(RP), public, save :: CNST_CP      = 1004.5_RP   ! Specific heat of air (consant pressure)
  real(RP), public, save :: CNST_CV                   ! Specific heat of air (consant volume)

  real(RP), public, save :: CNST_CPV     = 1846.0_RP   ! Specific heat of vapor (consant pressure)
  real(RP), public, save :: CNST_CVV                  ! Specific heat of vapor (consant volume)
  real(RP), public, save :: CNST_CL      = 4218.0_RP   ! Specific heat of water
  real(RP), public, save :: CNST_CI      = 2006.0_RP   ! Specific heat of ice

  !------ Density of water
  real(RP), public, save :: CNST_DWATR = 1000.0_RP
  !
  !------ Saturate pressure of water vapor at 0C
  real(RP), public, save :: CNST_PSAT0 = 610.7_RP
  !<----- unit : [Pa]
  !
  !------ Latent heat of vaporizaion at 0C
!  real(RP), public, save :: CNST_LH0   = 2.5008D+6 [mod] 20120704 H.Yashiro
  real(RP), public, save :: CNST_LH0   = 2.501E+6_RP
  !
  !------ Latent heat of sublimation at 0C
!  real(RP), public, save :: CNST_LHS0  = 2.8342E+6 [mod] 20120704 H.Yashiro
  real(RP), public, save :: CNST_LHS0  = 2.834E+6_RP
  !
  !------ Latent heat of melting
  real(RP), public, save :: CNST_EMELT = 3.40E+5_RP
  !
  !------ Melting temperature of water
  real(RP), public, save :: CNST_TMELT = 273.15_RP
  !
  !------ Freeze point of sea
  real(RP), public, save :: CNST_TFRZS  = 271.35_RP
  !
  !------ Wet-bulb temp. rain/snow
  real(RP), public, save :: CNST_TQICE = 273.15_RP
  !
  !------ Stefan-Boltzman constant
  real(RP), public, save :: CNST_STB   = 5.67E-8_RP
  !
  !------ Karman constant
  real(RP), public, save :: CNST_KARMAN = 0.4_RP
  !
  !------ Surface pressure
  real(RP), public, save :: CNST_PRES0    = 101325.0_RP
  !
  !------ Surface temperature
  real(RP), public, save :: CNST_TEMS0    = 300.0_RP
  !
  !------ Standard pressure
  real(RP), public, save :: CNST_PRE00    = 1.0E+5_RP
  !
  !------ Standard temperature
  real(RP), public, save :: CNST_TEM00    = 273.15_RP
  !
  !====== Misc. constants ======
  !
  !------ Definition of PI
  real(RP), public, save :: CNST_PI = 3.14159265358979323846_RP

  real(RP), public, save :: CNST_D2R

  !------ Allowable minimum value
  real(RP), public :: CNST_EPS_ZERO = 1.E-16_RP
  !
  !------ Allowable maximum value
  real(RP), public, parameter :: CNST_MAX_REAL = 1.E+30_RP
  !
  !------ Missing value
  real(RP), public, parameter :: CNST_VMISS    = 0.0_RP
  !
  !------ Undefined value
  real(DP), public, parameter :: CNST_UNDEF    = -99.9D+33
  !
  !------ Undefined value
  real(4), public, parameter :: CNST_UNDEF4   = -99.9E+33
  !
  !------ Undefined value
  integer(4), public, parameter :: CNST_UNDEF2   = -32768

  !-----------------------------------------------------------------------------

end module mod_cnst
!-------------------------------------------------------------------------------
module mod_adm
  implicit none

  !------ Character length of system control
  integer, public, parameter :: ADM_NSYS = 32
  !
  !------ Maximum length of file name
  integer, public, parameter :: ADM_MAXFNAME = 128
  !
  !====== Basic definition & information ======
  !
  !------ Log file ID & Control file ID
  integer, public, save      :: ADM_LOG_FID = 6 ! default is STDOUT
  integer, public, parameter :: ADM_CTL_FID = 35
  !
  !------ Identifier for single computation or parallel computation
  integer, public, parameter :: ADM_SINGLE_PRC = 0
  integer, public, parameter :: ADM_MULTI_PRC  = 1
  !
  !------ Identifiers of directions of region edges
  integer, public, parameter :: ADM_SW = 1
  integer, public, parameter :: ADM_NW = 2
  integer, public, parameter :: ADM_NE = 3
  integer, public, parameter :: ADM_SE = 4
  !
  !------ Identifiers of directions of region vertices
  integer, public, parameter :: ADM_W = 1
  integer, public, parameter :: ADM_N = 2
  integer, public, parameter :: ADM_E = 3
  integer, public, parameter :: ADM_S = 4
  !
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
  integer, public, parameter :: ADM_VNONE = 1
  !
  !------ Identifier of poles (north pole or south pole)
  integer, public, parameter :: ADM_NPL = 1
  integer, public, parameter :: ADM_SPL = 2
  !
  !------ Fist colomn on the table for region and direction
  integer, public, parameter :: ADM_RID = 1
  integer, public, parameter :: ADM_DIR = 2
  !
  real(8), public, parameter :: ADM_VMISS = 1.D0

  !
  !====== Information for processes ======
  !
  !------ Communication world for NICAM
  integer, public, save      :: ADM_COMM_WORLD
  logical, public, save      :: ADM_MPI_alive = .false.
  !
  !------ Master process
  integer, public, parameter :: ADM_prc_run_master = 1
  !
  !------ Total number of process
  integer, public, save      :: ADM_prc_all
  !
  !------ My process ID
  integer, public, save      :: ADM_prc_me
  !
  !------ Process ID which manages the pole regions.
  integer, public, save      :: ADM_prc_pl
  !
  !------ Process ID which have the pole regions.
  integer, public, save      :: ADM_prc_npl
  integer, public, save      :: ADM_prc_spl
  integer, public, save      :: ADM_prc_nspl(ADM_NPL:ADM_SPL)
  logical, public, save      :: ADM_have_pl = .true.

  !
  !====== Information for processes-region relationship ======
  !

  !
  !====== Information for regions ======
  !
  !------ Region division level
  integer, public, save      :: ADM_rlevel
  !
  !------ Total number of regular regions managed by all process
  integer, public, save      :: ADM_rgn_nmax
  !
  !------ Maximum number of pole regions
  integer, public, parameter :: ADM_rgn_nmax_pl = 2
  !
  !------ Local region number
  integer, public, save      :: ADM_lall = 1
  integer, public, parameter :: ADM_lall_para = 1
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  !
  !------ Present Local region number ! 2010.4.26 M.Satoh
  integer, public, save      :: ADM_l_me

  logical, public, save      :: ADM_have_sgp(ADM_lall_para) = .true. ! region have singlar point?

  !
  !====== Grid resolution informations  ======
  !
  !------ Grid division level
  integer, public, save      :: ADM_glevel
  !
  !------ Horizontal grid numbers
  integer, public, save      :: ADM_gmin = 2
  integer, public, save      :: ADM_gmax = 129
  integer, public, save      :: ADM_gall_1d = 130
#ifdef SINGLE_SIM
  integer, public, save      :: ADM_gall_1d_in = 65
  integer, public, save      :: ADM_gall = 16900/2
#else
  integer, public, save      :: ADM_gall_1d_in = 130
  integer, public, save      :: ADM_gall = 16900
#endif
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
  integer, public, save      :: ADM_vlayer
  integer, public, save      :: ADM_kmin = 2
!  integer, public, save      :: ADM_kmax = 95
!  integer, public, save      :: ADM_kall = 96
  integer, public, save      :: ADM_kmax = 47
  integer, public, save      :: ADM_kall = 48

end module mod_adm

module mod_runconf
  integer, public :: TRC_vmax   =  0 ! total number of tracers
  logical, public :: THUBURN_LIM = .true.  ! [add] 20130613 R.Yoshida
end module mod_runconf
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
       ADM_have_sgp, &
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

    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    use mod_adm, only: ADM_prc_me
  implicit none
  private

  public :: horizontal_limiter_thuburn
!  public :: H_Adv_limiter_snap_write
!  public :: H_Adv_limiter_snap_read


contains

  !-----------------------------------------------------------------------------
  !> Miura(2004)'s scheme with Thuburn(1996) limiter
  subroutine horizontal_limiter_thuburn( &
       q_a,    q_a_pl,  &
!fj       q,      q_pl,    &
       q_tmp,      q_pl,    &
       d,      d_pl,    &
       ch,     ch_pl,   &
       cmask,  cmask_pl )

!cx    use mod_comm, only: &
!cx       COMM_data_transfer

    implicit none

!fj    real(RP), intent(inout) :: q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: q       (  ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
    real(RP), intent(in)    :: q_tmp   (  ADM_gall   ,ADM_kall,ADM_lall   )
!fj<
    real(RP), intent(in)    :: q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
!fj    real(RP), intent(in)    :: cmask   (3,ADM_gall   ,ADM_kall,ADM_lall   )
!fj<
!!!    real(RP), intent(inout) :: q_a     (6,ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_a     (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall,6   )
    real(RP), intent(in)    :: d       (  ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
!!!    real(RP), intent(in)    :: ch      (6,ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: ch      (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall,6   )
    real(RP), intent(in)    :: cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!!    real(RP), intent(in)    :: cmask   (3,ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: cmask   (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall, 3   )

    real(RP) :: q_min_AI, q_min_AIJ, q_min_AJ, q_min_pl
    real(RP) :: q_max_AI, q_max_AIJ, q_max_AJ, q_max_pl

!fj    real(8) :: qnext_min   (ADM_gall), qnext_min_pl
!fj    real(8) :: qnext_max   (ADM_gall), qnext_max_pl
!fj    real(8) :: Cin_sum     (ADM_gall), Cin_sum_pl
!fj    real(8) :: Cout_sum    (ADM_gall), Cout_sum_pl
!fj    real(8) :: CQin_max_sum(ADM_gall), CQin_max_sum_pl
!fj    real(8) :: CQin_min_sum(ADM_gall), CQin_min_sum_pl
    !fj>
    real(RP) :: qnext_min   , qnext_min_pl
    real(RP) :: qnext_max   , qnext_max_pl
    real(RP) :: Cin_sum     , Cin_sum_pl
    real(RP) :: Cout_sum    , Cout_sum_pl
    real(RP) :: CQin_max_sum, CQin_max_sum_pl
    real(RP) :: CQin_min_sum, CQin_min_sum_pl
!fj    real(RP) :: q       (  -ADM_gall_1d+1:ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: q       (0:ADM_gall_1d_in+1,0:ADM_gall_1d+1   ,ADM_kall,ADM_lall   )
    !fj<

    integer, parameter :: I_min = 1
    integer, parameter :: I_max = 2
!fj    real(8) :: Qin    (6,ADM_gall   ,ADM_kall,ADM_lall   ,2)
!fj>
!fj    real(RP) :: Qin1    (ADM_gall      ,2)
!fj    real(RP) :: Qin2    (ADM_gall      ,2)
!fj    real(RP) :: Qin3    (ADM_gall      ,2)
!fj    real(RP) :: Qin4    (ADM_gall      ,2)
!fj    real(RP) :: Qin5    (ADM_gall      ,2)
    real(RP) :: Qin6    (ADM_gall_1d_in,ADM_gall_1d      ,2)
    real(RP) :: Qin1    (ADM_gall_1d_in,ADM_gall_1d      ,2)
    real(RP) :: Qin2    (ADM_gall_1d_in,ADM_gall_1d      ,2)
    real(RP) :: Qin3    (ADM_gall_1d_in,ADM_gall_1d      ,2)
    real(RP) :: Qin4    (ADM_gall_1d_in,ADM_gall_1d      ,2)
    real(RP) :: Qin5    (ADM_gall_1d_in,ADM_gall_1d      ,2)
!fj<
    real(RP) :: Qin_pl (2,ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
!fj    real(RP) :: Qout   (  ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(RP) :: Qout_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(RP) :: Qout   (  ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   ,2)

    real(RP) :: zerosw
    !fj>
    real(RP) :: tmp,cmaskch1,cmaskch2,cmaskch3,cmaskch4,cmaskch5,cmaskch6
    !fj<

!fj    integer :: ij
!fj    integer :: ip1j, ijp1, ip1jp1, ip2jp1
!fj    integer :: im1j, ijm1

!fj    integer :: nstart, nend
!fj    integer :: n, k, l, v
    integer :: k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d_in * ((j)-1) + (i)
    !---------------------------------------------------------------------------

!    call DEBUG_rapstart('____Horizontal_Adv_limiter')

!cx    if (mark_for_snapshots(ip_Horizontal_Adv_limiter).ne.0) then
!cx    call snapshot_seq_bin_open ("Horizontal_Adv_limiter", ADM_prc_me)
!cx    call H_Adv_limiter_snap_write &
!cx                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Horizontal_Adv_limiter) = 0
!cx    endif

!    call DEBUG_rapstart('loop1')
!    call timer_sta(2)

    !---< (i) define inflow bounds, eq.(32)&(33) >---
!OCL PREFETCH_SEQUENTIAL(SOFT),PREFETCH_CACHE_LEVEL(1),PREFETCH_STRONG
    do l = 1, ADM_lall
!$omp parallel
!$omp do
    do k = 1, ADM_kall

    j = 0
    do i = 0,ADM_gall_1d_in+1
    q(i,j,k,l) = 0.0_RP
    enddo

    j = ADM_gall_1d+1
    do i = 0,ADM_gall_1d_in+1
    q(i,j,k,l) = 0.0_RP
    enddo

    i = 0
    do j = 0,ADM_gall_1d+1
    q(i,j,k,l) = 0.0_RP
    enddo

    i = ADM_gall_1d+1
    do j = 0,ADM_gall_1d+1
    q(i,j,k,l) = 0.0_RP
    enddo

    do j = 1, ADM_gall_1d
    do i = 1, ADM_gall_1d_in
    q(i,j,k,l) = q_tmp(suf(i,j),k,l)
    enddo
    enddo

    j = 1 
    do i = 1, ADM_gall_1d_in
    Qout(i,j,k,l,I_min) = 0.0_RP
    Qout(i,j,k,l,I_max) = 0.0_RP
    enddo

    j = ADM_gmax+1 
    do i = 1, ADM_gall_1d_in
    Qout(i,j,k,l,I_min) = 0.0_RP
    Qout(i,j,k,l,I_max) = 0.0_RP
    enddo

    i = 1 
    do j = 1, ADM_gall_1d
    Qout(i,j,k,l,I_min) = 0.0_RP
    Qout(i,j,k,l,I_max) = 0.0_RP
    enddo

    i = ADM_gmax+1 
    do j = 1, ADM_gall_1d
    Qout(i,j,k,l,I_min) = 0.0_RP
    Qout(i,j,k,l,I_max) = 0.0_RP
    enddo
    enddo
!$omp end do
!$omp end parallel
    enddo

!    call timer_end(2)
    !call DEBUG_rapend  ('loop1')

    !call DEBUG_rapstart('loop2')
!    call timer_sta(3)

!OCL PREFETCH_SEQUENTIAL(SOFT),PREFETCH_CACHE_LEVEL(1),PREFETCH_STRONG
    do l = 1, ADM_lall
!$omp parallel private(q_min_AI,q_max_AI,q_min_AIJ, &
!$omp          q_max_AIJ,q_min_AJ,q_max_AJ,i,j,k,&
!$omp          cmaskch1,cmaskch2,cmaskch3,cmaskch4, &
!$omp          cmaskch5,cmaskch6,Cin_sum,Cout_sum,CQin_min_sum,CQin_max_sum, &
!$omp          qnext_min,qnext_max,tmp,zerosw)
    do k = 1, ADM_kall
!$omp do
       do j = ADM_gmin-1,ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin-1,ADM_gmax/2
#else
       do i = ADM_gmin-1,ADM_gmax
#endif
!OCL PREFETCH_READ(q(i,j,k+2,l),level=2,strong=1)
!OCL PREFETCH_READ(q(i,j+1,k+2,l),level=2,strong=1)
!!!!OCL PREFETCH_READ(cmask(1,i,j,k+2,l),level=2,strong=1)
          q_min_AI  = min( q(i,j,k,l), q(i,j-1,k,l), q(i+1,j,k,l), q(i+1,j+1,k,l) )
          q_max_AI  = max( q(i,j,k,l), q(i,j-1,k,l), q(i+1,j,k,l), q(i+1,j+1,k,l) )
          q_min_AIJ = min( q(i,j,k,l), q(i+1,j,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l) )
          q_max_AIJ = max( q(i,j,k,l), q(i+1,j,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l) )
          q_min_AJ  = min( q(i,j,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l), q(i-1,j,k,l) )
          q_max_AJ  = max( q(i,j,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l), q(i-1,j,k,l) )

          Qin1(i,j,    I_min) = (      cmask(i,j,k,l,1) ) * q_min_AI         &
                             + ( 1.0_RP-cmask(i,j,k,l,1) ) * CNST_MAX_REAL
          Qin4(i+1,j,    I_min) = (      cmask(i,j,k,l,1) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(i,j,k,l,1) ) * q_min_AI
          Qin1(i,j,    I_max) = (      cmask(i,j,k,l,1) ) * q_max_AI         &
                             + ( 1.0_RP-cmask(i,j,k,l,1) ) * (-CNST_MAX_REAL)
          Qin4(i+1,j,    I_max) = (      cmask(i,j,k,l,1) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(i,j,k,l,1) ) * q_max_AI

          Qin2(i,j,    I_min) = (      cmask(i,j,k,l,2) ) * q_min_AIJ        &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * CNST_MAX_REAL
          Qin5(i+1,j+1,    I_min) = (      cmask(i,j,k,l,2) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * q_min_AIJ
          Qin2(i,j,    I_max) = (      cmask(i,j,k,l,2) ) * q_max_AIJ        &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * (-CNST_MAX_REAL)
          Qin5(i+1,j+1,    I_max) = (      cmask(i,j,k,l,2) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * q_max_AIJ

          Qin3(i,j,    I_min) = (      cmask(i,j,k,l,3) ) * q_min_AJ         &
                             + ( 1.0_RP-cmask(i,j,k,l,3) ) * CNST_MAX_REAL
          Qin6(i,j+1,    I_min) = (      cmask(i,j,k,l,3) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(i,j,k,l,3) ) * q_min_AJ
          Qin3(i,j,    I_max) = (      cmask(i,j,k,l,3) ) * q_max_AJ         &
                             + ( 1.0_RP-cmask(i,j,k,l,3) ) * (-CNST_MAX_REAL)
          Qin6(i,j+1,    I_max) = (      cmask(i,j,k,l,3) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(i,j,k,l,3) ) * q_max_AJ
       enddo
       enddo


       if ( ADM_have_sgp(l) ) then
!$omp master

          j=ADM_gmin-1
          i=ADM_gmin-1

          q_min_AIJ = min( q(i,j,k,l), q(i+1,j+1,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l) )
          q_max_AIJ = max( q(i,j,k,l), q(i+1,j+1,k,l), q(i+1,j+1,k,l), q(i,j+1,k,l) )

          Qin2(i,j,    I_min) = (      cmask(i,j,k,l,2) ) * q_min_AIJ        &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * CNST_MAX_REAL
          Qin5(i+1,j+1,    I_min) = (      cmask(i,j,k,l,2) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * q_min_AIJ
          Qin2(i,j,    I_max) = (      cmask(i,j,k,l,2) ) * q_max_AIJ        &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * (-CNST_MAX_REAL)
          Qin5(i+1,j+1,    I_max) = (      cmask(i,j,k,l,2) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(i,j,k,l,2) ) * q_max_AIJ

!$omp end master
!$omp barrier
       endif


!ocl simd
!$omp do
      do j = ADM_gmin,ADM_gmax

#ifdef SINGLE_SIM
       do i = ADM_gmin-1,ADM_gmax/2
#else
       do i = ADM_gmin-1,ADM_gmax
#endif
       if( i == ADM_gmin-1 .and. j == ADM_gmin ) then
        cycle
       endif

!!!!OCL PREFETCH_READ(ch(1,i,j,k+2,l),level=2,strong=1)
!OCL PREFETCH_READ(q(i,j,k+2,l),level=2,strong=1)
!OCL PREFETCH_READ(d(i,j,k+2,l),level=2,strong=1)
!OCL PREFETCH_WRITE(Qout(i,j,k+2,l,I_min),level=2,strong=1)
!OCL PREFETCH_WRITE(Qout(i,j,k+2,l,I_max),level=2,strong=1)

       cmaskch1=min(0.0_RP,ch(i,j,k,l,1))
       cmaskch2=min(0.0_RP,ch(i,j,k,l,2))
       cmaskch3=min(0.0_RP,ch(i,j,k,l,3))
       cmaskch4=min(0.0_RP,ch(i,j,k,l,4))
       cmaskch5=min(0.0_RP,ch(i,j,k,l,5))
       cmaskch6=min(0.0_RP,ch(i,j,k,l,6))
          Cin_sum     =       cmaskch1 &
                      +       cmaskch2 &
                      +       cmaskch3 &
                      +       cmaskch4 &
                      +       cmaskch5 &
                      +       cmaskch6
          Cout_sum    =  ch(i,j,k,l,1)-cmaskch1 &
                      +  ch(i,j,k,l,2)-cmaskch2 &
                      +  ch(i,j,k,l,3)-cmaskch3 &
                      +  ch(i,j,k,l,4)-cmaskch4 &
                      +  ch(i,j,k,l,5)-cmaskch5 &
                      +  ch(i,j,k,l,6)-cmaskch6
          CQin_min_sum    = cmaskch1 * Qin1(i,j,I_min) &
                          + cmaskch2 * Qin2(i,j,I_min) &
                          + cmaskch3 * Qin3(i,j,I_min) &
                          + cmaskch4 * Qin4(i,j,I_min) &
                          + cmaskch5 * Qin5(i,j,I_min) &
                          + cmaskch6 * Qin6(i,j,I_min)

          CQin_max_sum    = cmaskch1 * Qin1(i,j,I_max) &
                          + cmaskch2 * Qin2(i,j,I_max) &
                          + cmaskch3 * Qin3(i,j,I_max) &
                          + cmaskch4 * Qin4(i,j,I_max) &
                          + cmaskch5 * Qin5(i,j,I_max) &
                          + cmaskch6 * Qin6(i,j,I_max)
          qnext_min    = min(Qin1(i,j,I_min),Qin2(i,j,I_min),Qin3(i,j,I_min),Qin4(i,j,I_min),&
                             Qin5(i,j,I_min),Qin6(i,j,I_min) )
          qnext_max    = max(Qin1(i,j,I_max),Qin2(i,j,I_max),Qin3(i,j,I_max),Qin4(i,j,I_max),&
                             Qin5(i,j,I_max),Qin6(i,j,I_max) )

          tmp=qnext_min
          qnext_min=q(i,j,k,l)
          if( tmp .ne.  CNST_MAX_REAL ) qnext_min = tmp
          tmp=qnext_max
          qnext_max=q(i,j,k,l)
          if( tmp .ne. -CNST_MAX_REAL ) qnext_max = tmp
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum)-CNST_EPS_ZERO) ! if Cout_sum = 0, sw = 1

          Qout(i,j,k,l,I_min) = ( q(i,j,k,l) - CQin_max_sum - qnext_max*(1.0_RP-Cin_sum-Cout_sum+d(i,j,k,l)) ) &
                        / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                        + q(i,j,k,l) * zerosw
          Qout(i,j,k,l,I_max) = ( q(i,j,k,l) - CQin_min_sum - qnext_min*(1.0_RP-Cin_sum-Cout_sum+d(i,j,k,l)) ) &
                        / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                        + q(i,j,k,l) * zerosw
       enddo ! i loop
       enddo ! j loop
!$omp end do nowait
!$omp do
       do j = ADM_gmin-1,ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin-1,ADM_gmax/2
#else
       do i = ADM_gmin-1,ADM_gmax
#endif
!!!!OCL PREFETCH_WRITE(q_a(1,i,j,k+2,l),level=2,strong=1)
          q_a(i,j,k,l,1) = (      cmask(i,j,k,l,1) ) * min(max(q_a(i,j,k,l,1), Qin1(i,j   ,I_min)), Qin1(i,j    ,I_max)) &
                       + ( 1.0_RP-cmask(i,j,k,l,1) ) * min(max(q_a(i,j,k,l,1), Qin4(i+1,j   ,I_min)), Qin4(i+1,j    ,I_max))

          q_a(i,j,k,l,2) = (      cmask(i,j,k,l,2) ) * min(max(q_a(i,j,k,l,2), Qin2(i,j,I_min)), Qin2(i,j,I_max)) &
                       + ( 1.0_RP-cmask(i,j,k,l,2) ) * min(max(q_a(i,j,k,l,2), Qin5(i+1,j+1,I_min)), Qin5(i+1,j+1,I_max))
          q_a(i,j,k,l,3) = (      cmask(i,j,k,l,3) ) * min(max(q_a(i,j,k,l,3), Qin3(i,j    ,I_min)), Qin3(i,j    ,I_max)) &
                       + ( 1.0_RP-cmask(i,j,k,l,3) ) * min(max(q_a(i,j,k,l,3), Qin6(i,j+1    ,I_min)), Qin6(i,j+1    ,I_max))

       enddo
       enddo
    enddo
!$omp end parallel
    enddo

!    call timer_end(3)
    !call DEBUG_rapend  ('loop2')

    !call DEBUG_rapstart('loop3')
!    call timer_sta(4)
    do l = 1, ADM_lall
!$omp parallel
!$omp do
    do k = 1, ADM_kall
       do j = ADM_gmin-1,ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin-1,ADM_gmax/2
#else
       do i = ADM_gmin-1,ADM_gmax
#endif

          q_a(i,j,k,l,1) = (      cmask(i,j,k,l,1) ) &
                       * max(min(q_a(i,j,k,l,1), Qout(  i+1,j  ,k,l,I_max)), Qout(  i+1,j  ,k,l,I_min)) &
                       + ( 1.0_RP-cmask(i,j,k,l,1) ) &
                       * max(min(q_a(i,j,k,l,1), Qout(  i,j    ,k,l,I_max)), Qout(  i,j    ,k,l,I_min))
          q_a(i+1,j,k,l,4) = q_a(i,j,k,l,1)

          q_a(i,j,k,l,2) = (      cmask(i,j,k,l,2) ) &
                       * max(min(q_a(i,j,k,l,2), Qout(  i+1,j+1,k,l,I_max)), Qout(  i+1,j+1,k,l,I_min)) &
                       + ( 1.0_RP-cmask(i,j,k,l,2) ) &
                       * max(min(q_a(i,j,k,l,2), Qout(  i,j    ,k,l,I_max)), Qout(  i,j    ,k,l,I_min))
          q_a(i+1,j+1,k,l,5) = q_a(i,j,k,l,2)

          q_a(i,j,k,l,3) = (      cmask(i,j,k,l,3) ) &
                       * max(min(q_a(i,j,k,l,3), Qout(  i,j+1  ,k,l,I_max)), Qout(  i,j+1  ,k,l,I_min)) &
                       + ( 1.0_RP-cmask(i,j,k,l,3) ) &
                       * max(min(q_a(i,j,k,l,3), Qout(  i,j    ,k,l,I_max)), Qout(  i,j    ,k,l,I_min))
          q_a(i,j+1,k,l,6) = q_a(i,j,k,l,3)
       enddo
       enddo
    enddo
!$omp end parallel
    enddo
!    call timer_end(4)
    !call DEBUG_rapend  ('loop3')

!    call DEBUG_rapend  ('____Horizontal_Adv_limiter')

    return
  end subroutine horizontal_limiter_thuburn

end module mod_src_tracer
