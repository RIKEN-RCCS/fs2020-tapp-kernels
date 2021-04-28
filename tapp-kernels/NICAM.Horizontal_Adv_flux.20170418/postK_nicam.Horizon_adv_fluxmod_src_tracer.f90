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

  !
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
  integer, public, parameter      :: ADM_lall_para = 1
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  !
  !------ Present Local region number ! 2010.4.26 M.Satoh
  integer, public, save      :: ADM_l_me

  logical, public, save :: ADM_have_sgp(ADM_lall_para) = .true. ! region have singlar point?

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
#else
  integer, public, save      :: ADM_gall_1d_in = 130
#endif
  integer, public, parameter :: ADM_gall_1d_para = 130
#ifdef SINGLE_SIM
  integer, public, parameter :: ADM_gall_1d_para_in = 65
#else
  integer, public, parameter :: ADM_gall_1d_para_in = 130
#endif
#ifdef SINGLE_SIM
  integer, public, save      :: ADM_gall = 16900 / 2
  integer, public, parameter :: ADM_gall_para = 16900 / 2
#else
  integer, public, save      :: ADM_gall = 16900 
  integer, public, parameter :: ADM_gall_para = 16900
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
  integer, public, save      :: ADM_kmax = 47
!  integer, public, save      :: ADM_kmax = 95
  integer, public, save      :: ADM_kall = 48
!  integer, public, save      :: ADM_kall = 96

end module mod_adm

module mod_grd
  use mod_precision
  use mod_adm, only: &
      ADM_gall,    &
      ADM_gall_para,    &
      ADM_gall_pl, &
      ADM_gall_1d,    &
      ADM_gall_1d_in,    &
      ADM_gall_1d_para_in,    &
      ADM_gall_1d_para,    &
      ADM_lall,    &
      ADM_lall_para,    &
      ADM_lall_pl, &
      ADM_kall
  use mod_adm, only: &
      AI  => ADM_AI,  &
      AJ  => ADM_AJ,  &
      K0  => ADM_KNONE
  implicit none

  !------ Indentifiers for the directions in the Cartesian coordinate.
  integer, public, parameter :: GRD_XDIR=1
  integer, public, parameter :: GRD_YDIR=2
  integer, public, parameter :: GRD_ZDIR=3

  !------ Grid points ( CELL ARC )
  real(RP), public, save :: GRD_xr   (ADM_gall_para ,K0, ADM_lall_para ,AI:AJ, GRD_XDIR:GRD_ZDIR)

  real(RP), public, save :: GRD_xr_ij (ADM_gall_para ,K0, ADM_lall_para,AI:AJ,GRD_XDIR:GRD_ZDIR)

end module mod_grd

module mod_gmtr
  use mod_precision
  use mod_adm, only: &
       ADM_TI,           &
       ADM_TJ,           &
       ADM_AI,           &
       ADM_AJ,           &
       ADM_gmin,         &
       ADM_gmax,         &
       ADM_gall,         &
       ADM_gall_para,         &
       ADM_gall_1d,      &
       ADM_gall_1d_para,      &
       ADM_gall_1d_in,      &
       ADM_gall_1d_para_in,      &
       ADM_gall_pl,      &
       ADM_lall,         &
       ADM_lall_para,         &
       ADM_lall_pl,      &
       K0 => ADM_KNONE
  implicit none

  integer, public, parameter :: GMTR_T_nmax_var = 7

  integer, public, parameter :: GMTR_T_AREA  = 1
  integer, public, parameter :: GMTR_T_RAREA = 2
  integer, public, parameter :: GMTR_T_W1    = 3
  integer, public, parameter :: GMTR_T_W2    = 4
  integer, public, parameter :: GMTR_T_W3    = 5
  integer, public, parameter :: GMTR_T_LAT   = 6
  integer, public, parameter :: GMTR_T_LON   = 7

  integer, public, parameter :: GMTR_P_nmax_var = 10

  integer, public, parameter :: GMTR_P_AREA  = 1
  integer, public, parameter :: GMTR_P_RAREA = 2
  integer, public, parameter :: GMTR_P_IX    = 3
  integer, public, parameter :: GMTR_P_IY    = 4
  integer, public, parameter :: GMTR_P_IZ    = 5
  integer, public, parameter :: GMTR_P_JX    = 6
  integer, public, parameter :: GMTR_P_JY    = 7
  integer, public, parameter :: GMTR_P_JZ    = 8
  integer, public, parameter :: GMTR_P_LAT   = 9
  integer, public, parameter :: GMTR_P_LON   = 10

  integer, public, parameter :: GMTR_A_nmax_var    = 12
  integer, public, parameter :: GMTR_A_nmax_var_pl = 18

  integer, public, parameter :: GMTR_A_HNX  = 1
  integer, public, parameter :: GMTR_A_HNY  = 2
  integer, public, parameter :: GMTR_A_HNZ  = 3
  integer, public, parameter :: GMTR_A_HTX  = 4
  integer, public, parameter :: GMTR_A_HTY  = 5
  integer, public, parameter :: GMTR_A_HTZ  = 6
  integer, public, parameter :: GMTR_A_TNX  = 7
  integer, public, parameter :: GMTR_A_TNY  = 8
  integer, public, parameter :: GMTR_A_TNZ  = 9
  integer, public, parameter :: GMTR_A_TTX  = 10
  integer, public, parameter :: GMTR_A_TTY  = 11
  integer, public, parameter :: GMTR_A_TTZ  = 12

  integer, public, parameter :: GMTR_A_TN2X = 13
  integer, public, parameter :: GMTR_A_TN2Y = 14
  integer, public, parameter :: GMTR_A_TN2Z = 15
  integer, public, parameter :: GMTR_A_TT2X = 16
  integer, public, parameter :: GMTR_A_TT2Y = 17
  integer, public, parameter :: GMTR_A_TT2Z = 18

  real(RP), public, save :: GMTR_T_var(ADM_gall_para,K0,ADM_lall_para,ADM_TI:ADM_TJ,GMTR_T_nmax_var) = 0.D0

end module mod_gmtr

module mod_cnst
  use mod_precision
  implicit none
  private
  !++ Public procedure

  !++ Public parameters & variables
  real(RP), public, save :: CNST_ERADIUS = 6.37122E6_RP ! Radius of the Earth[m]
  real(RP), public, save :: CNST_EOHM    = 7.292E-5_RP   ! Angular velocity of the Earth [/s]
  real(RP), public, save :: CNST_EGRAV   = 9.80616_RP  ! Gravitational accerlaration of the Earth [m/s2]

  real(RP), public, save :: CNST_RAIR    =  287.0_RP   ! Gas constant of air
  real(RP), public, save :: CNST_RVAP    =  461.5_RP   ! Gas constant of vapor

  real(RP), public, save :: CNST_CP      = 1004.5_RP   ! Specific heat of air(consant pressure)
  real(RP), public, save :: CNST_CV                   ! Specific heat of air(consant volume)

  real(RP), public, save :: CNST_CPV     = 1846.0_RP   ! Specific heat of vapor(consant pressure)
  real(RP), public, save :: CNST_CVV                  ! Specific heat of vapor(consant volume)
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
  real(RP), public, save :: CNST_LH0   = 2.501E+6_RP
  !
  !------ Latent heat of sublimation at 0C
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
  !
  !====== Misc. constants ======
  !
  !------ Definition of PI
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
!-------------------------------------------------------------------------------
!>
!! Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators.
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
!-------------------------------------------------------------------------------

module mod_oprt
  use mod_precision
  use mod_adm, only: &
    ADM_have_pl,    &   !cx logical
    ADM_have_sgp,   &   !cx logical array
    ADM_lall,       &
    ADM_lall_para,       &
    ADM_lall_pl,    &
    ADM_gall,       &
    ADM_gall_para,       &
    ADM_gall_pl,    &
    ADM_kall,       &
    ADM_gall_1d,    &
    ADM_gall_1d_para,    &
    ADM_gall_1d_in,    &
    ADM_gall_1d_para_in,    &
    ADM_gmin,       &
    ADM_gmax,       &
    ADM_gslf_pl,    &   !cx parameter
    ADM_gmin_pl,    &   !cx parameter
    ADM_gmax_pl,    &
    TI  => ADM_TI,  &
    TJ  => ADM_TJ,  &
    AI  => ADM_AI,  &
    AIJ => ADM_AIJ, &
    AJ  => ADM_AJ,  &
    K0  => ADM_KNONE

  use mod_adm, only: my_id => ADM_prc_me

  use mod_gmtr, only: &
     GMTR_T_var, &
     GMTR_T_nmax_var, &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z

  implicit none

  !++ Public parameters & variables
  integer, public, save :: OPRT_nstart
  integer, public, save :: OPRT_nend

  ! < for diffusion operator >
  real(RP), public, save :: cinterp_TN (AI:AJ,1:3,ADM_gall_para,ADM_lall_para)
  real(RP), public, save :: cinterp_HN (AI:AJ,1:3,ADM_gall_para,ADM_lall_para)
  real(RP), public, save :: cinterp_TRA(TI:TJ,ADM_gall_para,ADM_lall_para)
  real(RP), public, save :: cinterp_PRA(ADM_gall_para,ADM_lall_para)

  real(RP), public, save :: cinterp_HN_ij(ADM_gall_para,ADM_lall_para,AI:AJ,1:3)
  real(RP), public, save :: cinterp_PRA_ij(ADM_gall_para,ADM_lall_para)

contains

  subroutine OPRT_snap_read ( scl, scl_pl, kh, kh_pl, mfact )
    implicit none
    integer :: i

    real(RP), intent(out)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: mfact
    integer :: j,l,n,m,g

!cx Better To Do: first touch the following arrays

        cinterp_TN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_HN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_TRA(TI:TJ,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_PRA(      1:ADM_gall,1:ADM_lall) = 1.0
        scl    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        scl_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        kh     (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        kh_pl  (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        mfact = 1.0
        return

  end subroutine OPRT_snap_read

end module mod_oprt



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
module mod_src_tracer
  !-----------------------------------------------------------------------------

  use mod_precision
  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE
  use mod_grd, only: &
     XDIR => GRD_XDIR, &
     YDIR => GRD_YDIR, &
     ZDIR => GRD_ZDIR
  use mod_gmtr, only: &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z

    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gall_1d_in,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax
    use mod_grd, only: &
       GRD_xr_ij,   &
       GRD_xr   
    use mod_gmtr, only: &
       GMTR_T_nmax_var, &
       GMTR_T_var
    use mod_oprt, only: &
       cinterp_HN_ij,  &
       cinterp_PRA_ij, &
!       cinterp_HN,  &
       cinterp_PRA
    use mod_cnst, only: &
       CNST_EPS_ZERO
    use mod_adm, only: ADM_prc_me

  implicit none
  private

  public :: horizontal_flux
  public :: H_Adv_flux_snap_read


contains


  !-----------------------------------------------------------------------------
  !> prepare horizontal advection trem: mass flux, GRD_xc
  subroutine horizontal_flux( &
       flx_h,  flx_h_pl,  &
       GRD_xc, GRD_xc_pl, &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt,GRD_xr_ij_l,cinterp_HN_ij_l )

    implicit none

    real(RP), intent(out) :: flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(RP), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: dt
    real(RP), intent(out) :: flx_h    (ADM_gall ,ADM_kall,ADM_lall,6)               ! horizontal mass flux
    real(RP), intent(out) :: GRD_xc   (ADM_gall ,ADM_kall,ADM_lall,AI:AJ,XDIR:ZDIR )! mass centroid position
    real(RP), intent(in)  :: rho      (ADM_gall ,ADM_kall,ADM_lall   )              ! rho at cell center
    real(RP), intent(in)  :: rhovx    (ADM_gall ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovy    (ADM_gall ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovz    (ADM_gall ,ADM_kall,ADM_lall   )
    real(RP)  :: GRD_xr_ij_l    (ADM_gall ,K0, ADM_lall,AI:AJ,XDIR:ZDIR)
    real(RP)  :: cinterp_HN_ij_l(ADM_gall,ADM_lall,AI:AJ,1:3)

    real(RP) :: tmp_rhovxt
    real(RP) :: tmp_rhovyt
    real(RP) :: tmp_rhovzt
    real(RP) :: flux
    real(RP) :: rrhoa2
    integer :: i,j,ij,iijj, k, l
    integer :: nstart,nend
    integer :: n
    integer :: suf

    real(RP) :: flux1,flux2,flux3,rrhoa21,rrhoa22,rrhoa23
!    integer(kind=4),parameter :: block_size=4224 !  4block
!    integer(kind=4),parameter :: block_size=2112 !  8block
!    integer(kind=4),parameter :: block_size=1056 ! 16block
!    integer(kind=4),parameter :: block_size= 528 ! 32block
    integer(kind=4),parameter :: block_size= 128 ! SIMD noswp,nounroll,striping=2/striping=4
!    integer(kind=4),parameter :: block_size=  64 ! SIMD noswp,nounroll,striping=2/striping=4
!    integer(kind=4),parameter :: block_size=  32 ! SIMD noswp,nounroll,striping=2
!    integer(kind=4),parameter :: block_size=  16 ! SIMD noswp,nounroll
    real(RP) :: tmp_rhovxt1(block_size), tmp_rhovyt1(block_size), tmp_rhovzt1(block_size)
    real(RP) :: tmp_rhovxt2(block_size), tmp_rhovyt2(block_size), tmp_rhovzt2(block_size)
    real(RP) :: tmp_rhovxt3(block_size), tmp_rhovyt3(block_size), tmp_rhovzt3(block_size)
    !      local_array           xxx:max(block_size,ADM_gall_1d) -> block_size+ADM_gall_1d
    real(RP) :: rhot1  (           0:block_size+ADM_gall_1d) ! rho at cell vertex
    real(RP) :: rhovxt1(           0:block_size+ADM_gall_1d)
    real(RP) :: rhovyt1(           0:block_size+ADM_gall_1d)
    real(RP) :: rhovzt1(           0:block_size+ADM_gall_1d)
    real(RP) :: rhot2  (-ADM_gall_1d:block_size+ADM_gall_1d) ! rho at cell vertex
    real(RP) :: rhovxt2(-ADM_gall_1d:block_size+ADM_gall_1d)
    real(RP) :: rhovyt2(-ADM_gall_1d:block_size+ADM_gall_1d)
    real(RP) :: rhovzt2(-ADM_gall_1d:block_size+ADM_gall_1d)

    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------
!cx    if (mark_for_snapshots(ip_Horizontal_Adv_flux).ne.0) then
!cx    call snapshot_seq_bin_open ("Horizontal_Adv_flux", ADM_prc_me)
!cx    call H_Adv_flux_snap_write ( &
!cx       rho,    rho_pl,    &
!cx       rhovx,  rhovx_pl,  &
!cx       rhovy,  rhovy_pl,  &
!cx       rhovz,  rhovz_pl,  &
!cx       dt                 )
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Horizontal_Adv_flux) = 0
!cx    endif
!    call DEBUG_rapstart('____Horizontal_Adv_flux')
!flux2-peeling-check
       if(ADM_gall_1d.le.2)        stop !fj error-check
       if(ADM_gmin   .lt.2)        stop !fj error-check
       if(ADM_gall_1d.le.ADM_gmin) stop !fj error-check
!$omp parallel default(none), &
!$omp private(nstart,nend,i,j,ij,iijj,k,l, &
!$omp         flux1,flux2,flux3,rrhoa21,rrhoa22,rrhoa23, &
!$omp         tmp_rhovxt1,tmp_rhovyt1,tmp_rhovzt1, &
!$omp         tmp_rhovxt2,tmp_rhovyt2,tmp_rhovzt2, &
!$omp         tmp_rhovxt3,tmp_rhovyt3,tmp_rhovzt3, &
!$omp         rhot1, rhovxt1,rhovyt1,rhovzt1,rhot2 ,rhovxt2,rhovyt2,rhovzt2, &
!$omp         tmp_rhovxt,tmp_rhovyt,tmp_rhovzt, flux, rrhoa2) &
!$omp shared(ADM_gmin,ADM_gmax,ADM_kall,ADM_lall,ADM_gall_1d, ADM_gall_1d_in,&
!$omp        flx_h,dt, GRD_xr_ij_l,GRD_xc, ADM_have_sgp, &
!$omp        CNST_EPS_ZERO,cinterp_PRA_ij,cinterp_HN_ij_l,GMTR_T_var, &
!$omp        rhovz,rhovy,rhovx,rho)
    do l = 1, ADM_lall
!$omp do
    do k = 1, ADM_kall
!-------------------------------------------------------------------------------------------
!!! flux1-peeling-loop
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       do ij = nstart, ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
          rhot1(ij)     = rho  (ij    ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                        + rho  (ij+1  ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                        + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhot2(ij)     = rho  (ij    ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                        + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                        + rho  (ij+ADM_gall_1d,k,l)   * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = nstart, ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
          rhovxt1(ij)   = rhovx(ij    ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                        + rhovx(ij+1  ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                        + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovxt2(ij)   = rhovx(ij    ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                        + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                        + rhovx(ij+ADM_gall_1d,k,l)   * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = nstart, ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
          rhovyt1(ij)   = rhovy(ij    ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                        + rhovy(ij+1  ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) & 
                        + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovyt2(ij)   = rhovy(ij    ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                        + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                        + rhovy(ij+ADM_gall_1d,k,l)   * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = nstart, ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
          rhovzt1(ij)   = rhovz(ij    ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                        + rhovz(ij+1  ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                        + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovzt2(ij)   = rhovz(ij    ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                        + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                        + rhovz(ij+ADM_gall_1d,k,l)   * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!-------------------------------------------------------------------------------------------
       if ( ADM_have_sgp(l) ) then
          rhot1  (suf(ADM_gmin-1,ADM_gmin-1)) = rhot2  (suf(ADM_gmin,ADM_gmin-1))
          rhovxt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovxt2(suf(ADM_gmin,ADM_gmin-1))
          rhovyt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovyt2(suf(ADM_gmin,ADM_gmin-1))
          rhovzt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovzt2(suf(ADM_gmin,ADM_gmin-1))
       endif
!-------------------------------------------------------------------------------------------
!!! flux2-peeling-loop
       do ij = ADM_gall_1d * (ADM_gmin-2) + ADM_gmin-1, ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2

          tmp_rhovxt=rhovxt1(ij)+rhovxt2(ij)
          tmp_rhovyt=rhovyt1(ij)+rhovyt2(ij)
          tmp_rhovzt=rhovzt1(ij)+rhovzt2(ij)

          flux = 0.5_RP * ( tmp_rhovxt * cinterp_HN_ij_l(ij,l,AIJ,1) &
                           + tmp_rhovyt * cinterp_HN_ij_l(ij,l,AIJ,2) &
                           + tmp_rhovzt * cinterp_HN_ij_l(ij,l,AIJ,3) )
          flx_h(ij    ,k,l,2) =  flux * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+1+ADM_gall_1d,k,l,5) = -flux * cinterp_PRA_ij(ij+1+ADM_gall_1d,l) * dt

          rrhoa2 = 1.0_RP / max( rhot1(ij) + rhot2(ij), CNST_EPS_ZERO ) ! doubled

          GRD_xc(ij,k,l,AIJ,XDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,XDIR) - tmp_rhovxt * rrhoa2 * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,YDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,YDIR) - tmp_rhovyt * rrhoa2 * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,ZDIR) - tmp_rhovzt * rrhoa2 * dt * 0.5_RP
       enddo
       do ij = ADM_gall_1d * (ADM_gmin-2) + ADM_gmin  , ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
          tmp_rhovxt=rhovxt2(ij)+rhovxt1(ij-1)
          tmp_rhovyt=rhovyt2(ij)+rhovyt1(ij-1)
          tmp_rhovzt=rhovzt2(ij)+rhovzt1(ij-1)

          flux = 0.5_RP * ( tmp_rhovxt * cinterp_HN_ij_l(ij,l,AJ ,1) &
                          + tmp_rhovyt * cinterp_HN_ij_l(ij,l,AJ ,2) &
                          + tmp_rhovzt * cinterp_HN_ij_l(ij,l,AJ ,3) )

          flx_h(ij  ,k,l,3) =  flux * cinterp_PRA_ij(ij  ,l) * dt
          flx_h(ij+ADM_gall_1d,k,l,6) = -flux * cinterp_PRA_ij(ij+ADM_gall_1d,l) * dt

          rrhoa2 = 1.0_RP / max( rhot2(ij) + rhot1(ij-1), CNST_EPS_ZERO ) ! doubled

          GRD_xc(ij,k,l,AJ,XDIR) = GRD_xr_ij_l(ij,K0,l,AJ,XDIR) - tmp_rhovxt * rrhoa2 * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ,YDIR) = GRD_xr_ij_l(ij,K0,l,AJ,YDIR) - tmp_rhovyt * rrhoa2 * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AJ,ZDIR) - tmp_rhovzt * rrhoa2 * dt * 0.5_RP
       enddo
!-------------------------------------------------------------------------------------------
!!! blocking-start
      rhot1  (0) = rhot1  (ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2)
      rhovxt1(0) = rhovxt1(ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2)
      rhovyt1(0) = rhovyt1(ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2)
      rhovzt1(0) = rhovzt1(ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2)
      do i=1, ADM_gall_1d
        nend=ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-2
        rhot2  (i-ADM_gall_1d) = rhot2  (i+nend-ADM_gall_1d)
        rhovxt2(i-ADM_gall_1d) = rhovxt2(i+nend-ADM_gall_1d)
        rhovyt2(i-ADM_gall_1d) = rhovyt2(i+nend-ADM_gall_1d)
        rhovzt2(i-ADM_gall_1d) = rhovzt2(i+nend-ADM_gall_1d)
      end do
!-------------------------------------------------------------------------------------------
      nend   = suf(ADM_gmax  ,ADM_gmax  )
      do iijj = ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-1, nend , block_size
      if(block_size .le. (nend - iijj + 1) )then
!!! flux1-loop-body
#include "horizontal_flux1_L1_prefetch02.h"
#include "horizontal_flux1_L2_prefetch01.h"
!       do ij = ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-1, nend
       do ij = iijj, iijj+block_size-1
          rhot1(ij-iijj+1) = rho  (ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                           + rho  (ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                           + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhot2(ij-iijj+1) = rho  (ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                           + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                           + rho  (ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
#include "horizontal_flux1_L1_prefetch03.h"
#include "horizontal_flux1_L2_prefetch02.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rhovxt1(ij-iijj+1) = rhovx(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovx(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                             + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovxt2(ij-iijj+1) = rhovx(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovx(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
#include "horizontal_flux1_L1_prefetch04.h"
#include "horizontal_flux1_L2_prefetch03.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rhovyt1(ij-iijj+1) = rhovy(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovy(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) & 
                             + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovyt2(ij-iijj+1) = rhovy(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovy(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
#include "horizontal_flux2_L1_prefetch01.h"
#include "horizontal_flux1_L2_prefetch04.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rhovzt1(ij-iijj+1) = rhovz(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovz(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                             + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovzt2(ij-iijj+1) = rhovz(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovz(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!-------------------------------------------------------------------------------------------
!!! flux2-loop-body
       !--- calculate flux and mass centroid position
       nend   = suf(ADM_gmax  ,ADM_gmax  )
!       do ij = ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-1, nend
#include "horizontal_flux2_L1_prefetch02.h"
#include "horizontal_flux2_L2_prefetch01.h"
       do ij = iijj, iijj+block_size-1
          tmp_rhovxt1(ij-iijj+1)=rhovxt2(ij-iijj+1-ADM_gall_1d)+rhovxt1(ij-iijj+1)
          tmp_rhovxt2(ij-iijj+1)=rhovxt1(ij-iijj+1)            +rhovxt2(ij-iijj+1)
          tmp_rhovxt3(ij-iijj+1)=rhovxt2(ij-iijj+1)            +rhovxt1(ij-iijj+1-1)
       enddo
#include "horizontal_flux2_L1_prefetch03.h"
#include "horizontal_flux2_L2_prefetch02.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          tmp_rhovyt1(ij-iijj+1)=rhovyt2(ij-iijj+1-ADM_gall_1d)+rhovyt1(ij-iijj+1)
          tmp_rhovyt2(ij-iijj+1)=rhovyt1(ij-iijj+1)            +rhovyt2(ij-iijj+1)
          tmp_rhovyt3(ij-iijj+1)=rhovyt2(ij-iijj+1)            +rhovyt1(ij-iijj+1-1)
       enddo
#include "horizontal_flux2_L1_prefetch04.h"
#include "horizontal_flux2_L2_prefetch03.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          tmp_rhovzt1(ij-iijj+1)=rhovzt2(ij-iijj+1-ADM_gall_1d)+rhovzt1(ij-iijj+1)
          tmp_rhovzt2(ij-iijj+1)=rhovzt1(ij-iijj+1)            +rhovzt2(ij-iijj+1)
          tmp_rhovzt3(ij-iijj+1)=rhovzt2(ij-iijj+1)            +rhovzt1(ij-iijj+1-1)
       enddo
#include "horizontal_flux2_L1_prefetch05.h"
#include "horizontal_flux2_L2_prefetch04.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          flux1 = 0.5_RP * ( tmp_rhovxt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,1) &
                           + tmp_rhovyt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,2) &
                           + tmp_rhovzt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,3) )
          flx_h(ij    ,k,l,1) =  flux1 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+1  ,k,l,4) = -flux1 * cinterp_PRA_ij(ij+1  ,l) * dt
       enddo
#include "horizontal_flux2_L1_prefetch06.h"
#include "horizontal_flux2_L2_prefetch05.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          flux2 = 0.5_RP * ( tmp_rhovxt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,1) &
                           + tmp_rhovyt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,2) &
                           + tmp_rhovzt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,3) )
          flx_h(ij    ,k,l,2) =  flux2 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+1+ADM_gall_1d,k,l,5) = -flux2 * cinterp_PRA_ij(ij+1+ADM_gall_1d,l) * dt
       enddo
#include "horizontal_flux2_L1_prefetch07.h"
#include "horizontal_flux2_L2_prefetch06.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          flux3 = 0.5_RP * ( tmp_rhovxt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,1) &
                           + tmp_rhovyt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,2) &
                           + tmp_rhovzt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,3) )
          flx_h(ij    ,k,l,3) =  flux3 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+ADM_gall_1d  ,k,l,6) = -flux3 * cinterp_PRA_ij(ij+ADM_gall_1d  ,l) * dt
       enddo
#include "horizontal_flux2_L1_prefetch08.h"
#include "horizontal_flux2_L2_prefetch07.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rrhoa21 = 1.0_RP / max( rhot2(ij-iijj+1-ADM_gall_1d) + rhot1(ij-iijj+1)  , CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt1(ij-iijj+1) = tmp_rhovxt1(ij-iijj+1) * rrhoa21
          tmp_rhovyt1(ij-iijj+1) = tmp_rhovyt1(ij-iijj+1) * rrhoa21
          tmp_rhovzt1(ij-iijj+1) = tmp_rhovzt1(ij-iijj+1) * rrhoa21
       enddo
#include "horizontal_flux2_L1_prefetch09.h"
#include "horizontal_flux2_L2_prefetch08.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rrhoa22 = 1.0_RP / max( rhot1(ij-iijj+1)             + rhot2(ij-iijj+1)  , CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt2(ij-iijj+1) = tmp_rhovxt2(ij-iijj+1) * rrhoa22
          tmp_rhovyt2(ij-iijj+1) = tmp_rhovyt2(ij-iijj+1) * rrhoa22
          tmp_rhovzt2(ij-iijj+1) = tmp_rhovzt2(ij-iijj+1) * rrhoa22
       enddo
#include "horizontal_flux2_L1_prefetch10.h"
#include "horizontal_flux2_L2_prefetch09.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          rrhoa23 = 1.0_RP / max( rhot2(ij-iijj+1)             + rhot1(ij-iijj+1-1), CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt3(ij-iijj+1) = tmp_rhovxt3(ij-iijj+1) * rrhoa23
          tmp_rhovyt3(ij-iijj+1) = tmp_rhovyt3(ij-iijj+1) * rrhoa23
          tmp_rhovzt3(ij-iijj+1) = tmp_rhovzt3(ij-iijj+1) * rrhoa23
       enddo
       ! copy: -129 <- -1
       rhot2  (-ADM_gall_1d+1) = rhot2  (-1)
       rhovxt2(-ADM_gall_1d+1) = rhovxt2(-1)
       rhovyt2(-ADM_gall_1d+1) = rhovyt2(-1)
       rhovzt2(-ADM_gall_1d+1) = rhovzt2(-1)
       ! copy: -128 <- 0
       rhot2  (-ADM_gall_1d+2) = rhot2  ( 0)
       rhovxt2(-ADM_gall_1d+2) = rhovxt2( 0)
       rhovyt2(-ADM_gall_1d+2) = rhovyt2( 0)
       rhovzt2(-ADM_gall_1d+2) = rhovzt2( 0)
#include "horizontal_flux2_L1_prefetch11.h"
#include "horizontal_flux2_L2_prefetch10.h"
!ocl norecurrence
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          GRD_xc(ij,k,l,AI ,XDIR) = GRD_xr_ij_l(ij,K0,l,AI ,XDIR) - tmp_rhovxt1(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AI ,YDIR) = GRD_xr_ij_l(ij,K0,l,AI ,YDIR) - tmp_rhovyt1(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AI ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AI ,ZDIR) - tmp_rhovzt1(ij-iijj+1) * dt * 0.5_RP
          rhot2  (ij-iijj+1-block_size)=rhot2  (ij-iijj+1)
          rhovxt2(ij-iijj+1-block_size)=rhovxt2(ij-iijj+1)
       enddo
#include "horizontal_flux2_L1_prefetch12.h"
#include "horizontal_flux2_L2_prefetch11.h"
!ocl norecurrence
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          GRD_xc(ij,k,l,AIJ,XDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,XDIR) - tmp_rhovxt2(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,YDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,YDIR) - tmp_rhovyt2(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,ZDIR) - tmp_rhovzt2(ij-iijj+1) * dt * 0.5_RP
          rhovyt2(ij-iijj+1-block_size)=rhovyt2(ij-iijj+1)
          rhovzt2(ij-iijj+1-block_size)=rhovzt2(ij-iijj+1)
       enddo
#include "horizontal_flux1_L1_prefetch01.h"
#include "horizontal_flux2_L2_prefetch12.h"
!ocl loop_nofusion
       do ij = iijj, iijj+block_size-1
          GRD_xc(ij,k,l,AJ ,XDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,XDIR) - tmp_rhovxt3(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ ,YDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,YDIR) - tmp_rhovyt3(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,ZDIR) - tmp_rhovzt3(ij-iijj+1) * dt * 0.5_RP
       enddo
       ! copy:
       rhot1  (0) = rhot1  (block_size-1 )
       rhovxt1(0) = rhovxt1(block_size-1 )
       rhovyt1(0) = rhovyt1(block_size-1 )
       rhovzt1(0) = rhovzt1(block_size-1 )
!!! blocking-end
!!! -------------------------------------------------------------------------------------------------------------------
      else 
!!! blocking-mod-start
       do ij = iijj, nend
          rhot1(ij-iijj+1) = rho  (ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                           + rho  (ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                           + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhot2(ij-iijj+1) = rho  (ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                           + rho  (ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                           + rho  (ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rhovxt1(ij-iijj+1) = rhovx(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovx(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                             + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovxt2(ij-iijj+1) = rhovx(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovx(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovx(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rhovyt1(ij-iijj+1) = rhovy(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovy(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) & 
                             + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovyt2(ij-iijj+1) = rhovy(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovy(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovy(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rhovzt1(ij-iijj+1) = rhovz(ij              ,k,l) * GMTR_T_var(ij,K0,l,TI,W1) &
                             + rhovz(ij+1            ,k,l) * GMTR_T_var(ij,K0,l,TI,W2) &
                             + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TI,W3)
          rhovzt2(ij-iijj+1) = rhovz(ij              ,k,l) * GMTR_T_var(ij,K0,l,TJ,W1) &
                             + rhovz(ij+1+ADM_gall_1d,k,l) * GMTR_T_var(ij,K0,l,TJ,W2) &
                             + rhovz(ij+ADM_gall_1d  ,k,l) * GMTR_T_var(ij,K0,l,TJ,W3)
       enddo
!-------------------------------------------------------------------------------------------
       !--- calculate flux and mass centroid position
       nend   = suf(ADM_gmax  ,ADM_gmax  )
!       do ij = ADM_gall_1d * (ADM_gmin-1) + ADM_gmin-1, nend
       do ij = iijj, nend
          tmp_rhovxt1(ij-iijj+1)=rhovxt2(ij-iijj+1-ADM_gall_1d)+rhovxt1(ij-iijj+1)
          tmp_rhovxt2(ij-iijj+1)=rhovxt1(ij-iijj+1)            +rhovxt2(ij-iijj+1)
          tmp_rhovxt3(ij-iijj+1)=rhovxt2(ij-iijj+1)            +rhovxt1(ij-iijj+1-1)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          tmp_rhovyt1(ij-iijj+1)=rhovyt2(ij-iijj+1-ADM_gall_1d)+rhovyt1(ij-iijj+1)
          tmp_rhovyt2(ij-iijj+1)=rhovyt1(ij-iijj+1       )     +rhovyt2(ij-iijj+1)
          tmp_rhovyt3(ij-iijj+1)=rhovyt2(ij-iijj+1)            +rhovyt1(ij-iijj+1-1)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          tmp_rhovzt1(ij-iijj+1)=rhovzt2(ij-iijj+1-ADM_gall_1d)+rhovzt1(ij-iijj+1)
          tmp_rhovzt2(ij-iijj+1)=rhovzt1(ij-iijj+1)            +rhovzt2(ij-iijj+1)
          tmp_rhovzt3(ij-iijj+1)=rhovzt2(ij-iijj+1)            +rhovzt1(ij-iijj+1-1)
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          flux1 = 0.5_RP * ( tmp_rhovxt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,1) &
                           + tmp_rhovyt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,2) &
                           + tmp_rhovzt1(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AI ,3) )
          flx_h(ij    ,k,l,1) =  flux1 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+1  ,k,l,4) = -flux1 * cinterp_PRA_ij(ij+1  ,l) * dt
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          flux2 = 0.5_RP * ( tmp_rhovxt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,1) &
                           + tmp_rhovyt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,2) &
                           + tmp_rhovzt2(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AIJ,3) )
          flx_h(ij    ,k,l,2) =  flux2 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+1+ADM_gall_1d,k,l,5) = -flux2 * cinterp_PRA_ij(ij+1+ADM_gall_1d,l) * dt
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          flux3 = 0.5_RP * ( tmp_rhovxt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,1) &
                           + tmp_rhovyt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,2) &
                           + tmp_rhovzt3(ij-iijj+1) * cinterp_HN_ij_l(ij,l,AJ ,3) )
          flx_h(ij    ,k,l,3) =  flux3 * cinterp_PRA_ij(ij    ,l) * dt
          flx_h(ij+ADM_gall_1d  ,k,l,6) = -flux3 * cinterp_PRA_ij(ij+ADM_gall_1d  ,l) * dt
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rrhoa21 = 1.0_RP / max( rhot2(ij-iijj+1-ADM_gall_1d) + rhot1(ij-iijj+1)  , CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt1(ij-iijj+1) = tmp_rhovxt1(ij-iijj+1) * rrhoa21
          tmp_rhovyt1(ij-iijj+1) = tmp_rhovyt1(ij-iijj+1) * rrhoa21
          tmp_rhovzt1(ij-iijj+1) = tmp_rhovzt1(ij-iijj+1) * rrhoa21
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rrhoa22 = 1.0_RP / max( rhot1(ij-iijj+1)             + rhot2(ij-iijj+1)  , CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt2(ij-iijj+1) = tmp_rhovxt2(ij-iijj+1) * rrhoa22
          tmp_rhovyt2(ij-iijj+1) = tmp_rhovyt2(ij-iijj+1) * rrhoa22
          tmp_rhovzt2(ij-iijj+1) = tmp_rhovzt2(ij-iijj+1) * rrhoa22
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          rrhoa23 = 1.0_RP / max( rhot2(ij-iijj+1)             + rhot1(ij-iijj+1-1), CNST_EPS_ZERO ) ! doubled
          tmp_rhovxt3(ij-iijj+1) = tmp_rhovxt3(ij-iijj+1) * rrhoa23
          tmp_rhovyt3(ij-iijj+1) = tmp_rhovyt3(ij-iijj+1) * rrhoa23
          tmp_rhovzt3(ij-iijj+1) = tmp_rhovzt3(ij-iijj+1) * rrhoa23
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          GRD_xc(ij,k,l,AI ,XDIR) = GRD_xr_ij_l(ij,K0,l,AI ,XDIR) - tmp_rhovxt1(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AI ,YDIR) = GRD_xr_ij_l(ij,K0,l,AI ,YDIR) - tmp_rhovyt1(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AI ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AI ,ZDIR) - tmp_rhovzt1(ij-iijj+1) * dt * 0.5_RP
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          GRD_xc(ij,k,l,AIJ,XDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,XDIR) - tmp_rhovxt2(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,YDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,YDIR) - tmp_rhovyt2(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AIJ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AIJ,ZDIR) - tmp_rhovzt2(ij-iijj+1) * dt * 0.5_RP
       enddo
!ocl loop_nofusion
       do ij = iijj, nend
          GRD_xc(ij,k,l,AJ ,XDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,XDIR) - tmp_rhovxt3(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ ,YDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,YDIR) - tmp_rhovyt3(ij-iijj+1) * dt * 0.5_RP
          GRD_xc(ij,k,l,AJ ,ZDIR) = GRD_xr_ij_l(ij,K0,l,AJ ,ZDIR) - tmp_rhovzt3(ij-iijj+1) * dt * 0.5_RP
       enddo
!!! blocking-mod-end
      endif
      enddo
!-------------------------------------------------------------------------------------------
       if ( ADM_have_sgp(l) ) then
          flx_h(suf(ADM_gmin,ADM_gmin),k,l,6) = 0.0_RP
       endif
!-------------------------------------------------------------------------------------------
    enddo
!$omp end do
    enddo
!$omp end parallel

!    call DEBUG_rapend  ('____Horizontal_Adv_flux')
    return
  end subroutine horizontal_flux


  subroutine H_Adv_flux_snap_read ( &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

    implicit none
    real(RP), intent(inout)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: dt
    integer i
    integer j,g,l,m,n,k,t

!cx Better To Do: first touch the following arrays

#ifndef OUTPUTCHECK
        rho      (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rho_pl   (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovx    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovx_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovy    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovy_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovz    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovz_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        dt = 1.0
        GRD_xr(:,:,:,:,:)=1.0
!        cinterp_HN(:,:,:,:)=1.0
        cinterp_PRA(:,:)=1.0
#else
        dt = 1.0
        do l = 1 , ADM_lall
        do k = 1 , ADM_kall
        do g = 1 , ADM_gall
          rho  (g,k,l) =  1.0 * g * k * l / 1000 * 1
          rhovx(g,k,l) =  1.0 * g * k * l / 1000 * 2
          rhovy(g,k,l) =  1.0 * g * k * l / 1000 * 3
          rhovz(g,k,l) =  1.0 * g * k * l / 1000 * 4
        enddo
        enddo
        enddo
        do l = 1, ADM_lall
          do j = 1, ADM_gall_1d
            do i = 1, ADM_gall_1d_in
              g = i + ( (j - 1)  * ADM_gall_1d_in)
              do m = 1,3
                do n = AI, AJ
                  cinterp_HN_ij(g,l,n,m) = 0.5  / (i*j) / n / m / l
                enddo
              enddo
              cinterp_PRA_ij(g,l) = 0.1 / (i*j) / l 
            enddo
          enddo
        enddo
        do m = AI,AJ
          do l = 1, ADM_lall
            do j = 1, ADM_gall_1d
              do i = 1, ADM_gall_1d_in
                do n = XDIR,ZDIR
                  g = i + ( (j - 1)  * ADM_gall_1d_in)
                  GRD_xr_ij(g,K0,l,m,n) = 0.1 / (n*m) / i / j * l
                enddo
              enddo
            enddo
          enddo
        enddo

        do g = 1, GMTR_T_nmax_var
        do t = TI, TJ
        do l = 1, ADM_lall
        do j = 1, ADM_gall_1d
        do i = 1, ADM_gall_1d_in
          n = i + ( (j - 1)  * ADM_gall_1d_in)
          GMTR_T_var(n,K0,l,t,g) = 1.0 / (i*j*K0*l*t*g) 
        enddo
        enddo
        enddo
        enddo
        enddo

#endif
        return

  end subroutine H_Adv_flux_snap_read

end module mod_src_tracer
