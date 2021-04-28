!  This module is for the vertical implicit scheme of non-hydorostatic model.
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
  integer, public, parameter :: SP1 = kind(0.E0) ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = kind(0.D0) ! Double Precision: kind(0.D0)

#ifdef SINGLE
  integer, public, parameter :: RP = SP1
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

module mod_runconf
  integer, public :: NON_HYDRO_ALPHA    = 1 ! Nonhydrostatic/hydrostatic flag
  integer, public :: TRC_vmax   =  0 ! total number of tracers
  logical, public :: THUBURN_LIM = .true.  ! [add] 20130613 R.Yoshida
end module mod_runconf

module mod_cnst
  use mod_precision
  implicit none
  private
  !++ Public procedure

  !++ Public parameters & variables
  real(RP), public, save :: CNST_EGRAV   = 9.80616_RP  ! Gravitational accerlaration of the Earth [m/s2]

  real(RP), public, save :: CNST_RAIR    =  287.0_RP   ! Gas constant of air

  real(RP), public, save :: CNST_CP      = 1004.5_RP   ! Specific heat of air(consant pressure)
  real(RP), public, save :: CNST_CV      = 1004.5_RP - 287.0_RP !CNST_CP - CNST_RAIR             ! Specific heat of air(consant volume)

end module mod_cnst
!-------------------------------------------------------------------------------

module mod_adm
  implicit none

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
  logical, public, save      :: ADM_have_pl = .true.

  !------ Maximum number of pole regions
  integer, public, parameter :: ADM_rgn_nmax_pl = 2

  !------ Local region number
  integer, public, save      :: ADM_lall = 1
  integer, public, parameter      :: ADM_lall_para = 1
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  integer, public, parameter      :: ADM_lall_pl_para = ADM_rgn_nmax_pl
  !
  !------ Present Local region number ! 2010.4.26 M.Satoh

  logical, public, save :: ADM_have_sgp(ADM_lall_para) = .true. ! region have singlar point?

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
  integer, public, parameter      :: ADM_gmin_para = 2
  integer, public, parameter      :: ADM_gmax_para = 129
  integer, public, parameter      :: ADM_gall_1d_para = 130
#ifdef SINGLE_SIM
  integer, public, parameter      :: ADM_gall_1d_in_para = 65
#else
  integer, public, parameter      :: ADM_gall_1d_in_para = 130
#endif
  integer, public, parameter      :: ADM_gall_para = 16900
  !
  !
  !------ Identifiers of grid points around poles.
  integer, public, save :: ADM_gslf_pl = 1
  integer, public, save :: ADM_gmin_pl = 2
  integer, public, save      :: ADM_gmax_pl = 6  ! [mod] S.Iga 100607
  integer, public, save      :: ADM_gall_pl = 6  ! [mod] S.Iga 100607
  integer, public, parameter      :: ADM_gmax_pl_para = 6  ! [mod] S.Iga 100607
  integer, public, parameter      :: ADM_gall_pl_para = 6  ! [mod] S.Iga 100607
  !
  !------ Vertica grid numbers
  integer, public, save      :: ADM_kmin = 2
  integer, public, save      :: ADM_kmax = 95
  integer, public, save      :: ADM_kall = 96
  integer, public, parameter      :: ADM_kall_para = 96
  integer, public, parameter      :: ADM_kmin_para = 2
  integer, public, parameter      :: ADM_kmax_para = 95

end module mod_adm

module mod_grd
  use mod_precision
  use mod_adm, only: &
        ADM_kall_para

  implicit none
  real(RP), public, save :: GRD_rdgzh(ADM_kall_para)
  real(RP), public, save :: GRD_afac(ADM_kall_para)
  real(RP), public, save :: GRD_bfac(ADM_kall_para)

end module mod_grd

module mod_vmtr
  use mod_precision
  use mod_adm
  implicit none
  private
  !++ Public procedure

  !--- G^1/2 X Gamma^2 at the half level
  real(RP), public, save :: VMTR_GSGAM2H_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

  !--- 1 / Gamma at the integer level
  real(RP), public :: VMTR_RGAM_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

  !--- 1 / Gamma at the half level
  real(RP), public :: VMTR_RGAMH_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

  !--- 1 / (G^1/2 X Gamma^2) at the full level
  real(RP), public :: VMTR_RGSGAM2_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

  !--- 1 / (G^1/2 X Gamma^2) at the half level
  real(RP), public :: VMTR_RGSGAM2H_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

end module mod_vmtr
!-------------------------------------------------------------------------------


module mod_vi

    use mod_precision
    use mod_adm, only: &
        TI  => ADM_TI,  &
        TJ  => ADM_TJ,  &
        AI  => ADM_AI,  &
        AIJ => ADM_AIJ, &
        AJ  => ADM_AJ,  &
        K0  => ADM_KNONE
    use mod_adm, only: &
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
       ADM_lall_para,  &
       ADM_lall_pl_para,  &
       ADM_gall_para,  &
       ADM_gall_pl_para,  &
       ADM_kall_para,  &
       ADM_gall_1d_para,  &
       ADM_gall_1d_in_para,  &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_afac,  &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2_ij,     &
       VMTR_RGSGAM2H_ij,    &
       VMTR_RGAMH_ij,       &
       VMTR_RGAM_ij,        &
       VMTR_GSGAM2H_ij
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CVdry => CNST_CV

  implicit none
  private

  !++ Public procedure
  public :: vi_rhow_solver

  !++ Private parameters & variables
  real(RP), public, save :: Mc_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)
  real(RP), public, save :: Ml_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)
  real(RP), public, save :: Mu_ij(ADM_gall_1d_in_para,ADM_kall_para,ADM_gall_1d_para,ADM_lall_para)

contains

  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  & !--- [INOUT]
       rhogw0, rhogw0_pl, & !--- [IN]
       preg0,  preg0_pl,  & !--- [IN]
       rhog0,  rhog0_pl,  & !--- [IN]
       Sr,     Sr_pl,     & !--- [IN]
       Sw,     Sw_pl,     & !--- [IN]
       Sp,     Sp_pl,     & !--- [IN]
       dt                 )

!    use omp_lib
    implicit none
    external omp_get_num_threads
    integer omp_get_num_threads

    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt
    real(RP), intent(inout) :: rhogw    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(in)    :: rhogw0   (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: Sr       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Sw       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sp       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   ) ! source term for pres at the full level

    real(RP) :: Sall    (ADM_gall_1d_in,ADM_kall)
    real(RP) :: beta    (ADM_gall_1d_in   )
    real(RP) :: gamma   (ADM_gall_1d_in,ADM_kall)

    real(RP) :: alfa
    real(RP) :: CVovRt2   ! Cv / R / dt**2

    integer :: g, k, l
    integer :: nn
    integer :: num_threads
    integer :: onestep,onethread,mod_num
    integer :: nstart1,nstart2
    integer :: nend1,nend2
    integer :: i,j,istart,iend,jstart,jend
    !---------------------------------------------------------------------------
!cx    if (mark_for_snapshots(ip_vi_rhow_solver).ne.0) then
!cx    call snapshot_seq_bin_open ("vi_rhow_solver", ADM_prc_me)
!cx    call vi_rhow_snap_write ( &
!cx       rhogw,  rhogw_pl,  & !--- [INOUT]
!cx       rhogw0, rhogw0_pl, & !--- [IN]
!cx       preg0,  preg0_pl,  & !--- [IN]
!cx       rhog0,  rhog0_pl,  & !--- [IN]
!cx       Sr,     Sr_pl,     & !--- [IN]
!cx       Sw,     Sw_pl,     & !--- [IN]
!cx       Sp,     Sp_pl,     & !--- [IN]
!cx       dt                 )
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_vi_rhow_solver) = 0
!cx    endif

    alfa = real(NON_HYDRO_ALPHA,kind=8)
    CVovRt2 = CVdry / Rdry / (dt*dt)

!$omp parallel default(none), &
!$omp private(num_threads,onestep,onethread,mod_num, &
!$omp nstart1,nend1,nstart2,nend2,nn,i,j,g,k,l,jstart,jend,istart,iend, beta,Sall,gamma ), &
!$omp shared(ADM_gall,ADM_gall_1d,ADM_gall_1d_in,ADM_kall,ADM_lall,ADM_kmin,ADM_kmax, &
!$omp Sr, Sp, Sw, rhogw, rhogw0, rhog0, preg0, &
!$omp VMTR_RGAMH_ij, VMTR_RGSGAM2H_ij,VMTR_RGAM_ij, VMTR_RGSGAM2_ij, VMTR_GSGAM2H_ij,&
!$omp Ml_ij, Mu_ij, Mc_ij, &
!$omp CVovRt2, dt, alfa, GRD_rdgzh, GRAV, GRD_bfac, GRD_afac)
    num_threads = omp_get_num_threads()
!$omp do 
    do nn = 1,num_threads
       nstart1 = 1
       nend1   = ADM_gall_1d

       onestep = (nend1 - nstart1 + 1)
       onethread = onestep / num_threads
       mod_num = mod(onestep,num_threads)

       if(nn <= mod_num) then
         nstart1 = (onethread + 1) * (nn-1) + nstart1
         nend1 = nstart1 + onethread
       else
         nstart1 = (onethread + 1) * mod_num + (onethread * ((nn-1) - mod_num)) + nstart1
         nend1 = nstart1 + onethread - 1
       endif

    do l = 1, ADM_lall
       !--- < calc Sall > ---
       do j = nstart1, nend1
       do k  = ADM_kmin+1, ADM_kmax
       do i = 1, ADM_gall_1d_in  
!OCL PREFETCH_WRITE(Sall(i,k+3),level=2)
!OCL PREFETCH_READ(rhogw0(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(Sw(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(VMTR_RGAMH_ij(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(preg0(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(Sp(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(VMTR_RGSGAM2_ij(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(rhog0(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(Sr(i,k+3,j,l),level=2)
!OCL PREFETCH_READ(VMTR_RGAM_ij(i,k+3,j,l),level=2)
          Sall(i,k) = (   ( rhogw0(i,k,j,  l)*alfa + dt * Sw(i,k,j,  l) ) * VMTR_RGAMH_ij  (i,k,j,  l)**2            &
                      - ( ( preg0 (i,k,j,  l)      + dt * Sp(i,k,j,  l) ) * VMTR_RGSGAM2_ij(i,k,j,  l)               &
                        - ( preg0 (i,k-1,j,l)      + dt * Sp(i,k-1,j,l) ) * VMTR_RGSGAM2_ij(i,k-1,j,l)               &
                        ) * dt * GRD_rdgzh(k)                                                               &
                      - ( ( rhog0 (i,k,j,  l)      + dt * Sr(i,k,j,  l) ) * VMTR_RGAM_ij(i,k,j,  l)**2 * GRD_afac(k) &
                        + ( rhog0 (i,k-1,j,l)      + dt * Sr(i,k-1,j,l) ) * VMTR_RGAM_ij(i,k-1,j,l)**2 * GRD_bfac(k) &
                        ) * dt * 0.5_RP * GRAV                                                               &
                      ) * CVovRt2
       enddo
       enddo

       !--- boundary conditions
       do i = 1, ADM_gall_1d_in  
          rhogw(i,ADM_kmin,j,  l) = rhogw(i,ADM_kmin,j,  l) * VMTR_RGSGAM2H_ij(i,ADM_kmin,j,  l)
          rhogw(i,ADM_kmax+1,j,l) = rhogw(i,ADM_kmax+1,j,l) * VMTR_RGSGAM2H_ij(i,ADM_kmax+1,j,l)
          Sall (i,ADM_kmin+1)   = Sall (i,ADM_kmin+1) - Ml_ij(i,ADM_kmin+1,j,l) * rhogw(i,ADM_kmin,j,  l)
          Sall (i,ADM_kmax  )   = Sall (i,ADM_kmax  ) - Mu_ij(i,ADM_kmax,j,  l) * rhogw(i,ADM_kmax+1,j,l)
       enddo

       !--- < solve tri-daigonal matrix > ---

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do i = 1, ADM_gall_1d_in  
          beta (i)     = Mc_ij(i,k,j,l)
          rhogw(i,k,j,l) = Sall(i,k) / beta(i)
       enddo

       !--- forward
       do k = ADM_kmin+2, ADM_kmax
       do i = 1, ADM_gall_1d_in  
          gamma(i,k)   = Mu_ij(i,k-1,j,l) / beta(i)
          beta (i)     = Mc_ij(i,k,j,l) - Ml_ij(i,k,j,l) * gamma(i,k) ! update beta
          rhogw(i,k,j,l) = ( Sall(i,k) - Ml_ij(i,k,j,l) * rhogw(i,k-1,j,l) ) / beta(i)
       enddo
       enddo

       !--- backward
       do k = ADM_kmax-1, ADM_kmin+1, -1
       do i = 1, ADM_gall_1d_in  
!OCL PREFETCH_WRITE(rhogw(i,k-1,j,l),level=2)
!OCL PREFETCH_READ(VMTR_GSGAM2H_ij(i,k-1,j,l),level=2)
!OCL PREFETCH_READ(gamma(i,k),level=2)
!OCL PREFETCH_WRITE(rhogw(i+32,k,j,l),level=1)
!OCL PREFETCH_READ(VMTR_GSGAM2H_ij(i+32,k,j,l),level=1)
!OCL PREFETCH_READ(gamma(i+32,k+1),level=1)
          rhogw(i,k,j,l) = rhogw(i,k,j,l) - gamma(i,k+1) * rhogw(i,k+1,j,l)
          rhogw(i,k+1,j,l) = rhogw(i,k+1,j,l) * VMTR_GSGAM2H_ij(i,k+1,j,l)
       enddo
       enddo

       !--- return value ( G^1/2 x gam2 )
       do i = 1, ADM_gall_1d_in
          rhogw(i,adm_kmax+1,j,l) = rhogw(i,adm_kmax+1,j,l) * VMTR_GSGAM2H_ij(i,adm_kmax+1,j,l)
          rhogw(i,adm_kmin,j,l) = rhogw(i,adm_kmin,j,l) * VMTR_GSGAM2H_ij(i,adm_kmin,j,l)
          rhogw(i,adm_kmin+1,j,l) = rhogw(i,adm_kmin+1,j,l) * VMTR_GSGAM2H_ij(i,adm_kmin+1,j,l)
       enddo
       enddo
    enddo
    enddo
!$omp end do
!$omp end parallel

    return
  end subroutine vi_rhow_solver

end module mod_vi
