!-------------------------------------------------------------------------------
!> module streamlike
!!
!! @par Description
!!          common computational pattern in NICAM
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_streamlike
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: streamlike_pattern1
  public :: streamlike_pattern2
  public :: streamlike_pattern3

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine streamlike_pattern1( &
       gall,    &
       kall,    &
       lall,    &
       qall,    &
       metrics, &
       PROG,    &
       PROGq,   &
       rho,     &
       vx,      &
       vy,      &
       vz,      &
       ein,     &
       q        )
    implicit none

    integer,  intent(in)  :: gall
    integer,  intent(in)  :: kall
    integer,  intent(in)  :: lall
    integer,  intent(in)  :: qall
    real(RP), intent(in)  :: metrics(gall,kall,lall)
    real(RP), intent(in)  :: PROG   (gall,kall,lall,6)
    real(RP), intent(in)  :: PROGq  (gall,kall,lall,qall)
    real(RP), intent(out) :: rho    (gall,kall,lall)
    real(RP), intent(out) :: vx     (gall,kall,lall)
    real(RP), intent(out) :: vy     (gall,kall,lall)
    real(RP), intent(out) :: vz     (gall,kall,lall)
    real(RP), intent(out) :: ein    (gall,kall,lall)
    real(RP), intent(out) :: q      (gall,kall,lall,qall)

    integer  :: g, k, l, iq
    !---------------------------------------------------------------------------

    do l = 1, lall
    !$OMP PARALLEL DO
    do k = 1, kall
    do g = 1, gall
       rho(g,k,l) = PROG(g,k,l,1) / metrics(g,k,l)
       vx (g,k,l) = PROG(g,k,l,2) / PROG(g,k,l,1)
       vy (g,k,l) = PROG(g,k,l,3) / PROG(g,k,l,1)
       vz (g,k,l) = PROG(g,k,l,4) / PROG(g,k,l,1)
       ein(g,k,l) = PROG(g,k,l,5) / PROG(g,k,l,1)
    enddo
    enddo
    enddo

    do iq = 1, qall
    do l  = 1, lall
    !$OMP PARALLEL DO
    do k  = 1, kall
    do g  = 1, gall
       q(g,k,l,iq) = PROGq(g,k,l,iq) / PROG(g,k,l,1)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine streamlike_pattern1

  !-----------------------------------------------------------------------------
  subroutine streamlike_pattern2( &
       gall,        &
       kall,        &
       lall,        &
       tendency_0,  &
       tendency_11, &
       tendency_12, &
       tendency_13, &
       tendency_14, &
       tendency_15, &
       tendency_16, &
       tendency_21, &
       tendency_22, &
       tendency_23, &
       tendency_24, &
       tendency_25, &
       tendency_26, &
       tendency     )
    implicit none

    integer,  intent(in)  :: gall
    integer,  intent(in)  :: kall
    integer,  intent(in)  :: lall
    real(RP), intent(in)  :: tendency_0 (gall,kall,lall,6)
    real(RP), intent(in)  :: tendency_11(gall,kall,lall)
    real(RP), intent(in)  :: tendency_12(gall,kall,lall)
    real(RP), intent(in)  :: tendency_13(gall,kall,lall)
    real(RP), intent(in)  :: tendency_14(gall,kall,lall)
    real(RP), intent(in)  :: tendency_15(gall,kall,lall)
    real(RP), intent(in)  :: tendency_16(gall,kall,lall)
    real(RP), intent(in)  :: tendency_21(gall,kall,lall)
    real(RP), intent(in)  :: tendency_22(gall,kall,lall)
    real(RP), intent(in)  :: tendency_23(gall,kall,lall)
    real(RP), intent(in)  :: tendency_24(gall,kall,lall)
    real(RP), intent(in)  :: tendency_25(gall,kall,lall)
    real(RP), intent(in)  :: tendency_26(gall,kall,lall)
    real(RP), intent(out) :: tendency   (gall,kall,lall,6)

    integer  :: g, k, l
    !---------------------------------------------------------------------------

    do l = 1, lall
    !$OMP PARALLEL DO
    do k = 1, kall
    do g = 1, gall
       tendency(g,k,l,1) = tendency_0(g,k,l,1) + tendency_11(g,k,l) + tendency_21(g,k,l)
       tendency(g,k,l,2) = tendency_0(g,k,l,2) + tendency_12(g,k,l) + tendency_22(g,k,l)
       tendency(g,k,l,3) = tendency_0(g,k,l,3) + tendency_13(g,k,l) + tendency_23(g,k,l)
       tendency(g,k,l,4) = tendency_0(g,k,l,4) + tendency_14(g,k,l) + tendency_24(g,k,l)
       tendency(g,k,l,5) = tendency_0(g,k,l,5) + tendency_15(g,k,l) + tendency_25(g,k,l)
       tendency(g,k,l,6) = tendency_0(g,k,l,6) + tendency_16(g,k,l) + tendency_26(g,k,l)
    enddo
    enddo
    enddo

!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,1) = tendency_0(g,k,l,1) + tendency_11(g,k,l) + tendency_21(g,k,l)
!    enddo
!    enddo
!    enddo
!
!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,2) = tendency_0(g,k,l,2) + tendency_12(g,k,l) + tendency_22(g,k,l)
!    enddo
!    enddo
!    enddo
!
!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,3) = tendency_0(g,k,l,3) + tendency_13(g,k,l) + tendency_23(g,k,l)
!    enddo
!    enddo
!    enddo
!
!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,4) = tendency_0(g,k,l,4) + tendency_14(g,k,l) + tendency_24(g,k,l)
!    enddo
!    enddo
!    enddo
!
!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,5) = tendency_0(g,k,l,5) + tendency_15(g,k,l) + tendency_25(g,k,l)
!    enddo
!    enddo
!    enddo
!
!    do l = 1, lall
!    !$OMP PARALLEL DO
!    do k = 1, kall
!    do g = 1, gall
!       tendency(g,k,l,6) = tendency_0(g,k,l,6) + tendency_16(g,k,l) + tendency_26(g,k,l)
!    enddo
!    enddo
!    enddo

    return
  end subroutine streamlike_pattern2

  !-----------------------------------------------------------------------------
  subroutine streamlike_pattern3( &
       gall,     &
       kall,     &
       lall,     &
       tendency, &
       value,    &
       fraction  )
    implicit none

    integer,  intent(in)    :: gall
    integer,  intent(in)    :: kall
    integer,  intent(in)    :: lall
    real(RP), intent(in)    :: tendency(gall,kall,lall,6)
    real(RP), intent(inout) :: value   (gall,kall,lall,6)
    real(RP), intent(in)    :: fraction

    integer  :: g, k, l
    !---------------------------------------------------------------------------

    do l = 1, lall
    !$OMP PARALLEL DO
    do k = 1, kall
    do g = 1, gall
       value(g,k,l,1) = value(g,k,l,1) + tendency(g,k,l,1) * fraction
       value(g,k,l,2) = value(g,k,l,2) + tendency(g,k,l,2) * fraction
       value(g,k,l,3) = value(g,k,l,3) + tendency(g,k,l,3) * fraction
       value(g,k,l,4) = value(g,k,l,4) + tendency(g,k,l,4) * fraction
       value(g,k,l,5) = value(g,k,l,5) + tendency(g,k,l,5) * fraction
       value(g,k,l,6) = value(g,k,l,6) + tendency(g,k,l,6) * fraction
    enddo
    enddo
    enddo

    return
  end subroutine streamlike_pattern3

end module mod_streamlike
