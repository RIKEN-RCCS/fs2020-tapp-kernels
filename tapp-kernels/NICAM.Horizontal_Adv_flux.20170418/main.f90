#include "profiler.h"
program main
    call dynamics_step
    stop
end program main

  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
        H_Adv_flux_snap_read, &
        horizontal_flux
    use mod_adm, only: &
        ADM_gall,    &
        ADM_gall_pl, &
        ADM_gall_1d,    &
        ADM_gall_1d_in,    &
        ADM_gmin,    &
        ADM_gmax,    &
        ADM_lall,    &
        ADM_lall_pl, &
        ADM_kall
    use mod_adm, only: AI=>ADM_AI, AJ=>ADM_AJ, K0=>ADM_KNONE
    use mod_grd, only: XDIR=>GRD_XDIR, YDIR=>GRD_YDIR, ZDIR=>GRD_ZDIR,GRD_xr_ij
    use mod_oprt, only:  cinterp_HN_ij

    implicit none
    integer :: i

    real(RP) flx_h_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,XDIR:ZDIR)
    real(RP) rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) dt

!    real(RP) flx_h_ij    (6,ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall   )
!    real(RP) flx_h_ij    (6,ADM_gall,ADM_kall,ADM_lall   )
    real(RP) flx_h_ij    (ADM_gall,ADM_kall,ADM_lall,6)
!    real(RP) GRD_xc_ij (XDIR:ZDIR,ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall,AI:AJ)
!    real(RP) GRD_xc_ij (ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall,AI:AJ,XDIR:ZDIR)
    real(RP) GRD_xc_ij (ADM_gall,ADM_kall,ADM_lall,AI:AJ,XDIR:ZDIR)
    real(RP) rho_ij      (ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall   )
    real(RP) rhovx_ij    (ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall   )
    real(RP) rhovy_ij    (ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall   )
    real(RP) rhovz_ij    (ADM_gall_1d_in,ADM_gall_1d,ADM_kall,ADM_lall   )
#ifdef OUTPUTCHECK
    real(8) :: tmp_flx_h = 0.0_8
    real(8) :: tmp_flx_h_pl = 0.0_8
    real(8) :: tmp_GRD_xc = 0.0_8
    real(8) :: tmp_GRD_xc_pl = 0.0_8
    real(8) :: ref_flx_h = 0.4168435964941339_8
    real(8) :: ref_GRD_xc = 4135.705156248743_8
    integer :: j,k,l,g,n,m
#endif
#ifdef USE_TIMER
    real*8 t_all_0, t_all, t_kernel_0, t_kernel
    t_all_0 = 0.0
    t_all = 0.0
    t_kernel_0 = 0.0
    t_kernel = 0.0
#endif
    PROF_INIT
    PROF_START_ALL
#ifdef USE_TIMER
    call gettod(t_all_0)
#endif

    flx_h_ij = 0.0_RP
    GRD_xc_ij = 0.0_RP
    flx_h_pl = 0.0_RP
    GRD_xc_pl = 0.0_RP

    call H_Adv_flux_snap_read ( &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

#ifdef OUTPUTCHECK
    do l=1,ADM_lall
      do k=1,ADM_kall
        do j=1,ADM_gall_1d
          do i=1,ADM_gall_1d_in
            g = i + ( (j - 1)  * ADM_gall_1d_in)
            rho_ij(i,j,k,l)   = rho(g,k,l)
            rhovx_ij(i,j,k,l) = rhovx(g,k,l)
            rhovy_ij(i,j,k,l) = rhovy(g,k,l)
            rhovz_ij(i,j,k,l) = rhovz(g,k,l)
          end do
        end do
      end do
    end do
#endif
    PROF_START("horizontal_flux")
#ifdef USE_TIMER
    call gettod(t_kernel_0)
#endif
!    do i=1,10
    do i=1,1
    call horizontal_flux( &
       flx_h_ij,  flx_h_pl,  &
       GRD_xc_ij, GRD_xc_pl, &
!       rho_ij,    rho_pl,    &
!       rhovx_ij,  rhovx_pl,  &
!       rhovy_ij,  rhovy_pl,  &
!       rhovz_ij,  rhovz_pl,  &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt,GRD_xr_ij,cinterp_HN_ij )
    end do
#ifdef USE_TIMER
    call gettod(t_kernel)
    t_kernel = t_kernel - t_kernel_0
#endif
    PROF_STOP("horizontal_flux")
#ifdef OUTPUTCHECK
    do l=1, ADM_lall
     do k=1, ADM_kall
      do j=ADM_gmin, ADM_gmax
      do i=ADM_gmin, ADM_gmax
        g = i + ( (j - 1)  * ADM_gall_1d_in)
        do n = 1, 6
!          tmp_flx_h = tmp_flx_h + (flx_h_ij(n,i,j,k,l) ** 2)
!          tmp_flx_h = tmp_flx_h + (flx_h_ij(n,g,k,l) ** 2)
          tmp_flx_h = tmp_flx_h + (flx_h_ij(g,k,l,n) ** 2)
        enddo
        do n = AI,AJ
         do m = XDIR,ZDIR
           !tmp_GRD_xc = tmp_GRD_xc + (GRD_xc_ij(m,i,j,k,l,n) ** 2)
           tmp_GRD_xc = tmp_GRD_xc + (GRD_xc_ij(g,k,l,n,m) ** 2)
         enddo
       enddo
      enddo
      enddo
     enddo
    enddo
    do l=1, ADM_lall_pl
     do k=1, ADM_kall
      do g=1, ADM_gall_pl
        tmp_flx_h_pl = tmp_flx_h_pl + (flx_h_pl(g,k,l) ** 2)
         do n = XDIR,ZDIR
           tmp_GRD_xc_pl = tmp_GRD_xc_pl + (GRD_xc_pl(g,k,l,n) ** 2)
         enddo
      enddo
     enddo
    enddo

!   write(*,*) "OUTPUTCHECK RMS flx_h=",sqrt(tmp_flx_h),"flx_h_pl=", &
!               sqrt(tmp_flx_h_pl),"GRD_xc=",sqrt(tmp_GRD_xc),"GRD_xc_pl=",sqrt(tmp_GRD_xc_pl)

    tmp_flx_h = sqrt(tmp_flx_h)
    tmp_GRD_xc = sqrt(tmp_GRD_xc)

    call report_validation(tmp_flx_h, ref_flx_h, 3e-6_8)
    call report_validation(tmp_GRD_xc, ref_GRD_xc, 3e-6_8)

#endif

#ifdef USE_TIMER
    call gettod(t_all)
    t_all = t_all - t_all_0
#endif
  PROF_STOP_ALL
  PROF_FINALIZE 
#ifdef USE_TIMER
    write(*,*)"t_kernel ",t_kernel
    write(*,*)"t_all ",t_all
#endif
  
    return
  end subroutine dynamics_step




