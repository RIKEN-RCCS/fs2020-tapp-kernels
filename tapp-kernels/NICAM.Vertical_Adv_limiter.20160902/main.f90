#include "profiler.h"

program main

    call dynamics_step
    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
        vertical_limiter_thuburn
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_1d,    &
       ADM_gall_1d_in,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    integer :: i

    real(RP) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(RP) :: q_h_ij   (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,   ADM_lall   )
    real(RP) :: q_ij     (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,   ADM_lall   )
    real(RP) :: d_ij     (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,   ADM_lall   )
    real(RP) :: ck_ij    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,   ADM_lall   ,2)
#ifdef USE_TIMER
    real*8 t_all_0, t_all, t_kernel_0, t_kernel
#endif
#ifdef OUTPUTCHECK
    integer :: j,k,l,g
    real(8) :: tmp_q_h = 0.0
    real(8) :: ref_q_h = 1.729361169265623_8
#endif

    PROF_INIT
    PROF_START_ALL
#ifdef USE_TIMER
    call gettod(t_all_0)
#endif

#ifdef OUTPUTCHECK
    do l=1, ADM_lall
      do j=1, ADM_gall_1d
       do k=1, ADM_kall
        do i=1, ADM_gall_1d_in
           g = ADM_gall_1d_in * (j -1) + i
           q_h_ij(i,k,j,l) = 1.1 / g / k / l
           q_ij(i,k,j,l) = 1.2 / g / k / l
           d_ij(i,k,j,l) = 1.3 / g / k / l
           ck_ij(i,k,j,l,1) = 1.4 / g / k / l
           ck_ij(i,k,j,l,2) = 1.5 / g / k / l
        enddo
       enddo
      enddo
    enddo
#endif

!    call vertical_limiter_thuburn (q_h_ij, q_h_pl, q_ij, q_pl, d_ij, d_pl, ck_ij, ck_pl)

    PROF_START("vertical_limiter_thuburn")
#ifdef USE_TIMER
    call gettod(t_kernel_0)
#endif
!    do i=1,2020
    do i=1,2
    call vertical_limiter_thuburn (q_h_ij, q_h_pl, q_ij, q_pl, d_ij, d_pl, ck_ij, ck_pl)
    end do
#ifdef USE_TIMER
    call gettod(t_kernel)
    t_kernel= t_kernel - t_kernel_0
#endif
    PROF_STOP("vertical_limiter_thuburn")

#ifdef OUTPUTCHECK
!ocl serial
    do l=1, ADM_lall
      do k=1, ADM_kall
!        do g=1, ADM_gall
         do j=1,ADM_gall_1d
          do i=1,ADM_gall_1d_in
           tmp_q_h = tmp_q_h + dble((q_h_ij(i,k,j,l)) ** 2)
        enddo
        enddo
      enddo
    enddo
!    write(*,*) "OUTPUTCHECK q_h=",sqrt(tmp_q_h)
    tmp_q_h = sqrt(tmp_q_h)
    call report_validation(tmp_q_h,ref_q_h, 0.000002_8)
#endif
#ifdef USE_TIMER
    call gettod(t_all)
    t_all= t_all - t_all_0
#endif
    PROF_STOP_ALL
    PROF_FINALIZE
#ifdef USE_TIMER
    write(*,*)"t_kernel ",t_kernel
    write(*,*)"t_all ",t_all
#endif


    return
  end subroutine dynamics_step

