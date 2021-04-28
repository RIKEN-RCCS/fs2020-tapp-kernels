#include "profiler.h"

program main
    call dynamics_step
    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
!        H_Adv_limiter_snap_write, &
!        H_Adv_limiter_snap_read, &
        horizontal_limiter_thuburn
    use mod_adm, only: &
!fj 2016 >
        ADM_gall_1d,    &
        ADM_gall_1d_in,    &
        ADM_gmin,    &
        ADM_gmax,    &
!fj 2016 <
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    integer :: i

!!!    real(RP) q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) q_a     (ADM_gall   ,ADM_kall,ADM_lall,6   )
    real(RP) q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) q       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!!!    real(RP) ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) ch      (ADM_gall   ,ADM_kall,ADM_lall,6   )
    real(RP) ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP) cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
!!!    real(RP) cmask   (3,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) cmask   (ADM_gall   ,ADM_kall,ADM_lall, 3   )
!fj<
    real(RP) cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer :: l,g,k
!fj 2016 >
    integer :: g1,g2
#ifdef OUTPUTCHECK
    real(8) :: tmp_q_a = 0.0_8
    real(8) :: ref_q_a = 0.2534228124694846_8
#endif
!fj 2016 <
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
#ifdef OUTPUTCHECK

    do i = 1,6
      do l=1, ADM_lall
        do k=1, ADM_kall
          do g=1, ADM_gall
           q_a(g,k,l,i) = 1.1 / g/ k / l / i
           ch(g,k,l,i) = 1.4 / g/ k / l / i
          enddo
        enddo
      enddo
    enddo

    do l=1, ADM_lall
      do k=1, ADM_kall
        do g=1, ADM_gall
!!!          do i = 1,6
!!!           q_a(i,g,k,l) = 1.1 / g/ k / l / i
!!!           ch(i,g,k,l) = 1.4 / g/ k / l / i
!!!          enddo
           q(g,k,l) = 1.2 / g/ k / l
           d(g,k,l) = 1.3 / g/ k / l
        enddo
      enddo
    enddo
    do i= 1,3
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do g = 1,ADM_gall
          cmask(g,k,l,i)=0.0_RP
       enddo
    enddo
    enddo
    enddo
#endif

!    call H_Adv_limiter_snap_read &
!                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
    call horizontal_limiter_thuburn &
                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
     
    PROF_START("horizontal_limiter_thuburn")
#ifdef USE_TIMER
    call gettod(t_kernel_0)
#endif

    !cx call fapp_start( '____Horizontal_Adv_limiter', 1, 1 )
!    do i=1,1010
    do i=1,1
    call horizontal_limiter_thuburn &
                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
    end do
    !cx call fapp_stop( '____Horizontal_Adv_limiter', 1, 1 )
#ifdef USE_TIMER
    call gettod(t_kernel)
    t_kernel = t_kernel - t_kernel_0
#endif
    PROF_STOP("horizontal_limiter_thuburn")

!fj 2016 output check q_a
#ifdef OUTPUTCHECK
    do l=1, ADM_lall
     do k=1, ADM_kall
!      do g=1, ADM_gall
      do g1=ADM_gmin, ADM_gmax
#ifdef SINGLE_SIM
      do g2=ADM_gmin, ADM_gmax/2
#else
      do g2=ADM_gmin, ADM_gmax
#endif
        g = ADM_gall_1d_in * (g1 -1) + g2
        do i = 1, 6
!!!          tmp_q_a = tmp_q_a + dble(q_a(i,g,k,l) ** 2)
          tmp_q_a = tmp_q_a + dble(q_a(g,k,l,i) ** 2)
        enddo
      enddo
      enddo
     enddo
    enddo

!    write(*,*) "OUTPUTCHECK q_a=",sqrt(tmp_q_a)
    tmp_q_a = sqrt(tmp_q_a)
    call report_validation(tmp_q_a, ref_q_a, 9e-7_8)
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

