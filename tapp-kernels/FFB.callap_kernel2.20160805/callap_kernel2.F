#include "profiler.h"

      subroutine callap_kernel2(MP, NP, ME, NE, IENP,
     *                          FXYZ, S, DNXP, DNYP, DNZP)
      implicit none

      ! argument
      integer(4),intent(in) :: MP, NP, ME, NE
      integer(4),intent(in),dimension(8,MP) :: IENP
      real(4),intent(inout),dimension(3,MP) :: FXYZ
      real(4),intent(in),dimension(0:ME) :: S
      real(4),intent(in),dimension(8,MP) :: DNXP
      real(4),intent(in),dimension(8,MP) :: DNYP
      real(4),intent(in),dimension(8,MP) :: DNZP

      ! local
      integer(4) :: IP, IE, I
      real(4)    :: SWRK
      real(4)    :: FXBUF, FYBUF, FZBUF

      integer(kind=4)      :: NP3
      real(kind=4)         :: RES, ANS
      real(kind=8)         :: IRE, IAN, IER

      PROF_START("1")
!$omp parallel default(none) private(IP,I,FXBUF,FYBUF,FZBUF,IE,SWRK) 
!$omp&shared(NP,IENP,S,DNXP,DNYP,DNZP,FXYZ)
!$omp do
!ocl simd
!ocl swp
      do IP=1,NP
         FXBUF = 0.0
         FYBUF = 0.0
         FZBUF = 0.0
         do I=1,8
            IE = IENP(I,IP)
            SWRK = S(IE)
         
            FXBUF = FXBUF - SWRK * DNXP(I,IP)
            FYBUF = FYBUF - SWRK * DNYP(I,IP)
            FZBUF = FZBUF - SWRK * DNZP(I,IP)
         enddo
         FXYZ(1,IP) = FXYZ(1,IP) + FXBUF
         FXYZ(2,IP) = FXYZ(2,IP) + FYBUF
         FXYZ(3,IP) = FXYZ(3,IP) + FZBUF
      enddo
!$omp end do
!$omp end parallel
      PROF_STOP("1")

      NP3 = NP * 3
      call get_ss_r4(FXYZ, NP3, RES)
#if 0
      write (*,*) "get_ss_r4=", RES
#else
      ANS = 1.56653343E+10
      IRE = DBLE(RES)
      IAN = DBLE(ANS)
      IER = 1.d0
      call report_validation(IRE, IAN, IER)
#endif

      end subroutine callap_kernel2
      
