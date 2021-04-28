#include "profiler.h"
program main
      implicit none

#include "parameters.fh"

      ! matrix: num_nonzero x rows      
!fj      real*4,   allocatable :: A(:,:)
!fj      integer*4,allocatable :: LIST(:,:)

      ! input vector, result vector
!fj      real*4, allocatable :: X(:,:)   ! input
!fj      real*4, allocatable :: AX(:,:)  ! result

#ifdef CLEAN_CACHE
      ! for cache clear
      real*4 :: BUF(1572864)
      real*4 :: TMP
#endif
      integer*4 :: ITR

      ! for reading list
!fj      real*4, allocatable :: TA(:)
!fj      real*4, allocatable :: S(:)
!fj      real*4, allocatable :: AS(:)
!fj      real*4, allocatable :: TS(:)
!fj      integer*4, allocatable :: ITPCRS(:)
      integer*4 :: I, J, ICRS

!fj>
      real(kind=4),dimension(NZ,NP) :: A
      COMMON /a/A
!fj      integer(kind=4),dimension(NZ,NP) :: LIST
      real(kind=4),dimension(NV,0:NP_ORG) :: X
      COMMON /x/X
      real(kind=4),dimension(NV,NP) :: AX
      COMMON /ax/AX

      integer(kind=4),dimension(NZ*NP) :: ITPCRS
      COMMON /itpcrs/ITPCRS

      real(kind=4)         :: RES, ANS
      real(kind=8)         :: IRE, IAN, IER

    !cx call write_data_file()
    !cx stop

    call read_data_file()

!fj<
      !---------------------------------
      ! read list file
      !---------------------------------
!fj      open(101,file='inputarray.dat', form='unformatted', action='read')
!fj      read(101) NP              ! read array size
!fj      write(*,"('NP = ',I8)") NP

!fj      allocate(A(NZ,NP))
!fj      allocate(X(NV,NP))
!fj      allocate(AX(NV,NP))
      A(:,:)  = 1.0
      X(:,:)  = 2.0
      AX(:,:) = 3.0
!fj      allocate(LIST(NZ,NP))
!fj      allocate(TA(NP*NZ))       ! allocate array TA
!fj      allocate(ITPCRS(NP*NZ))   ! allocate array ITPCRS
!     allocate(S(NP))           ! allocate array S

!fj      read(101) TA(1:NP*NZ)     ! read array TA
!fj      read(101) ITPCRS(1:NP*NZ) ! read array ITPCRS
!     read(101) S(NP)           ! read array S
!fj      close(101)
      
      !---------------------------------
      ! convert list into suitable form
      !---------------------------------
!fj      ICRS=0
!fj!ocl serial      
!fj      do I=1,NP
!fj         do J=1,NZ
!fj            ICRS=ICRS+1
!fj            LIST(J,I)=ITPCRS(ICRS)
!fj         enddo
!fj      enddo
!fj#ifdef IDEAL
!fj      LIST(:,:) = 1
!fj#endif      

    PROF_INIT
    PROF_START_ALL

      !---------------------------------
      ! do the kernel
      !---------------------------------
!fj      do ITR=1,100
!fj>
      do ITR=1,NITR
!fj<
#ifdef CLEAN_CACHE
         call clear_cache(BUF, TMP)
#endif

         if(NV .EQ. 8) then
!fj            call ax04_8(A, AX, X, NP, LIST)
!fj>
            call ax04_8(A, AX, X, NP, NP_ORG, ITPCRS)
!fj<

         else if(NV .EQ. 6) then
!fj            call ax04_6(A, AX, X, NP, LIST)
!fj>
            call ax04_6(A, AX, X, NP, NP_ORG, ITPCRS)
!fj<
            
         else if(NV .EQ. 4) then
!fj            call ax04_4(A, AX, X, NP, LIST)
!fj>
            call ax04_4(A, AX, X, NP, NP_ORG, ITPCRS)
!fj<

         else if(NV .EQ. 2) then
!fj            call ax04_2(A, AX, X, NP, LIST)
!fj>
            call ax04_2(A, AX, X, NP, NP_ORG, ITPCRS)
!fj<

         else if(NV .EQ. 1) then
!fj            call ax04_1(A, AX, X, NP, LIST)
!fj>
            call ax04_1(A, AX, X, NP, NP_ORG, ITPCRS)
!fj<

         endif
      enddo
#ifdef CLEAN_CACHE
      write(*,*) TMP
#endif

      PROF_STOP_ALL
      PROF_FINALIZE

      !cx   ANS = 8.421424128000000e+09
      !cx   IRE = DBLE(RES)
      !cx   IAN = DBLE(ANS)
      !cx   IER = 1.d0
      !cx   call report_validation(IRE, IAN, IER)
      call validation()

end program main

#ifdef CLEAN_CACHE
!-----------------------------------------------------------
! clear cache       
!-----------------------------------------------------------      
subroutine clear_cache(BUF, TMP)
      implicit none

      ! argumet
      real*4 :: BUF(1572864)
      real*4 :: TMP

      ! local
      integer*4 :: I

      TMP = 0.0
      do I=1,1572864
         TMP = TMP + BUF(I)
      enddo

end subroutine clear_cache
#endif
