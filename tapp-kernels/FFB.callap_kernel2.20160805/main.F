CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "profiler.h"

      PROGRAM LES3X
      implicit none

      ! num node and element
      integer(4) :: MP, NP  ! node(max and actual)
      integer(4) :: ME, NE  ! elem(max and actual)

      ! misc
      integer(4) :: I
#ifdef CLEAN_CACHE
      integer(4) :: NITR
#endif

      parameter(MP = 1093727)
      parameter(NP = 1092727)
      parameter(ME = 1062208)
      parameter(NE = 1061208)
#ifdef CLEAN_CACHE
      parameter(NITR = 1) ! 11880
#endif

      ! element value
      real(4),dimension(0:ME) :: S
      real(4),dimension(8,MP) :: DNXP,DNYP,DNZP

      ! node value
      real(4),dimension(3,MP) :: FXYZ                 ! calculat value

#ifdef CLEAN_CACHE
      ! L2 cache cleaning
      integer(4) :: SZDUMMY
      parameter(SZDUMMY=1572864)  ! 6MB
      integer(4),dimension(SZDUMMY) :: DUMMY
#endif

      INTEGER(4),DIMENSION(8,MP) :: ELEM1
      COMMON /data_elem1/ELEM1

      character*20 :: filename
      integer :: iunit, is_ok
      iunit=77
      write(filename,'(a)') "data_file_ELEM1"

      !cx   open(iunit, file=filename, form="formatted", status="new",
      open(iunit, file=filename, form="formatted", status="old",
     &      iostat=is_ok)
      if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
      endif
      !cx   call sub_idata_write(iunit, 8*MP, ELEM1)
      call sub_idata_read(iunit, 8*MP, ELEM1)
      close(iunit)

      PROF_INIT
      PROF_START_ALL

      FXYZ = 0.0
      DNXP = 1.0
      DNYP = 2.0
      DNZP = 3.0
      S    = 4.0
      
#ifdef CLEAN_CACHE
      ! dummy array to clear cache for each iteration
      DUMMY = 0

      do I=1, NITR
         DUMMY(:) = 0.0 ! clear cache memory by reading 6MB array
#endif

         ! do the kernel
         call callap_kernel2(MP, NP, ME, NE, ELEM1, 
     *                       FXYZ, S, DNXP, DNYP, DNZP)
#ifdef CLEAN_CACHE
      end do
#endif

      PROF_STOP_ALL
      PROF_FINALIZE
      
      END program les3x
