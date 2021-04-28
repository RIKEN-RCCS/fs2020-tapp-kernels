#include "profiler.h"
      program main
      integer,parameter:: np = 12
      integer,parameter:: n = 565285
      integer,parameter:: ne = 408384
      integer,parameter:: kd = 2
      integer,parameter:: n_size = 565285
      integer,parameter:: ne_size = 408384
      integer,parameter:: nab = 0
      integer,parameter:: ncolor = 1

      ! in
      integer cny(10,ne_size),num(ne_size),npl(np+1)
      integer nnum(np),nlist(n,np),nabc(6,nab)
      integer nab_ista(0:np-1),nab_iend(0:np-1)
      integer color_ind(np*MAXCOLOR)
      REAL_4 ug(3,n_size)
      REAL_4 coor(3,n_size),rmat(kd,10),younglst(kd,2)
      REAL_4 coe_merge(4,ne_size),alp,bet,dt,me(10,10)
      REAL_4 abc(3,6,6,nab)
      logical addglobalvector
      ! out
      REAL_4 rg(3,n_size)
      ! temp
      REAL_4 rrtmp(3,n_size,np)
      logical true(15)
c$$$      character(len=7) fname
      common /cblk_cny/ cny
      common /cblk_num/ num

c$$$      write(fname,'(a4,i3.3)') 'xxx.',10
c$$$      open(unit=999,file=fname,form='UNFORMATTED')
c$$$      read(unit=999) cny,num
c$$$      close(999)
c$$$      write(fname,'(a4,i3.3)') 'xxx.',10
c$$$      open(unit=999,file=fname,form='UNFORMATTED')
c$$$      read(unit=999) cny,num,npl,nnum,nlist,nabc,color_ind
c$$$      close(999)

      character*20 :: filename
      integer :: iunit, is_ok
      iunit=77
      write(filename,'(a)') "data_file_77"

      open(iunit, file=filename, form="formatted", status="old",
     &      iostat=is_ok)
      if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
      endif
!cx      call sub_idata_write(iunit, ne_size, num)
!cx      call sub_idata_write(iunit, 10*ne_size, cny)
      call sub_idata_read(iunit, ne_size, num)
      call sub_idata_read(iunit, 10*ne_size, cny)
      close(iunit)

      PROF_INIT
      PROF_START_ALL

      nab_ista = 1
      nab_iend = 0
      alp = 0.1
      bet = 0.1
      dt = 2.0
      ug = 0.1
!      coor = 1.0
!ocl nosimd
      do i=1,n
      coor(1,i)=sin(real(mod(2*i,56000)))
      coor(2,i)=sin(real(mod(3*i,57000)))
      coor(3,i)=sin(real(mod(5*i,58000)))
      enddo
      rmat = 0.1
      younglst = 2.0
      coe_merge = 0.1
      me = 0.1
      abc = 0.1
      true = .true.

c$$$      color_ind(1) = 0
c$$$      color_ind(2) = 28768
c$$$      color_ind(3) = 54860
c$$$      color_ind(4) = 80798
c$$$      color_ind(5) = 107792
c$$$      color_ind(6) = 133844
c$$$      color_ind(7) = 160736
c$$$      color_ind(8) = 190608
c$$$      color_ind(9) = 215791
c$$$      color_ind(10) = 240888
c$$$      color_ind(11) = 267798
c$$$      color_ind(12) = 293683
c$$$      color_ind(13) = 322375

      color_ind(1)  = 0     
      color_ind(2)  = 17848 
      color_ind(3)  = 34034 
      color_ind(4)  = 50126 
      color_ind(5)  = 66874 
      color_ind(6)  = 83036 
      color_ind(7)  = 99720 
      color_ind(8)  = 118252
      color_ind(9)  = 133876
      color_ind(10) = 149446
      color_ind(11) = 166140
      color_ind(12) = 182200
      color_ind(13) = 200000

      PROF_START('TIMER_COMP_MATVEC_IF')
!ocl serial
      do i=1, 1000
         call cal_amat_s_color
     &        (ncolor,color_ind,n_size,ne_size,np,n,ne,kd,nab,
     &        ug,rg,
     &        cny,num,npl,
     &        coor,rmat,younglst,coe_merge,alp,bet,dt,me,
     &        nnum,nlist,nabc,abc,
     &        nab_ista,nab_iend,rrtmp,addglobalvector,true)
      enddo
      PROF_STOP('TIMER_COMP_MATVEC_IF')
#ifdef RESULT_CHECK
      call result_check(rg)
#endif
      PROF_STOP_ALL
      PROF_FINALIZE
      
      end
      
#ifdef RESULT_CHECK
      subroutine result_check(rg)
      integer,parameter:: n_size = 565285
      REAL_4 rg(3,n_size)
      real*8 res
      real*8 :: computed_result, expected_result, error_norm 

      res = 0
      do j=1, n_size
         do i=1, 3
            res = res+dble(rg(i,j)*rg(i,j))
         enddo
      enddo
      res = res / dble(3*n_size)
      res = sqrt(res)
!cx !      call report_validation(res, 0.0007517793566033252_8, 0.001_8)
!cx       call report_validation(res, 0.0007232123085724831_8, 0.001_8)

        expected_result =0.0007232123085724831_8
        computed_result = res
        error_norm = computed_result/expected_result

        if ( 0.999 < error_norm .and. error_norm < 1.001 ) then
            write(*,*) "the computed result seems to be OK."
        else
        !cx    write(*,*) "the computed result is not close enough to the expected value." 
        write(*,*) "the computed result is not close enough." 
        endif

      end subroutine
#endif
