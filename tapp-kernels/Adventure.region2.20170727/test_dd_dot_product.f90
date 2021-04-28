#include "profiler.h"
program main
  implicit none
  integer(kind=  4),parameter :: N = 393216
  integer(kind=  4),parameter :: M = 24
  integer(kind=  4),parameter :: RKP = 8
  integer(kind=4),parameter :: BLEN = 10000
!  real   (kind=RKP),dimension(:,:),allocatable :: matrix_H
!  real   (kind=RKP),dimension(:,:),allocatable :: matrix_L
!  real   (kind=RKP),dimension(:  ),allocatable :: vector_i
!  real   (kind=RKP),dimension(:  ),allocatable :: vector_o
  real   (kind=RKP),dimension(N,M) :: matrix_H
  real   (kind=RKP),dimension(N,M) :: matrix_L
  real   (kind=RKP),dimension(N  ) :: vector_i
  real   (kind=RKP),dimension(  M) :: vector_o
  real   (kind=RKP),dimension(  M) :: vector_oL
  real   (kind=RKP)                :: ref_vector_o = 1.233550469751905e+31_RKP
  real   (kind=RKP)                :: ref_vector_oL = -7.680281892695775e+13_RKP
  integer(kind=  4) :: itr,i,j
  integer(kind=  8) :: iaddr1, iaddr2, iaddr3
  real   (kind=RKP)                :: tmpH,tmpL
  tmpH = 0.0_8
  tmpL = 0.0_8
  !
!  allocate(matrix_H(N,M))
!  allocate(matrix_L(N,M))
!  allocate(vector_i(N  ))
!  allocate(vector_o(  M))
  !
  PROF_INIT
  PROF_START_ALL
  do j=1,M
     do i=1,N
        if(i.eq.1) then
           matrix_H(i,j)=real( 1.0E+00_8,kind=RKP)
           matrix_L(i,j)=real( 0.0E+00_8,kind=RKP)
        else if(i.eq.N) then
           matrix_H(i,j)=real(-1.0E+00_8,kind=RKP)
           matrix_L(i,j)=real( 0.0E+00_8,kind=RKP)
        else
           matrix_H(i,j)=real( 1.0E-18_8,kind=RKP)
           matrix_L(i,j)=real( 0.0E+00_8,kind=RKP)
           matrix_H(i,j) = matrix_H(i,j) * mod(i,5)
        end if
     end do
  end do
  do i=1,N
     vector_i(i)=real(1.0E-18_8,kind=RKP)
     vector_i(i) = vector_i(i)  * mod(i,5)
  end do
   vector_i(1) = real( 1.0E+00_8,kind=RKP)
  !
  do itr=1,2
     !$omp parallel default(shared), private(j)
     !$omp do
     do j=1,M
!        call DD_InnerProduct(N,matrix_H(:,j),matrix_L(:,j),vector_i(:),vector_o(j))
        call DD_InnerProduct(N,matrix_H(:,j),matrix_L(:,j),vector_i(:),vector_o(j),vector_oL(j),BLEN,itr)
     end do
     !$omp end do
     !$omp end parallel
  end do
  !
!  write(6,'(e30.20)') real(vector_o(1),kind=8)

  !write(6,'(e30.20)') vector_o(1),ref_vector_o
  !write(6,'(e30.20)') vector_oL(1)
  call report_validation(vector_o(1),ref_vector_o, 1e-13_8)
  call report_validation(vector_oL(1),ref_vector_oL, 1e-13_8)
  !
  iaddr1=loc(matrix_H(1,1))
  iaddr2=loc(matrix_L(1,1))
  iaddr3=loc(vector_i(1)  )
  !
!  call iaddr_print(iaddr1)
!  call iaddr_print(iaddr2)
!  call iaddr_print(iaddr3)
  !
!  deallocate(matrix_H,matrix_L)
!  deallocate(vector_i,vector_o)
  !
  PROF_STOP_ALL
  PROF_FINALIZE
end program main
!
!
!
subroutine DD_InnerProduct(N,matrix_H,matrix_L,vector_i,vector_o,vector_oL,BLEN,itr)
  implicit none
  integer(kind=4),             intent(IN)  :: N
  real   (kind=8),dimension(N),intent(IN)  :: matrix_H, matrix_L
  real   (kind=8),dimension(N),intent(IN)  :: vector_i
  real   (kind=8)             ,intent(OUT) :: vector_o
  real   (kind=8)             ,intent(OUT) :: vector_oL
  integer(kind=4),intent(IN) :: BLEN,itr
  !
  !
  real   (kind=8),parameter :: QD_SPLITTER = 134217729.0E0_8
  !
  real   (kind=8),dimension(BLEN    +1):: t1_H,t1_L
  !
  real   (kind=8) :: vector_H,vector_L
  real   (kind=8) :: p3,temp0
  real   (kind=8) :: aa_hi,aa_lo
  real   (kind=8) :: temp1
  real   (kind=8) :: bb_hi,bb_lo
  real   (kind=8) :: p1_L1,p1_L2
  real   (kind=8) :: p2_L1,p2_L2,p2_L3
  real   (kind=8) :: a0,a1
  real   (kind=8) :: b0,b1,bb0
  real   (kind=8) :: e_L1,e_L2,e_L3
  real   (kind=8) :: s0,s1,s_L1,s_L2
  integer(kind=4) :: h,i,k,ii,nn
  !
  vector_H = 0.0E0_8
  vector_L = 0.0E0_8
   do i=1,BLEN
      t1_H(i) = 1.0E-18_8
      t1_L(i) = 0.0E0_8
   end do
   t1_H(1) = 1.0_8
   if(itr > 1) then
!$omp barrier
      PROF_START("region2")
   endif
  do ii=1,N,BLEN
     nn=min(BLEN,N-ii+1)
     t1_H(nn+1) = vector_H
     t1_L(nn+1) = vector_L
     !
     h=nn+1
     do
        !
        k=h/2
!ocl norecurrence
        do i=1,k
           a0 = t1_H(    i)
           a1 = t1_L(    i)
           b0 = t1_H(h-k+i)
           b1 = t1_L(h-k+i)
           s0 = a0 + b0
           bb0 = s0 - a0
           e_L1 = ( a0 - (s0 - bb0) ) + (b0 - bb0)
           s_L1 = s0
           e_L2 = e_L1 + (a1 + b1)
           s1 = s_L1 + e_L2
           e_L3 = e_L2 - (s1 - s_L1)
           s_L2 = s1
           t1_H(i) = s_L2
           t1_L(i) = e_L3
        end do
        if(h-k*2.eq.1) then
           h=k+1
        else
           h=k
        end if
        if(h.eq.1) then
           vector_H = t1_H(1)
           vector_L = t1_L(1)
           exit
        end if
        !
     end do
     !
  end do
  if(itr > 1) then
!$omp barrier
     PROF_STOP("region2")
  endif
  
  vector_o = vector_H
  vector_oL = vector_L
  !
  return
end subroutine DD_InnerProduct


subroutine report_validation (computed_result, expected_result, error_bar)
    real(kind=8) :: computed_result, expected_result, error_bar
    real(kind=8) :: error_norm

    error_norm = computed_result/expected_result
    if ( 0.999 < error_norm .and. error_norm < 1.001 ) then
        write(*,*) "the computed result seems to be OK."
    else
     write(*,*) "the computed result is not close enough to the expected value." 
    endif
end subroutine

