#include "profiler.h"
!-----------------------------------------------------------
! ax-04 vec1
!-----------------------------------------------------------
      subroutine ax04_1(A, AX, X, NP, NP_ORG, LIST)
      implicit none

      ! argument
      integer*4 :: NP, NP_ORG
      real*4    :: A(30,NP)     ! matrix
      real*4    :: AX(NP)       ! result vector
      real*4    :: X(0:NP_ORG)    ! input vector
      integer*4 :: LIST(30,NP)  ! list

      ! local
      integer*4 :: IP, II, JJ, NN, IP2

      PROF_START("AX04_1")
!$omp parallel default(none) private(IP,II,JJ,IP2)
!$omp&shared(NP,LIST,AX,A,X)
!$omp do
      DO IP=1,NP
!ocl unroll('full')         
         DO II=1,30
            IP2=LIST(II,IP)
            AX(IP)=AX(IP)+A(II,IP)*X(IP2)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel
      PROF_STOP("AX04_1")

      return
      end subroutine ax04_1
      
!-----------------------------------------------------------
! ax-04 vec2
!-----------------------------------------------------------
      subroutine ax04_2(A, AX, X, NP, NP_ORG, LIST)
      implicit none

      ! argument
      integer*4 :: NP, NP_ORG
      real*4    :: A(30,NP)     ! matrix
      real*4    :: AX(2,NP)     ! result vector
      real*4    :: X(2,0:NP_ORG)  ! input vector
      integer*4 :: LIST(30,NP)  ! list

      ! local
      integer*4 :: IP, II, JJ, NN, IP2

      PROF_START("AX04_2")
!$omp parallel default(none) private(IP,II,JJ,IP2)
!$omp&shared(NP,LIST,AX,A,X)
!$omp do
      DO IP=1,NP
!ocl unroll('full')         
         DO II=1,30
            IP2=LIST(II,IP)
!ocl nounroll            
            DO JJ=1,2
               AX(JJ,IP)=AX(JJ,IP)+A(II,IP)*X(JJ,IP2)
            ENDDO
            
!           V=A(II,IP)
!           AX( 1,IP)=AX( 1,IP)+V*X( 1,IP2)
!           AX( 2,IP)=AX( 2,IP)+V*X( 2,IP2)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel
      PROF_STOP("AX04_2")

      return
      end subroutine ax04_2
      
!-----------------------------------------------------------
! ax-04 vec4
!-----------------------------------------------------------
      subroutine ax04_4(A, AX, X, NP, NP_ORG, LIST)
      implicit none

      ! argument
      integer*4 :: NP, NP_ORG
      real*4    :: A(30,NP)     ! matrix
      real*4    :: AX(4,NP)     ! result vector
      real*4    :: X(4,0:NP_ORG)  ! input vector
      integer*4 :: LIST(30,NP)  ! list

      ! local
      integer*4 :: IP, II, JJ, NN, IP2

      PROF_START("AX04_4")
!$omp parallel default(none) private(IP,II,JJ,IP2)
!$omp&shared(NP,LIST,AX,A,X)
!$omp do
      DO IP=1,NP
!ocl unroll('full')         
         DO II=1,30
            IP2=LIST(II,IP)
!ocl nounroll            
            DO JJ=1,4
               AX(JJ,IP)=AX(JJ,IP)+A(II,IP)*X(JJ,IP2)
            ENDDO
            
!           V=A(II,IP)
!           AX( 1,IP)=AX( 1,IP)+V*X( 1,IP2)
!           AX( 2,IP)=AX( 2,IP)+V*X( 2,IP2)
!           AX( 3,IP)=AX( 3,IP)+V*X( 3,IP2)
!           AX( 4,IP)=AX( 4,IP)+V*X( 4,IP2)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel
      PROF_STOP("AX04_4")

      return
      end subroutine ax04_4

!-----------------------------------------------------------
! ax-04 vec6
!-----------------------------------------------------------
      subroutine ax04_6(A, AX, X, NP, NP_ORG, LIST)
      implicit none

      ! argument
      integer*4 :: NP, NP_ORG
      real*4    :: A(30,NP)     ! matrix
      real*4    :: AX(6,NP)     ! result vector
      real*4    :: X(6,0:NP_ORG)  ! input vector
      integer*4 :: LIST(30,NP)  ! list

      ! local
      integer*4 :: IP, II, JJ, NN, IP2

      PROF_START("AX04_6")
!$omp parallel default(none) private(IP,II,JJ,IP2)
!$omp&shared(NP,LIST,AX,A,X)
!$omp do
      DO IP=1,NP
!ocl unroll('full')         
         DO II=1,30
            IP2=LIST(II,IP)
!ocl nounroll            
            DO JJ=1,6
               AX(JJ,IP)=AX(JJ,IP)+A(II,IP)*X(JJ,IP2)
            ENDDO
            
!           V=A(II,IP)
!           AX( 1,IP)=AX( 1,IP)+V*X( 1,IP2)
!           AX( 2,IP)=AX( 2,IP)+V*X( 2,IP2)
!           AX( 3,IP)=AX( 3,IP)+V*X( 3,IP2)
!           AX( 4,IP)=AX( 4,IP)+V*X( 4,IP2)
!           AX( 5,IP)=AX( 5,IP)+V*X( 5,IP2)
!           AX( 6,IP)=AX( 6,IP)+V*X( 6,IP2)
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel
      PROF_STOP("AX04_6")

      return
      end subroutine ax04_6

!-----------------------------------------------------------
! ax-04 vec8
!-----------------------------------------------------------
      subroutine ax04_8(A, AX, X, NP, NP_ORG, LIST)
      implicit none

      ! argument
      integer*4 :: NP, NP_ORG
      real*4    :: A(30,NP)     ! matrix
      real*4    :: AX(8,NP)     ! result vector
      real*4    :: X(8,0:NP_ORG)  ! input vector
      integer*4 :: LIST(30,NP)  ! list

      ! local
      integer*4 :: IP, II, JJ, NN, IP2

      PROF_START("AX04_8")
!$omp parallel default(none) private(IP,II,JJ,IP2)
!$omp&shared(NP,LIST,AX,A,X)
!$omp do
!ocl nounroll
!ocl swp
      DO IP=1,NP
!fj!ocl unroll('full')         
!fj         DO II=1,30
!fj            IP2=LIST(II,IP)
!ocl nounroll,loop_nofission
            DO JJ=1,8
               AX(JJ,IP  )=AX(JJ,IP  )+A( 1,IP  )*X(JJ,LIST( 1,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 2,IP  )*X(JJ,LIST( 2,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 3,IP  )*X(JJ,LIST( 3,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 4,IP  )*X(JJ,LIST( 4,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 5,IP  )*X(JJ,LIST( 5,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 6,IP  )*X(JJ,LIST( 6,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 7,IP  )*X(JJ,LIST( 7,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 8,IP  )*X(JJ,LIST( 8,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A( 9,IP  )*X(JJ,LIST( 9,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(10,IP  )*X(JJ,LIST(10,IP  ))

               AX(JJ,IP  )=AX(JJ,IP  )+A(11,IP  )*X(JJ,LIST(11,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(12,IP  )*X(JJ,LIST(12,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(13,IP  )*X(JJ,LIST(13,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(14,IP  )*X(JJ,LIST(14,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(15,IP  )*X(JJ,LIST(15,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(16,IP  )*X(JJ,LIST(16,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(17,IP  )*X(JJ,LIST(17,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(18,IP  )*X(JJ,LIST(18,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(19,IP  )*X(JJ,LIST(19,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(20,IP  )*X(JJ,LIST(20,IP  ))

               AX(JJ,IP  )=AX(JJ,IP  )+A(21,IP  )*X(JJ,LIST(21,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(22,IP  )*X(JJ,LIST(22,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(23,IP  )*X(JJ,LIST(23,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(24,IP  )*X(JJ,LIST(24,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(25,IP  )*X(JJ,LIST(25,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(26,IP  )*X(JJ,LIST(26,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(27,IP  )*X(JJ,LIST(27,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(28,IP  )*X(JJ,LIST(28,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(29,IP  )*X(JJ,LIST(29,IP  ))
               AX(JJ,IP  )=AX(JJ,IP  )+A(30,IP  )*X(JJ,LIST(30,IP  ))
!fj               AX(JJ,IP)=AX(JJ,IP)+A(II,IP)*X(JJ,IP2)
            ENDDO
            
!           V=A(II,IP)
!           AX( 1,IP)=AX( 1,IP)+V*X( 1,IP2)
!           AX( 2,IP)=AX( 2,IP)+V*X( 2,IP2)
!           AX( 3,IP)=AX( 3,IP)+V*X( 3,IP2)
!           AX( 4,IP)=AX( 4,IP)+V*X( 4,IP2)
!           AX( 5,IP)=AX( 5,IP)+V*X( 5,IP2)
!           AX( 6,IP)=AX( 6,IP)+V*X( 6,IP2)
!           AX( 7,IP)=AX( 7,IP)+V*X( 7,IP2)
!           AX( 8,IP)=AX( 8,IP)+V*X( 8,IP2)
!fj         ENDDO
      ENDDO
!$omp end do
!$omp end parallel
      PROF_STOP("AX04_8")

      return
      end subroutine ax04_8
