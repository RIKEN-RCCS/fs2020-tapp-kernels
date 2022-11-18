!#include "profiler.h"

program nonbond_kernel
    use gparameter
    use module_pointers
    implicit none
    integer  :: step
    type(s_genesis_kernel_param):: gparam

    call read_data_file(gparam)

!    PROF_INIT
!    PROF_START_ALL

!    do step = 1, 1000
    do step = 1, 1
      call kernel(gparam)
    end do

!    PROF_STOP_ALL
!    PROF_FINALIZE

    call result_validation(gparam%MaxAtom, gparam%maxcell, gparam%num_nb15_calc)

end program nonbond_kernel

