!#include "profiler.h"
program nonbond_kernel
    use gparameter
    use module_pointers
    implicit none
    integer  :: step
    type(s_genesis_kernel_param):: gparam

    call read_data_file(gparam)

    do step = 1, 1
      call kernel(gparam)
    end do

    call result_validation(gparam%MaxAtom, gparam%maxcell, gparam%nb15_list)

end program nonbond_kernel


