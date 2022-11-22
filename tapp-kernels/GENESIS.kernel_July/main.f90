!#include "profiler.h"
program nonbond_kernel
    use gparameter
    use module_pointers
    implicit none
    integer  :: step
    type(s_genesis_kernel_param):: gparam

    call read_data_file(gparam)

    do step = 1, 1000
      call kernel(gparam)
    end do

    call check_validation(gparam%MaxAtom,gparam%ncell,gparam%nthread,gparam%force)

end program nonbond_kernel


