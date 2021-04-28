#include "profiler.h"

program nonbond_kernel
    use gparameter
    use module_pointers
    implicit none
    integer  :: step
    type(s_genesis_kernel_param):: gparam

    call read_data_file()
    !cx call check_data_file()

    call readset_parameters(gparam)

    PROF_INIT
    PROF_START_ALL

    do step = 1, 1000
        call kernel(gparam)
    end do

    PROF_STOP_ALL
    PROF_FINALIZE

    call result_validation()

end program nonbond_kernel

