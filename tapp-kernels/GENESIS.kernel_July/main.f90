#include "profiler.h"
program nonbond_kernel
    use module_pointers

    call read_data_file()

    do step = 1, 1000
    call kernel()
    end do

    call check_force_validation()

end program nonbond_kernel


