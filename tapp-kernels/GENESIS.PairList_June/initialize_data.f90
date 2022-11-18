!#include "profiler.h"
module gparameter

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  integer,  public, parameter :: nxyz         = 3
  integer,  public, parameter :: ncq         = 4
  real(wp), public, parameter :: expected_result=   3057402.00


  type,public  :: s_genesis_kernel_param

    integer                     :: maxcell_near
    integer                     :: maxcell
    integer                     :: ncell
    integer                     :: ncell_local
    integer                     :: MaxAtom
    integer                     :: MaxNb15
    integer                     :: nthread
    real(wp)                    :: pairdist2
    integer, allocatable        :: natom(:)
    integer, allocatable        :: cell_pairlist1(:,:)
    integer, allocatable        :: num_nb15_calc(:,:)
    integer, allocatable        :: nb15_calc_list(:,:,:)
    integer(1), allocatable     :: exclusion_mask1(:,:,:)
    integer(1), allocatable     :: exclusion_mask(:,:,:)
    real(wp) , allocatable      :: coord(:,:,:)
    real(wp) , allocatable      :: coord_pbc_x(:,:)
    real(wp) , allocatable      :: coord_pbc_y(:,:)
    real(wp) , allocatable      :: coord_pbc_z(:,:)
    real(wp) , allocatable      :: trans1(:,:,:)

  end type s_genesis_kernel_param

end module gparameter

module module_pointers
use gparameter
contains

subroutine read_data_file (gparam)

    type(s_genesis_kernel_param), intent(inout) :: gparam
    integer   :: alloc_stat
    integer   :: maxcell, maxcell_near, ncell, ncell_local, &
                 MaxAtom, MaxNb15, nthread, omp_get_num_threads

    character(20) :: filename
    integer :: iunit, is_ok

    !$omp parallel shared(nthread)
    gparam%nthread = omp_get_num_threads()
    !$omp end parallel
    iunit=77
    write(filename,'(a)') "data_pairlist_June"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

!    read(1001) gparam%maxcell_near, gparam%maxcell, gparam%ncell,  &
!               gparam%ncell_local, gparam%MaxAtom, gparam%MaxNb15
    call sub_i4_data_read(iunit, gparam%maxcell_near,     1 )
    call sub_i4_data_read(iunit, gparam%maxcell,          1 )
    call sub_i4_data_read(iunit, gparam%ncell,            1 )
    call sub_i4_data_read(iunit, gparam%ncell_local,      1 )
    call sub_i4_data_read(iunit, gparam%MaxAtom,          1 )
    call sub_i4_data_read(iunit, gparam%MaxNb15,          1 )

    alloc_stat = 0
    maxcell      = gparam%maxcell
    maxcell_near = gparam%maxcell_near
    ncell        = gparam%ncell
    ncell_local  = gparam%ncell_local
    MaxAtom      = gparam%MaxAtom
    MaxNb15      = gparam%MaxNb15
    nthread      = gparam%nthread

    allocate(                                                     &
             gparam%coord(1:3,1:MaxAtom,1:ncell),                 &
             gparam%coord_pbc_x(1:MaxAtom,1:ncell),             &
             gparam%coord_pbc_y(1:MaxAtom,1:ncell),             &
             gparam%coord_pbc_z(1:MaxAtom,1:ncell),             &
             gparam%trans1(1:3,1:MaxAtom,1:ncell),                &
             gparam%natom(1:ncell),                               &
             gparam%cell_pairlist1(1:2,1:maxcell),                &
             gparam%num_nb15_calc(1:MaxAtom,1:maxcell),             &
             gparam%nb15_calc_list(1:MaxNb15,1:MaxAtom,1:maxcell),  &
             gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local),  &
             gparam%exclusion_mask (1:MaxAtom,1:MaxAtom,1:maxcell_near), &
             stat = alloc_stat)

    if (alloc_stat /= 0) then
       write(0,*) 'error: allocation'
       stop
    endif

    call sub_i4_data_read(iunit, gparam%natom,           (ncell) )
    call sub_i4_data_read(iunit, gparam%cell_pairlist1,  (2*maxcell) )
    call sub_i1_data_read(iunit, gparam%exclusion_mask1, (MaxAtom*MaxAtom*ncell_local) )
    call sub_i1_data_read(iunit, gparam%exclusion_mask,  (MaxAtom*MaxAtom*maxcell_near) )
    call sub_r4_data_read(iunit, gparam%coord,           (nxyz*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%trans1,          (nxyz*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%pairdist2,        1)

end subroutine

end module module_pointers


subroutine result_validation (MaxAtom, maxcell, c_num_nb15_calc)
    use gparameter

    integer        :: c_num_nb15_calc(1:MaxAtom,1:maxcell)
    integer        :: i1,i2
    real(dp)       :: val, check_result

    val = 0
    do i2 = 1, maxcell
     do i1 = 1, MaxAtom
        val = val+real(c_num_nb15_calc(i1,i2),dp)
     enddo
    enddo

!cx call report_validation(dble(val), 1.5581600d+05, 0.000001)

    check_result = val / expected_result

    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif

end subroutine

