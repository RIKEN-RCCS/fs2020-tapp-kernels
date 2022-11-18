module gparameter

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  real(wp), public, parameter :: expected_result=   8179114.00

  type,public  :: s_genesis_kernel_param

    integer                     :: maxcell_near
    integer                     :: maxcell
    integer                     :: ncell
    integer                     :: ncell_local
    integer                     :: MaxAtom
    integer                     :: nthread
    real(wp)                    :: pairdist2
    real(wp)                    :: system_size(3)
    integer, allocatable        :: natom(:)
    integer, allocatable        :: cell_pairlist1(:,:)
    integer, allocatable        :: nb15_cell(:)
    integer, allocatable        :: nb15_list(:,:)
    integer(1), allocatable     :: exclusion_mask1(:,:,:)
    integer(1), allocatable     :: exclusion_mask(:,:,:)
    real(wp) , allocatable      :: cell_move(:,:,:)
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
                 MaxAtom, nthread, omp_get_num_threads


    character(20) :: filename
    integer :: iunit, is_ok

    iunit=77
    write(filename,'(a)') "data_pairlist_July"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif
    !$omp parallel shared(nthread)
    gparam%nthread = omp_get_num_threads()
    !$omp end parallel
    call sub_i4_data_read(iunit, gparam%maxcell_near,     1 )
    call sub_i4_data_read(iunit, gparam%maxcell,          1 )
    call sub_i4_data_read(iunit, gparam%ncell,            1 )
    call sub_i4_data_read(iunit, gparam%ncell_local,      1 )
    call sub_i4_data_read(iunit, gparam%MaxAtom,          1 )

    alloc_stat = 0
    maxcell      = gparam%maxcell
    maxcell_near = gparam%maxcell_near
    ncell        = gparam%ncell
    ncell_local  = gparam%ncell_local
    MaxAtom      = gparam%MaxAtom
    nthread      = gparam%nthread
    nxyz         = 3

    allocate(                                                     &
             gparam%coord(1:3,1:MaxAtom,1:ncell),                 &
             gparam%coord_pbc_x(1:MaxAtom,1:ncell),             &
             gparam%coord_pbc_y(1:MaxAtom,1:ncell),             &
             gparam%coord_pbc_z(1:MaxAtom,1:ncell),             &
             gparam%trans1(1:3,1:MaxAtom,1:ncell),                &
             gparam%natom(1:ncell),                               &
             gparam%cell_pairlist1(1:2,1:maxcell),                &
             gparam%nb15_cell(1:maxcell),                         &
             gparam%nb15_list(1:2*MaxAtom,1:maxcell),            &
             gparam%cell_move(1:3,1:ncell,1:ncell),               &
             gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local),  &
             gparam%exclusion_mask (1:MaxAtom,1:MaxAtom,1:maxcell_near), &
             stat = alloc_stat)

    call sub_i4_data_read(iunit, gparam%natom,           (ncell) )
    call sub_i4_data_read(iunit, gparam%cell_pairlist1,  (2*maxcell) )
    call sub_i1_data_read(iunit, gparam%exclusion_mask1, (MaxAtom*MaxAtom*ncell_local) )
    call sub_i1_data_read(iunit, gparam%exclusion_mask,  (MaxAtom*MaxAtom*maxcell_near) )
    call sub_r4_data_read(iunit, gparam%cell_move,       (nxyz*ncell*ncell) )
    call sub_r4_data_read(iunit, gparam%coord,           (nxyz*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%trans1,          (nxyz*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%pairdist2,        1)
    call sub_r4_data_read(iunit, gparam%system_size,      3)

end subroutine

end module module_pointers

subroutine result_validation (MaxAtom, maxcell, nb15_list)
    use gparameter
 
    integer   :: nb15_list(1:2*MaxAtom,1:maxcell)
    integer   :: i, j
    real(dp)   :: val,  check_result
 
    val=0.d0
    do j=1,maxcell
      do i=1,2*MaxAtom
        val=val+real(nb15_list(i,j),dp)
      enddo
    enddo

!cx    call report_validation(val,814655.0_8,1d-1)
 
    check_result = val / expected_result
 
    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif
 
end subroutine

