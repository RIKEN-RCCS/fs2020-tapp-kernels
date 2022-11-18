module gparameter
  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  real(wp), public, parameter :: expected_result(1:3)=[-108.248734,282.639008,-239.026306]

  type,public  :: s_genesis_kernel_param

  integer                     :: maxcell_near
  integer                     :: maxcell
  integer                     :: ncell
  integer                     :: ncell_local
  integer                     :: num_atom_cls
  integer                     :: MaxAtom
  integer                     :: nthread
  integer                     :: ntable
  real(wp)                    :: density
  real(wp)                    :: cutoff
  integer, allocatable        :: natom(:)
  integer, allocatable        :: atmcls(:,:)
  integer, allocatable        :: cell_pairlist1(:,:)
  integer, allocatable        :: virial_check(:,:)
  integer, allocatable        :: nb15_cell(:)
  integer, allocatable        :: nb15_list(:,:)
  integer(1), allocatable     :: exclusion_mask1(:,:,:)
  integer(1), allocatable     :: exclusion_mask(:,:,:)
  real(wp) , allocatable      :: system_size(:)
  real(wp) , allocatable      :: cell_move(:,:,:)
  real(wp) , allocatable      :: charge(:,:)
  real(wp) , allocatable      :: coord(:,:,:)
  real(wp) , allocatable      :: coord_pbc(:,:,:)
  real(wp) , allocatable      :: trans1(:,:,:)
  real(wp) , allocatable      :: table_grad(:)
  real(wp) , allocatable      :: force(:,:,:,:)
  real(wp) , allocatable      :: virial(:,:)
  real(wp) , allocatable      :: nonb_lj12(:,:)
  real(wp) , allocatable      :: nonb_lj6(:,:)

  end type s_genesis_kernel_param

end module gparameter

module module_pointers
use gparameter
contains

subroutine read_data_file (gparam)

  type(s_genesis_kernel_param), intent(inout) :: gparam

  character*20 :: filename
  integer                     :: omp_get_num_threads
  integer                     :: maxcell_near
  integer                     :: maxcell
  integer                     :: ncell
  integer                     :: ncell_local
  integer                     :: num_atom_cls
  integer                     :: MaxAtom
  integer                     :: nthread
  integer                     :: ntable
  real(wp)                    :: density
  real(wp)                    :: cutoff
  integer   :: alloc_stat
  integer :: iunit, is_ok

  !$omp parallel 
  gparam%nthread = omp_get_num_threads()
  !$omp end parallel

  iunit=77
  write(filename,'(a)') "data_kernel_July"

  open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
  if (is_ok.ne.0) then
  write(*,'(a,a)') "*** Error. failed to open file: ", filename
  endif

  call sub_i4_data_read(iunit, gparam%maxcell, 1 )
  call sub_i4_data_read(iunit, gparam%maxcell_near, 1 )
  call sub_i4_data_read(iunit, gparam%ncell, 1 )
  call sub_i4_data_read(iunit, gparam%ncell_local, 1 )
  call sub_i4_data_read(iunit, gparam%num_atom_cls, 1 )
  call sub_i4_data_read(iunit, gparam%MaxAtom, 1 )
  call sub_i4_data_read(iunit, gparam%ntable, 1 )

  maxcell      = gparam%maxcell
  maxcell_near = gparam%maxcell_near
  ncell        = gparam%ncell
  ncell_local  = gparam%ncell_local
  num_atom_cls = gparam%num_atom_cls
  MaxAtom      = gparam%MaxAtom
  nthread      = gparam%nthread
  ntable       = gparam%ntable

  allocate(gparam%charge(1:MaxAtom,1:ncell),                    &
           gparam%atmcls(1:MaxAtom,1:ncell),                    &
           gparam%nonb_lj12(1:num_atom_cls, 1:num_atom_cls),    &
           gparam%nonb_lj6(1:num_atom_cls, 1:num_atom_cls),     &
           gparam%coord(1:3,1:MaxAtom,1:ncell),                 &
           gparam%coord_pbc(1:MaxAtom,1:3,1:ncell),             &
           gparam%trans1(1:3,1:MaxAtom,1:ncell),                &
           gparam%table_grad(1:ntable*6),                       &
           gparam%force(1:MaxAtom,1:3,1:ncell,1:nthread),       &
           gparam%virial(1:3,1:maxcell),                        &
           gparam%natom(1:ncell),                               &
           gparam%cell_pairlist1(1:2,1:maxcell),                &
           gparam%nb15_cell(1:maxcell),                         &
           gparam%nb15_list(1:2*MaxAtom,1:maxcell),             &
           gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local), &
           gparam%exclusion_mask(1:MaxAtom,1:MaxAtom,1:maxcell_near), &
           gparam%virial_check(1:ncell,1:ncell),                &
           gparam%cell_move(1:3,1:ncell,1:ncell),               &
           gparam%system_size(1:3),                             &
           stat = alloc_stat)

  call sub_i4_data_read(iunit, gparam%natom, (ncell) )
  call sub_i4_data_read(iunit, gparam%cell_pairlist1, (2*maxcell) )
  call sub_i4_data_read(iunit, gparam%nb15_cell, (maxcell) )
  call sub_i4_data_read(iunit, gparam%nb15_list, (2*MaxAtom*maxcell) )
  call sub_i1_data_read(iunit, gparam%exclusion_mask1, (MaxAtom*MaxAtom*ncell_local) )
  call sub_i1_data_read(iunit, gparam%exclusion_mask, (MaxAtom*MaxAtom*maxcell_near) )
  call sub_i4_data_read(iunit, gparam%atmcls, (MaxAtom*ncell) )
  call sub_r4_data_read(iunit, gparam%nonb_lj12, (num_atom_cls*num_atom_cls) )
  call sub_r4_data_read(iunit, gparam%nonb_lj6, (num_atom_cls*num_atom_cls) )
  call sub_i4_data_read(iunit, gparam%virial_check, (ncell*ncell) )
  call sub_r4_data_read(iunit, gparam%coord, (3*MaxAtom*ncell) )
  call sub_r4_data_read(iunit, gparam%trans1, (3*MaxAtom*ncell) )
  call sub_r4_data_read(iunit, gparam%charge, (MaxAtom*ncell) )
  call sub_r4_data_read(iunit, gparam%cell_move, (3*ncell*ncell) )
  call sub_r4_data_read(iunit, gparam%system_size, (3) )
  call sub_r4_data_read(iunit, gparam%table_grad, (ntable*6) )
  call sub_r4_data_read(iunit, gparam%density, 1)
  call sub_r4_data_read(iunit, gparam%cutoff, 1)
close(iunit)

end subroutine

subroutine check_validation (MaxAtom, ncell, nthread, force)
    use gparameter

    real(wp) :: force(1:MaxAtom,1:3,1:ncell,1:nthread)
    real(wp) :: val, check_result

    do i=1,3
      val=0.0_wp
      do l=1,nthread
        do k=1,ncell
          do j=1,MaxAtom, 5 ! since sum of all forces should zero. 5 is required.
            val=val+force(j,i,k,l)
          enddo
        enddo
      enddo

      check_result = val / expected_result(i)

      if ( 0.999 < check_result .and. check_result < 1.001 ) then
          write(*,*) "val=", val, "    expected_result=", expected_result(i)
          write(*,*) "the computed result seems to be OK."
      else
          write(*,*) "val=", val, "    expected_result=", expected_result(i)
          write(*,*) "the computed result is not close enough to the expected value."
      endif
    enddo

end subroutine

end module 
