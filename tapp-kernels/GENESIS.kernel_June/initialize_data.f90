module gparameter
  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp
  real(wp), public, parameter :: expected_result(1:3)=[-249.517136,269.371796,-305.023773]

  type,public  :: s_genesis_kernel_param

  integer                     :: maxcell
  integer                     :: ncell
  integer                     :: ncell_local
  integer                     :: MaxAtom
  integer                     :: MaxAtomCls
  integer                     :: MaxNb15
  integer                     :: nthread
  integer                     :: ntable
  real(wp)                    :: density
  real(wp)                    :: cutoff
  real(wp)                    :: cutoff2
  real(wp)                    :: inv_MaxAtom
  integer, allocatable        :: natom(:)
  integer, allocatable        :: atmcls(:,:)
  integer, allocatable        :: num_nb15_calc(:,:)
  integer, allocatable        :: nb15_calc_list(:,:,:)
  real(wp) , allocatable      :: charge(:,:)
  real(wp) , allocatable      :: coord(:,:,:)
  real(wp) , allocatable      :: coord_pbc(:,:,:)
  real(wp) , allocatable      :: trans1(:,:,:)
  real(wp) , allocatable      :: table_grad(:)
  real(wp) , allocatable      :: force(:,:,:,:)
  real(wp) , allocatable      :: lj_coef(:,:)

  end type s_genesis_kernel_param

end module gparameter

module module_pointers
use gparameter
contains

subroutine read_data_file (gparam)

  type(s_genesis_kernel_param), intent(inout) :: gparam

    integer                  :: omp_get_num_threads
    integer                  :: maxcell,ncell,ncell_local
    integer                  :: alloc_stat
    integer                  :: nthread
    integer                  :: MaxAtom,MaxAtomCls,MaxNb15
    real(wp)                  :: density, cutoff


    character*20 :: filename
    integer :: iunit, is_ok

  !$omp parallel 
  gparam%nthread = omp_get_num_threads()
  !$omp end parallel

    iunit=77
    write(filename,'(a)') "data_kernel_June"

    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_read(iunit,  gparam%maxcell     , 1)
    call sub_i4_data_read(iunit,  gparam%ncell       , 1)
    call sub_i4_data_read(iunit,  gparam%ncell_local , 1)
    call sub_i4_data_read(iunit,  gparam%MaxAtom     , 1)
    call sub_i4_data_read(iunit,  gparam%MaxAtomCls  , 1)
    call sub_i4_data_read(iunit,  gparam%MaxNb15     , 1)
    call sub_i4_data_read(iunit,  gparam%ntable     , 1)
    call sub_r4_data_read(iunit,  gparam%inv_MaxAtom  , 1)
    maxcell = gparam%maxcell
    ncell = gparam%ncell
    ncell_local = gparam%ncell_local
    MaxAtom = gparam%MaxAtom
    MaxAtomCls = gparam%MaxAtomCls
    MaxNb15 = gparam%MaxNb15
    ntable = gparam%ntable
    nthread = gparam%nthread

    allocate(gparam%charge(1:MaxAtom,1:ncell),                &
             gparam%atmcls(1:MaxAtom,1:ncell),                &
             gparam%lj_coef(1:2,1:MaxAtomCls),                &
             gparam%num_nb15_calc(1:MaxAtom,1:maxcell),       &
             gparam%nb15_calc_list(1:MaxNb15,1:MaxAtom,1:ncell),  &
             gparam%coord(1:3,1:MaxAtom,1:ncell),             &
             gparam%coord_pbc(1:3,1:MaxAtom,1:ncell),         &
             gparam%trans1(1:3,1:MaxAtom,1:ncell),            &
             gparam%table_grad(1:ntable*6),                   &
             gparam%force(1:3,1:MaxAtom,1:ncell,1:nthread),   &
             gparam%natom(1:ncell),                           &
             stat = alloc_stat)

    call sub_i4_data_read(iunit, gparam%natom        , (ncell) )
    call sub_r4_data_read(iunit, gparam%charge       , (MaxAtom*ncell) )
    call sub_i4_data_read(iunit, gparam%atmcls       , (MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%lj_coef      , (2*MaxAtomCls) )
    call sub_i4_data_read(iunit, gparam%num_nb15_calc  , (MaxAtom*maxcell) )
    call sub_i4_data_read(iunit, gparam%nb15_calc_list , (MaxNb15*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%coord          , (3*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%trans1         , (3*MaxAtom*ncell) )
    call sub_r4_data_read(iunit, gparam%table_grad   , (ntable*6) )
    call sub_r4_data_read(iunit, gparam%density      , 1)
    call sub_r4_data_read(iunit, gparam%cutoff       , 1)
    call sub_r4_data_read(iunit, gparam%cutoff2      , 1)

    close(iunit)

    write(*,*) "maxcell=", maxcell, "   loc(maxcell)=", loc(maxcell)

end subroutine


subroutine check_validation (MaxAtom,ncell,nthread,force)
    use gparameter
    real(wp) :: force(1:3,1:MaxAtom,1:ncell,1:nthread)
    real(wp) :: val(1:3),  check_result
    integer  :: j, k, l

    val(1:3)=0.0d0
    do l=1,nthread
      do k=1,ncell
        do j=1,MaxAtom, 5 ! since sum of all forces should zero. 5 is required.
          do i=1,3
            val(i)=val(i)+real(force(i,j,k,l),wp)
          enddo
        enddo
      enddo
    enddo

    do i=1,3
      check_result = val(i) / expected_result(i)

      if ( 0.999 < check_result .and. check_result < 1.001 ) then
          write(*,*) "val=", val(i), "    expected_result=", expected_result(i)
          write(*,*) "the computed result seems to be OK."
      else
          write(*,*) "val=", val(i), "    expected_result=", expected_result(i)
          write(*,*) "the computed result is not close enough to the expected value."
      endif
    end do

end subroutine

end module

