module conv

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp

  type,public  :: s_genesis_kernel_param

    integer                     :: maxcell_near
    integer                     :: maxcell
    integer                     :: ncell
    integer                     :: ncell_local
    integer                     :: MaxAtom
    integer                     :: MaxNb15
    integer                     :: nthread
    integer                     :: nxyz
    real(wp)                    :: pairdist2
    integer, allocatable        :: natom(:)
    integer, allocatable        :: cell_pairlist1(:,:)
    integer, allocatable        :: num_nb15_calc(:,:)
    integer, allocatable        :: nb15_calc_list(:,:,:)
    integer(1), allocatable     :: exclusion_mask1(:,:,:)
    integer(1), allocatable     :: exclusion_mask(:,:,:)
    real(wp) , allocatable      :: coord(:,:,:)
    real(wp) , allocatable      :: coord_pbc(:,:,:)
    real(wp) , allocatable      :: trans1(:,:,:)

  end type s_genesis_kernel_param

end module conv

program nonbond_kernel

  use conv

  implicit none

  type(s_genesis_kernel_param):: gparam

  call readset_parameters(gparam)
  call reset_parameters(gparam)

stop

contains

  subroutine readset_parameters(gparam)
    type(s_genesis_kernel_param),   intent(out)    :: gparam

    integer   :: alloc_stat
    integer   :: maxcell, maxcell_near, ncell, ncell_local, &
                 MaxAtom, MaxNb15, nthread, omp_get_num_threads, nxyz
    integer :: iunit, is_ok
    character(35) :: filename

!    !$omp parallel shared(nthread)
!    gparam%nthread = omp_get_num_threads()
!    !$omp end parallel
    gparam%nthread =1
    

!    open(1001, file='pairlist_june', status='old',  &
!         form='unformatted', err=900)
    open(1001, file='pairlist_june_jj_20170823', status='old',  &
         form='unformatted', err=900)

    read(1001) gparam%maxcell_near, gparam%maxcell, gparam%ncell,  &
               gparam%ncell_local, gparam%MaxAtom, gparam%MaxNb15
!    write(6,*)       ' maxcell_near maxcell    ncell  ncell_local   MaxAtom'// &
!                     ' MaxNb15 '
!    write(6,'(10I10)') gparam%maxcell_near, gparam%maxcell, gparam%ncell,  &
!                       gparam%ncell_local, gparam%MaxAtom, gparam%MaxNb15

    alloc_stat = 0
    maxcell      = gparam%maxcell
    maxcell_near = gparam%maxcell_near
    ncell        = gparam%ncell
    ncell_local  = gparam%ncell_local
    MaxAtom      = gparam%MaxAtom
    MaxNb15      = gparam%MaxNb15
    nthread      = gparam%nthread
    gparam%nxyz=3
    nxyz         = gparam%nxyz

    allocate(                                                     &
             gparam%coord(1:3,1:MaxAtom,1:ncell),                 &
             gparam%coord_pbc(1:3,1:MaxAtom,1:ncell),             &
             gparam%trans1(1:3,1:MaxAtom,1:ncell),                &
             gparam%natom(1:ncell),                               &
             gparam%cell_pairlist1(1:2,1:maxcell),                &
             gparam%num_nb15_calc(1:MaxAtom,1:ncell),             &
             gparam%nb15_calc_list(1:MaxNb15,1:MaxAtom,1:ncell),  &
             gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local),  & 
             gparam%exclusion_mask (1:MaxAtom,1:MaxAtom,1:maxcell_near), &      
             stat = alloc_stat)

    if (alloc_stat /= 0) then
       write(0,*) 'error: allocation'
       stop
    endif

    read(1001) gparam%natom(1:ncell)
    read(1001) gparam%cell_pairlist1(1:2,1:maxcell)
    read(1001) gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local)
    read(1001) gparam%exclusion_mask (1:MaxAtom,1:MaxAtom,1:maxcell_near)
    read(1001) gparam%coord(1:3,1:MaxAtom,1:ncell)
    read(1001) gparam%trans1(1:3,1:MaxAtom,1:ncell)
    read(1001) gparam%pairdist2

    gparam%coord_pbc(1:3,1:MaxAtom,1:ncell)=0.0_wp

    close(1001)

    iunit=77
    write(filename,'(a)') "data_pairlist_june_JJ_20170823-asci"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_write(iunit, gparam%maxcell_near,     1 , "maxcell_near        ")
    call sub_i4_data_write(iunit, gparam%maxcell,          1 , "maxcell             ")
    call sub_i4_data_write(iunit, gparam%ncell,            1 , "ncell               ")
    call sub_i4_data_write(iunit, gparam%ncell_local,      1 , "ncell_local         ")
    call sub_i4_data_write(iunit, gparam%MaxAtom,          1 , "MaxAtom             ")
    call sub_i4_data_write(iunit, gparam%MaxNb15,          1 , "MaxNb15             ")

    call sub_i4_data_write(iunit, gparam%natom,           (ncell),     "natom               ")
    call sub_i4_data_write(iunit, gparam%cell_pairlist1,  (2*maxcell), "cell_pairlist1      " )
    call sub_i1_data_write(iunit, gparam%exclusion_mask1, (MaxAtom*MaxAtom*ncell_local),"exclusion_mask1     " )
    call sub_i1_data_write(iunit, gparam%exclusion_mask,  (MaxAtom*MaxAtom*maxcell_near),"exclusion_mask      " )
    call sub_r4_data_write(iunit, gparam%coord,           (nxyz*MaxAtom*ncell) ,"coord               ")
    call sub_r4_data_write(iunit, gparam%trans1,          (nxyz*MaxAtom*ncell) ,"trans               ")
    call sub_r4_data_write(iunit, gparam%pairdist2,        1, "pairdist2           ")
    close(iunit)

    return
900 write(0, *) 'error: read'
    stop

  end subroutine readset_parameters

  subroutine reset_parameters(gparam)
    type(s_genesis_kernel_param),   intent(inout)    :: gparam

    integer   :: dealloc_stat

    dealloc_stat = 0
    deallocate(                             &
               gparam%coord,                &
               gparam%coord_pbc,            &
               gparam%trans1,               &
               gparam%natom,                &
               gparam%cell_pairlist1,       &
               gparam%num_nb15_calc,        &
               gparam%nb15_calc_list,       &
               stat = dealloc_stat)
    if (dealloc_stat /= 0) then
       write(0,*) 'error: deallocation'
       stop
    endif
    return
  end subroutine reset_parameters

end program nonbond_kernel

