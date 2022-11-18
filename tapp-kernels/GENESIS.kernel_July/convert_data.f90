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
    integer                     :: num_atom_cls
    integer                     :: MaxAtom
    integer                     :: nthread
    integer                     :: ntable
    integer                     :: nxyz
    real(wp)                    :: density
    real(wp)                    :: cutoff
    real(wp)                    :: cutoff2
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
                 MaxAtom,  nthread, omp_get_num_threads, nxyz, ntable, &
                 num_atom_cls
    integer :: iunit, is_ok
    character(35) :: filename

!    !$omp parallel shared(nthread)
!    gparam%nthread = omp_get_num_threads()
!    !$omp end parallel
    gparam%nthread =1
    

    open(1001, file='kernel_july', status='old',  &
         form='unformatted', err=900)

    read(1001) gparam%maxcell_near, gparam%maxcell, gparam%ncell, &
               gparam%ncell_local, gparam%num_atom_cls,           &
               gparam%MaxAtom
    read(1001) gparam%ntable

    alloc_stat = 0
    nthread      = gparam%nthread
    maxcell      = gparam%maxcell
    maxcell_near = gparam%maxcell_near
    ncell        = gparam%ncell
    ncell_local  = gparam%ncell_local
    num_atom_cls = gparam%num_atom_cls
    MaxAtom      = gparam%MaxAtom
    nthread      = gparam%nthread
    ntable       = gparam%ntable
    gparam%nxyz=3
    nxyz         = gparam%nxyz

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

    if (alloc_stat /= 0) then
       write(0,*) 'error: allocation'
       stop
    endif

    read(1001) gparam%natom(1:ncell)
    read(1001) gparam%cell_pairlist1(1:2,1:maxcell)
    read(1001) gparam%nb15_cell(1:maxcell)
    read(1001) gparam%nb15_list(1:2*MaxAtom,1:maxcell)
    read(1001) gparam%exclusion_mask1(1:MaxAtom,1:MaxAtom,1:ncell_local)
    read(1001) gparam%exclusion_mask(1:MaxAtom,1:MaxAtom,1:maxcell_near)
    read(1001) gparam%virial_check(1:ncell,1:ncell)
    read(1001) gparam%cell_move(1:3,1:ncell,1:ncell)
    read(1001) gparam%charge(1:MaxAtom,1:ncell)
    read(1001) gparam%atmcls(1:MaxAtom,1:ncell)
    read(1001) gparam%nonb_lj12(1:num_atom_cls,1:num_atom_cls)
    read(1001) gparam%nonb_lj6(1:num_atom_cls,1:num_atom_cls)
    read(1001) gparam%coord(1:3,1:MaxAtom,1:ncell)
    read(1001) gparam%trans1(1:3,1:MaxAtom,1:ncell)
    read(1001) gparam%table_grad(1:ntable*6)
    read(1001) gparam%density, gparam%cutoff, gparam%cutoff2
    read(1001) gparam%system_size(1:3)

    gparam%coord_pbc(1:3,1:MaxAtom,1:ncell)=0.0_wp

    close(1001)

    iunit=77
    write(filename,'(a)') "data_kernel_July"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_write(iunit, gparam%maxcell,          1 , "maxcell             ")
    call sub_i4_data_write(iunit, gparam%maxcell_near,     1 , "maxcell_near        ")
    call sub_i4_data_write(iunit, gparam%ncell,            1  ,"ncell               ")
    call sub_i4_data_write(iunit, gparam%ncell_local,      1  ,"ncell_local         ")
    call sub_i4_data_write(iunit, gparam%num_atom_cls,      1 ,"num_atom_cls        ")
    call sub_i4_data_write(iunit, gparam%MaxAtom,          1 , "MaxAtom             ")
    call sub_i4_data_write(iunit, gparam%ntable,          1 ,  "ntable              ")

    call sub_i4_data_write(iunit, gparam%natom, (ncell),                "natom               " )
    call sub_i4_data_write(iunit, gparam%cell_pairlist1, (2*maxcell),   "cell_pairlist1      " )
    call sub_i4_data_write(iunit, gparam%nb15_cell, (maxcell),          "nb15_cell           " )
    call sub_i4_data_write(iunit, gparam%nb15_list, (2*MaxAtom*maxcell),"nb15_list           " )
    call sub_i1_data_write(iunit, gparam%exclusion_mask1, (MaxAtom*MaxAtom*ncell_local),"exclusion_mask1     " )
    call sub_i1_data_write(iunit, gparam%exclusion_mask, (MaxAtom*MaxAtom*maxcell_near) ,"exclusion_mask     ")
    call sub_i4_data_write(iunit, gparam%atmcls, (MaxAtom*ncell),       "atomcls             " )
    call sub_r4_data_write(iunit, gparam%nonb_lj12, (num_atom_cls*num_atom_cls),"nonb_lj12           " )
    call sub_r4_data_write(iunit, gparam%nonb_lj6, (num_atom_cls*num_atom_cls),"nonb_lj6            " )
    call sub_i4_data_write(iunit, gparam%virial_check, (ncell*ncell),   "virial_check        " )
    call sub_r4_data_write(iunit, gparam%coord, (3*MaxAtom*ncell),      "coord               " )
    call sub_r4_data_write(iunit, gparam%trans1, (3*MaxAtom*ncell),     "trans1              " )
    call sub_r4_data_write(iunit, gparam%charge, (MaxAtom*ncell),       "charge              " )
    call sub_r4_data_write(iunit, gparam%cell_move, (3*ncell*ncell),    "cell_move           " )
    call sub_r4_data_write(iunit, gparam%system_size, (3),              "system_size         " )
    call sub_r4_data_write(iunit, gparam%table_grad, (ntable*6),        "table_grad          " )
    call sub_r4_data_write(iunit, gparam%density, 1,                    "density             ")
    call sub_r4_data_write(iunit, gparam%cutoff, 1,                     "cutoff              ")

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
               gparam%exclusion_mask1,      &
               gparam%exclusion_mask,       &
               gparam%nb15_cell,            &
               gparam%nb15_list,            &
               stat = dealloc_stat)
    if (dealloc_stat /= 0) then
       write(0,*) 'error: deallocation'
       stop
    endif
    return
  end subroutine reset_parameters

end program nonbond_kernel

