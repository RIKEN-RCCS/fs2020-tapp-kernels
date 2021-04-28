module module_pointers
contains

subroutine write_data_file ()

    ! local variables
    real(4)                  :: dij1,dij2,dij3, rij2
    real(4)                  :: R, lj12, lj6
    real(4)                  :: term_lj12, term_lj6, term_elec
    real(4)                  :: cutoff2, grad_coef
    real(4)                  :: work(1:3)
    real(4)                  :: rtmp(1:3), qtmp, jqtmp
    real(4)                  :: force_local(3)
    real(4)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                  :: i, ix, iy, j, k, ij, L, L1, ii
    integer                  :: id, ik, omp_get_thread_num
    integer                  :: iatmcls,jatmcls
    integer                  :: maxcell,ncell,ncell_local
    integer                  :: nthread
    integer                  :: MaxAtom,MaxAtomCls,MaxNb15
    real(4)                  :: density, cutoff
    real(4)                  :: inv_MaxAtom
    real(8)                  :: Val

    integer                  :: natom(1:256,1:4)
    integer                  :: atmcls(1:181,1:256,1:4)
    integer                  :: num_nb15_calc(1:181,1:4683,1:4)
    integer                  :: nb15_calc_list(1:950,1:100,1:256,1:4)
    real(4)                  :: charge(1:181,1:256,1:4)
    real(4)                  :: coord(1:3,1:181,1:256,1:4)
    real(4)                  :: coord_pbc(1:3,1:181,1:256,1:4)
    real(4)                  :: trans1(1:3,1:181,1:256,1:4)
    real(4)                  :: table_grad(1:7655*6,1:4)
    real(4)                  :: force(1:3,1:181,1:256,1:3,1:4)
    real(4)                  :: lj_coef(1:2,1:1000,1:4)

    common /cb_maxcell/ maxcell         ! integer
    common /cb_ncell/ ncell             ! integer
    common /cb_ncell_local/ ncell_local ! integer
    common /cb_maxatom/ MaxAtom         ! integer
    common /cb_maxatomcls/ MaxAtomCls   ! integer
    common /cb_maxnb15/ MaxNb15         ! integer
    common /cb_nthread/ nthread         ! integer
    common /cb_inv_maxatom/ inv_MaxAtom ! real(4)
    common /cb_natom/ natom             ! integer :: natom(1:256,1:4)
    common /cb_charge/ charge           ! real(4) :: charge(1:181,1:256,1:4)
    common /cb_atmcls/ atmcls           ! integer :: atmcls(1:181,1:256,1:4)
    common /cb_lj_coef/ lj_coef         ! real(4) :: lj_coef(1:2,1:1000,1:4)
    common /cb_num_nb15_calc/ num_nb15_calc   ! integer :: num_nb15_calc(1:181,1:4683,1:4)
    common /cb_nb15_calc_list/ nb15_calc_list ! integer :: nb15_calc_list(1:950,1:100,1:256,1:4)
    common /cb_coord/ coord             ! real(4) :: coord(1:3,1:181,1:256,1:4)
    common /cb_trans1/ trans1           ! real(4) :: trans1(1:3,1:181,1:256,1:4)
    common /cb_table_grad/ table_grad   ! real(4) :: table_grad(1:7655*6,1:4)
    common /cb_density/ density         ! real(4)
    common /cb_cutoff/ cutoff           ! real(4)
    common /cb_cutoff2/ cutoff2         ! real(4)

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_kernel_June"

    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_write(iunit,  maxcell     , 1)
    call sub_i4_data_write(iunit,  ncell       , 1)
    call sub_i4_data_write(iunit,  ncell_local , 1)
    call sub_i4_data_write(iunit,  MaxAtom     , 1)
    call sub_i4_data_write(iunit,  MaxAtomCls  , 1)
    call sub_i4_data_write(iunit,  MaxNb15     , 1)
    call sub_i4_data_write(iunit,  nthread     , 1)

    call sub_R4_data_write(iunit, inv_MaxAtom  , 1)

    call sub_i4_data_write(iunit, natom        , (256*4) )
    call sub_R4_data_write(iunit, charge       , (181*256*4) )
    call sub_i4_data_write(iunit, atmcls       , (181*256*4) )
    call sub_R4_data_write(iunit, lj_coef      , (2*1000*4) )
    call sub_i4_data_write(iunit, num_nb15_calc  , (181*4683*4) )
    call sub_i4_data_write(iunit, nb15_calc_list , (950*100*256*4) )
    call sub_R4_data_write(iunit, coord          , (3*181*256*4) )
    call sub_R4_data_write(iunit, trans1         , (3*181*256*4) )
    call sub_r4_data_write(iunit, table_grad   , (7655*6*4) )
    call sub_r4_data_write(iunit, density      , 1)
    call sub_r4_data_write(iunit, cutoff       , 1)
    call sub_r4_data_write(iunit, cutoff2      , 1)

    close(iunit)

end subroutine


subroutine read_data_file ()

    ! local variables
    real(4)                  :: dij1,dij2,dij3, rij2
    real(4)                  :: R, lj12, lj6
    real(4)                  :: term_lj12, term_lj6, term_elec
    real(4)                  :: cutoff2, grad_coef
    real(4)                  :: work(1:3)
    real(4)                  :: rtmp(1:3), qtmp, jqtmp
    real(4)                  :: force_local(3)
    real(4)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                  :: i, ix, iy, j, k, ij, L, L1, ii
    integer                  :: id, ik, omp_get_thread_num
    integer                  :: iatmcls,jatmcls
    integer                  :: maxcell,ncell,ncell_local
    integer                  :: nthread
    integer                  :: MaxAtom,MaxAtomCls,MaxNb15
    real(4)                  :: density, cutoff
    real(4)                  :: inv_MaxAtom
    real(8)                  :: Val

    integer                  :: natom(1:256,1:4)
    integer                  :: atmcls(1:181,1:256,1:4)
    integer                  :: num_nb15_calc(1:181,1:4683,1:4)
    integer                  :: nb15_calc_list(1:950,1:100,1:256,1:4)
    real(4)                  :: charge(1:181,1:256,1:4)
    real(4)                  :: coord(1:3,1:181,1:256,1:4)
    real(4)                  :: coord_pbc(1:3,1:181,1:256,1:4)
    real(4)                  :: trans1(1:3,1:181,1:256,1:4)
    real(4)                  :: table_grad(1:7655*6,1:4)
    real(4)                  :: force(1:3,1:181,1:256,1:3,1:4)
    real(4)                  :: lj_coef(1:2,1:1000,1:4)

    common /cb_maxcell/ maxcell         ! integer
    common /cb_ncell/ ncell             ! integer
    common /cb_ncell_local/ ncell_local ! integer
    common /cb_maxatom/ MaxAtom         ! integer
    common /cb_maxatomcls/ MaxAtomCls   ! integer
    common /cb_maxnb15/ MaxNb15         ! integer
    common /cb_nthread/ nthread         ! integer
    common /cb_inv_maxatom/ inv_MaxAtom ! real(4)
    common /cb_natom/ natom             ! integer :: natom(1:256,1:4)
    common /cb_charge/ charge           ! real(4) :: charge(1:181,1:256,1:4)
    common /cb_atmcls/ atmcls           ! integer :: atmcls(1:181,1:256,1:4)
    common /cb_lj_coef/ lj_coef         ! real(4) :: lj_coef(1:2,1:1000,1:4)
    common /cb_num_nb15_calc/ num_nb15_calc   ! integer :: num_nb15_calc(1:181,1:4683,1:4)
    common /cb_nb15_calc_list/ nb15_calc_list ! integer :: nb15_calc_list(1:950,1:100,1:256,1:4)
    common /cb_coord/ coord             ! real(4) :: coord(1:3,1:181,1:256,1:4)
    common /cb_trans1/ trans1           ! real(4) :: trans1(1:3,1:181,1:256,1:4)
    common /cb_table_grad/ table_grad   ! real(4) :: table_grad(1:7655*6,1:4)
    common /cb_density/ density         ! real(4)
    common /cb_cutoff/ cutoff           ! real(4)
    common /cb_cutoff2/ cutoff2         ! real(4)

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_kernel_June"

    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
        write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_read(iunit,  maxcell     , 1)
    call sub_i4_data_read(iunit,  ncell       , 1)
    call sub_i4_data_read(iunit,  ncell_local , 1)
    call sub_i4_data_read(iunit,  MaxAtom     , 1)
    call sub_i4_data_read(iunit,  MaxAtomCls  , 1)
    call sub_i4_data_read(iunit,  MaxNb15     , 1)
    call sub_i4_data_read(iunit,  nthread     , 1)

    call sub_r4_data_read(iunit, inv_MaxAtom  , 1)

    call sub_i4_data_read(iunit, natom        , (256*4) )
    call sub_r4_data_read(iunit, charge       , (181*256*4) )
    call sub_i4_data_read(iunit, atmcls       , (181*256*4) )
    call sub_r4_data_read(iunit, lj_coef      , (2*1000*4) )
    call sub_i4_data_read(iunit, num_nb15_calc  , (181*4683*4) )
    call sub_i4_data_read(iunit, nb15_calc_list , (950*100*256*4) )
    call sub_r4_data_read(iunit, coord          , (3*181*256*4) )
    call sub_r4_data_read(iunit, trans1         , (3*181*256*4) )
    call sub_r4_data_read(iunit, table_grad   , (7655*6*4) )
    call sub_r4_data_read(iunit, density      , 1)
    call sub_r4_data_read(iunit, cutoff       , 1)
    call sub_r4_data_read(iunit, cutoff2      , 1)

    close(iunit)

    write(*,*) "maxcell=", maxcell, "   loc(maxcell)=", loc(maxcell)

end subroutine


subroutine check_force_validation ()

    real(4) :: force(1:3,1:181,1:256,1:3,1:4)
    common /cx_force/ force
    real(8) :: val, expected_result, check_result

    val=0.0d0
    do k=1,256
      do j=1,100
        do i=1,3
          val=val+dble(force(i,j,k,1,1))
        enddo
      enddo
    enddo

    expected_result = -3285.611248243382
    check_result = val / expected_result

    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif

end subroutine

end module

